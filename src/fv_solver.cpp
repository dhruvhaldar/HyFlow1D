#include "fv_solver.hpp"
#include "numerics.hpp"
#include <stdexcept>
#include <string>

#if defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

void FiniteVolumeSolver::initialize(double start, double end, int n_elem) {
    if (n_elem <= 0) {
        throw std::invalid_argument("Number of elements must be positive.");
    }

    // Security: Prevent Denial of Service (DoS) via excessive memory allocation.
    // Limit to 50 million elements (approx 400MB per vector).
    if (n_elem > MAX_ELEMENTS) {
         throw std::length_error("Number of elements exceeds security limit: " + std::to_string(MAX_ELEMENTS));
    }

    if (start >= end) {
        throw std::invalid_argument("Domain start must be less than end.");
    }

    x_start = start;
    x_end = end;
    n_elements = n_elem;
    dx = (x_end - x_start) / n_elements;
    u.resize(n_elements);
    u_stage0.resize(n_elements);
    rhs.resize(n_elements);
    left_ghost = 0.0;
    right_ghost = 0.0;
    is_initialized = true;
}

void FiniteVolumeSolver::set_initial_condition(double (*func)(double)) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    for (int i = 0; i < n_elements; ++i) {
        double x_center = x_start + (i + 0.5) * dx;
        // For 1st order FV, cell average is approx value at center
        u[i] = func(x_center);
    }
}

void FiniteVolumeSolver::compute_rhs(double /*t*/, double a) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    double* RESTRICT rhs_ptr = rhs.data();
    const double* RESTRICT u_ptr = u.data();
    const double inv_dx = 1.0 / dx;

    if (equation_type == EquationType::LinearAdvection && a >= 0.0) {
        // 1st Order Upwind Scheme: u_i^{n+1} = u_i^n - a * dt/dx * (u_i - u_{i-1})
        // RHS = -a * (u_i - u_{i-1}) / dx
        const double coeff = -a * inv_dx;

        // Peel first iteration to remove branch from main loop
        if (n_elements > 0) {
            rhs_ptr[0] = coeff * (u_ptr[0] - left_ghost);
        }

        // Main loop: parallelized with OpenMP across elements
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < n_elements; ++i) {
            rhs_ptr[i] = coeff * (u_ptr[i] - u_ptr[i-1]);
        }
    } else {
        // Conservative Finite Volume formulation with Rusanov flux
        // d(u_i)/dt = - (F_{i+1/2} - F_{i-1/2}) / dx
        std::vector<double> num_flux(n_elements + 1);
        num_flux[0] = numerics::rusanov_flux(left_ghost, u_ptr[0], equation_type, a);
        num_flux[n_elements] = numerics::rusanov_flux(u_ptr[n_elements - 1], right_ghost, equation_type, a);

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < n_elements; ++i) {
            num_flux[i] = numerics::rusanov_flux(u_ptr[i-1], u_ptr[i], equation_type, a);
        }

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_elements; ++i) {
            rhs_ptr[i] = -inv_dx * (num_flux[i+1] - num_flux[i]);
        }
    }
}

void FiniteVolumeSolver::update_state(double dt) {
    blend_stage(0.0, 1.0, dt);
}

void FiniteVolumeSolver::save_stage() {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    u_stage0 = u;
}

void FiniteVolumeSolver::blend_stage(double alpha, double beta, double dt) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    double* RESTRICT u_ptr = u.data();
    const double* RESTRICT u0_ptr = u_stage0.data();
    const double* RESTRICT rhs_ptr = rhs.data();
    const size_t size = u.size();

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < size; ++i) {
        u_ptr[i] = alpha * u0_ptr[i] + beta * (u_ptr[i] + dt * rhs_ptr[i]);
    }
}

void FiniteVolumeSolver::step(double dt, double advection_speed) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    if (equation_type == EquationType::LinearAdvection && advection_speed >= 0.0) {
        // Optimization: Fused Step
        // Combines compute_rhs and update_state into a single pass to reduce memory bandwidth.
        double* RESTRICT u_ptr = u.data();
        const double coeff = advection_speed * dt / dx;
        double prev_u = left_ghost;

        for (int i = 0; i < n_elements; ++i) {
            double current_u = u_ptr[i];
            u_ptr[i] = current_u - coeff * (current_u - prev_u);
            prev_u = current_u;
        }
        return;
    }

    // Generic Fallback
    compute_rhs(0.0, advection_speed);
    update_state(dt);
}

std::vector<std::pair<double, double>> FiniteVolumeSolver::get_solution() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    std::vector<std::pair<double, double>> sol;
    sol.reserve(n_elements);
    for (int i = 0; i < n_elements; ++i) {
        sol.push_back({x_start + (i + 0.5) * dx, u[i]});
    }
    return sol;
}

double FiniteVolumeSolver::get_left_boundary_value() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    return u[0]; // Value in first cell
}

double FiniteVolumeSolver::get_right_boundary_value() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    return u.back(); // Value in last cell
}

void FiniteVolumeSolver::set_left_neighbor_value(double val) {
    left_ghost = val;
}

void FiniteVolumeSolver::set_right_neighbor_value(double val) {
    right_ghost = val;
}

double FiniteVolumeSolver::compute_total_mass() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    double total = 0.0;
    const double* RESTRICT u_ptr = u.data();
    const int n = n_elements;

    #pragma omp parallel for reduction(+:total) schedule(static)
    for (int i = 0; i < n; ++i) {
        total += u_ptr[i] * dx;
    }
    return total;
}

double FiniteVolumeSolver::compute_total_energy() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    double total = 0.0;
    const double* RESTRICT u_ptr = u.data();
    const int n = n_elements;

    #pragma omp parallel for reduction(+:total) schedule(static)
    for (int i = 0; i < n; ++i) {
        total += 0.5 * u_ptr[i] * u_ptr[i] * dx;
    }
    return total;
}
