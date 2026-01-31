#include "fv_solver.hpp"
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
    // 1st Order Upwind Scheme: u_i^{n+1} = u_i^n - a * dt/dx * (u_i - u_{i-1})
    // RHS = -a * (u_i - u_{i-1}) / dx
    
    // We assume a > 0 for this implementation
    
    // Optimization: Precompute coefficient to avoid repeated divisions
    // RHS = -a * (u_i - u_{i-1}) / dx = (-a/dx) * (u_i - u_{i-1})
    const double coeff = -a / dx;

    double* RESTRICT rhs_ptr = rhs.data();
    const double* RESTRICT u_ptr = u.data();

    // Peel first iteration to remove branch from main loop
    if (n_elements > 0) {
        rhs_ptr[0] = coeff * (u_ptr[0] - left_ghost);
    }

    // Main loop - no branching, vectorization friendly
    // Optimization: Loop structure allows auto-vectorization by compiler.
    // Modern compilers handle the adjacent loads (u[i] and u[i-1]) efficiently.
    for (int i = 1; i < n_elements; ++i) {
        rhs_ptr[i] = coeff * (u_ptr[i] - u_ptr[i-1]);
    }
}

void FiniteVolumeSolver::step(double dt, double /*t*/, double a) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    // Fused compute_rhs and update_state for performance.
    // Avoids writing to intermediate rhs vector and reading it back.
    // Standard 1st Order Upwind: u_i += dt * (-a/dx * (u_i - u_{i-1}))

    const double coeff = -a / dx;
    double* RESTRICT u_ptr = u.data();

    if (n_elements > 0) {
        // Handle first element (depends on ghost)
        double u_prev_val = u_ptr[0]; // Save old u[0]
        double rhs_0 = coeff * (u_prev_val - left_ghost);
        u_ptr[0] += dt * rhs_0;

        // Loop uses prev_u as the "old u[i-1]"
        double prev_u = u_prev_val;

        int i = 1;
        // 4-way unrolling to break dependency chain and enable vectorization
        for (; i <= n_elements - 4; i += 4) {
            double u0 = u_ptr[i];
            double u1 = u_ptr[i+1];
            double u2 = u_ptr[i+2];
            double u3 = u_ptr[i+3];

            u_ptr[i]   += dt * coeff * (u0 - prev_u);
            u_ptr[i+1] += dt * coeff * (u1 - u0);
            u_ptr[i+2] += dt * coeff * (u2 - u1);
            u_ptr[i+3] += dt * coeff * (u3 - u2);

            prev_u = u3;
        }

        // Tail cleanup
        for (; i < n_elements; ++i) {
            double curr_u = u_ptr[i];
            double rhs = coeff * (curr_u - prev_u);
            u_ptr[i] += dt * rhs;
            prev_u = curr_u;
        }
    }
}

void FiniteVolumeSolver::update_state(double dt) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    double* RESTRICT u_ptr = u.data();
    const double* RESTRICT rhs_ptr = rhs.data();
    size_t size = u.size();

    for (size_t i = 0; i < size; ++i) {
        u_ptr[i] += dt * rhs_ptr[i];
    }
}

std::vector<std::pair<double, double>> FiniteVolumeSolver::get_solution() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    std::vector<std::pair<double, double>> sol;
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
