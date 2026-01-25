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
    // Optimization: Cache u_prev to reduce memory loads (u[i] is loaded once instead of twice)
    double u_prev = u_ptr[0];
    for (int i = 1; i < n_elements; ++i) {
        double u_curr = u_ptr[i];
        rhs_ptr[i] = coeff * (u_curr - u_prev);
        u_prev = u_curr;
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
