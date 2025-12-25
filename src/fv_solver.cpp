#include "fv_solver.hpp"
#include <stdexcept>

void FiniteVolumeSolver::initialize(double start, double end, int n_elem) {
    if (n_elem <= 0) {
        throw std::invalid_argument("Number of elements must be positive.");
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
}

void FiniteVolumeSolver::set_initial_condition(double (*func)(double)) {
    for (int i = 0; i < n_elements; ++i) {
        double x_center = x_start + (i + 0.5) * dx;
        // For 1st order FV, cell average is approx value at center
        u[i] = func(x_center);
    }
}

void FiniteVolumeSolver::compute_rhs(double /*t*/, double a) {
    // 1st Order Upwind Scheme: u_i^{n+1} = u_i^n - a * dt/dx * (u_i - u_{i-1})
    // RHS = -a * (u_i - u_{i-1}) / dx
    
    // We assume a > 0 for this implementation
    
    for (int i = 0; i < n_elements; ++i) {
        double u_left = (i == 0) ? left_ghost : u[i - 1];
        // Flux at left face (i-1/2) is a * u_left
        // Flux at right face (i+1/2) is a * u[i]
        
        double flux_in = a * u_left;
        double flux_out = a * u[i];
        
        rhs[i] = -(flux_out - flux_in) / dx;
    }
}

void FiniteVolumeSolver::update_state(double dt) {
    for (size_t i = 0; i < u.size(); ++i) {
        u[i] += dt * rhs[i];
    }
}

std::vector<std::pair<double, double>> FiniteVolumeSolver::get_solution() const {
    std::vector<std::pair<double, double>> sol;
    for (int i = 0; i < n_elements; ++i) {
        sol.push_back({x_start + (i + 0.5) * dx, u[i]});
    }
    return sol;
}

double FiniteVolumeSolver::get_left_boundary_value() const {
    return u[0]; // Value in first cell
}

double FiniteVolumeSolver::get_right_boundary_value() const {
    return u.back(); // Value in last cell
}

void FiniteVolumeSolver::set_left_neighbor_value(double val) {
    left_ghost = val;
}

void FiniteVolumeSolver::set_right_neighbor_value(double val) {
    right_ghost = val;
}
