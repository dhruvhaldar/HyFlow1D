#include "dg_solver.hpp"
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <limits>

DiscontinuousGalerkinSolver::DiscontinuousGalerkinSolver(int p_order) 
    : poly_order(p_order), n_modes(p_order + 1) {
    if (p_order < 0) {
        throw std::invalid_argument("Polynomial order must be non-negative.");
    }

    // Safety check: numerics::gauss_legendre supports max 5 points
    // We need N = p_order + 2 points (p_order+1 modes, plus one extra for safety in quadrature?)
    // Actually the code uses n_modes + 1 = p_order + 2.
    // If p_order = 3, n_modes = 4, we request 5. Supported.
    // If p_order = 4, n_modes = 5, we request 6. Unsupported.
    // Max supported n_modes + 1 is 5 => Max n_modes is 4 => Max p_order is 3.
    if (n_modes + 1 > 5) {
        throw std::invalid_argument("Polynomial order too high. Max supported order is 3 (5 quadrature points).");
    }

    // Get quadrature for accurate integration of mass matrix and stiffness
    // We need to integrate basis*basis, which is order 2P. 
    // Gauss-Legendre with N points integrates 2N-1 exactly.
    // 2N-1 >= 2P -> 2N >= 2P+1 -> N >= P + 0.5. So N=P+1 is sufficient.
    auto qw = numerics::gauss_legendre(n_modes + 1); 
    quad_nodes = qw.first;
    quad_weights = qw.second;

    // Precompute basis
    basis_at_quad.resize(quad_nodes.size() * n_modes);
    d_basis_at_quad.resize(quad_nodes.size() * n_modes);
    weighted_d_basis_at_quad.resize(quad_nodes.size() * n_modes);

    for(size_t q=0; q<quad_nodes.size(); ++q) {
        double w = quad_weights[q];
        for(int k=0; k<n_modes; ++k) {
             size_t idx = q * n_modes + k;
             basis_at_quad[idx] = numerics::legendre(k, quad_nodes[q]);
             d_basis_at_quad[idx] = numerics::legendre_derivative(k, quad_nodes[q]);
             // Optimization: Pre-multiply weight into derivative basis
             weighted_d_basis_at_quad[idx] = w * d_basis_at_quad[idx];
        }
    }

    // Initialize scratch space
    volume_ints_scratch.resize(n_modes);
}

void DiscontinuousGalerkinSolver::initialize(double start, double end, int n_elem) {
    if (n_elem <= 0) {
        throw std::invalid_argument("Number of elements must be positive.");
    }
    // Check for integer overflow in state allocation (n_elem * n_modes)
    // If n_elem * n_modes > MAX_INT, iterating with int indices will wrap and cause OOB access.
    if (n_elem > std::numeric_limits<int>::max() / n_modes) {
        throw std::overflow_error("Number of elements too large, would cause integer overflow.");
    }

    if (start >= end) {
        throw std::invalid_argument("Domain start must be less than end.");
    }

    x_start = start;
    x_end = end;
    n_elements = n_elem;
    dx = (x_end - x_start) / n_elements;
    
    // Flattened storage
    u.assign(n_elements * n_modes, 0.0);
    rhs.assign(n_elements * n_modes, 0.0);
    
    // Precompute inverse mass matrix diagonal
    // M_kk = dx / (2k+1)
    // inv_M_kk = (2k+1) / dx
    inv_mass_matrix.resize(n_modes);
    for (int k = 0; k < n_modes; ++k) {
        inv_mass_matrix[k] = (2.0 * k + 1.0) / dx;
    }

    left_ghost = 0.0;
    right_ghost = 0.0;
}

void DiscontinuousGalerkinSolver::set_initial_condition(double (*func)(double)) {
    // Project initial condition onto basis
    // u_i^k = (2k+1)/2 * integral(f(x)*P_k(xi)) dx
    // Coordinate transform: x = x_center + xi * dx/2
    
    for (int i = 0; i < n_elements; ++i) {
        double x_center = x_start + (i + 0.5) * dx;
        
        for (int k = 0; k < n_modes; ++k) {
            double integral = 0.0;
            for (size_t q = 0; q < quad_nodes.size(); ++q) {
                double w = quad_weights[q];
                double x_phys = x_center + quad_nodes[q] * dx / 2.0;
                
                size_t basis_idx = q * n_modes + k;
                integral += w * func(x_phys) * basis_at_quad[basis_idx];
            }
            // Mass matrix diagonal term is 2/(2k+1) * (dx/2)
            // But we are working in standard element [-1, 1], so factor is just (2k+1)/2
            // Wait, careful.
            // Standard DG: M u_dot = RHS.
            // Basis P_k. Integral P_k P_m = 2/(2k+1) delta_km.
            // Physical space integral introduces Jacobian J = dx/2.
            // So mass matrix M_kk = (dx/2) * (2/(2k+1)).
            // RHS vector entry R_k = Integral(f * P_k) dx = (dx/2) * Integral(f(xi) * P_k(xi)) dxi.
            // u_k = R_k / M_kk = Integral(f(xi) * P_k(xi)) dxi * (2k+1)/2.
            
            u[i * n_modes + k] = integral * (2.0 * k + 1.0) / 2.0;
        }
    }
}

double DiscontinuousGalerkinSolver::evaluate_element(int element_idx, double xi) const {
    double val = 0.0;
    int base_idx = element_idx * n_modes;
    if (xi == 1.0) {
        for (int k = 0; k < n_modes; ++k) {
            val += u[base_idx + k]; // P_k(1) = 1
        }
    } else if (xi == -1.0) {
        for (int k = 0; k < n_modes; ++k) {
            val += u[base_idx + k] * ((k % 2 == 0) ? 1.0 : -1.0); // P_k(-1) = (-1)^k
        }
    } else {
        for (int k = 0; k < n_modes; ++k) {
            val += u[base_idx + k] * numerics::legendre(k, xi);
        }
    }
    return val;
}

void DiscontinuousGalerkinSolver::compute_rhs(double /*t*/, double a) {
    // DG Formulation:
    // (dx/2) * (2/(2k+1)) * du_k/dt = Volume_Integral - Surface_Terms
    // Volume_Integral = Integral(a * u * dP_k/dx) dx = a * Integral(u * dP_k/dxi * 2/dx) * dx/2 dxi
    //                 = a * Integral(u * dP_k/dxi) dxi
    // Surface_Terms   = [a * u * P_k] boundaries
    //                 = Flux_Right * P_k(1) - Flux_Left * P_k(-1)
    
    // Assume a > 0 (Upwind flux)
    // Flux at boundary is a * u_upwind.
    
    double prev_boundary_val = left_ghost;

    for (int i = 0; i < n_elements; ++i) {
        double u_left_boundary_val = prev_boundary_val;

        // Optimization: evaluate_element(i, 1.0) is effectively summing u[base_idx + k].
        // Implementing this inline avoids the function call and the generic loop in evaluate_element.
        double u_right_boundary_val = 0.0;
        int base_idx = i * n_modes;
        for (int k = 0; k < n_modes; ++k) {
             u_right_boundary_val += u[base_idx + k];
        }
        prev_boundary_val = u_right_boundary_val;
        
        // Upwind fluxes
        double flux_surf_left = a * u_left_boundary_val;
        double flux_surf_right = a * u_right_boundary_val; 

        // Compute volume integrals for all modes
        // Loop Fusion Optimization:
        // We calculate u(xi_q) and immediately use it to update volume integrals.
        // This eliminates the need for the intermediate 'u_at_quad_scratch' vector,
        // reducing memory traffic and iterating over quad nodes only once.

        std::fill(volume_ints_scratch.begin(), volume_ints_scratch.end(), 0.0);

        for (size_t q = 0; q < quad_nodes.size(); ++q) {
            size_t basis_row_start = q * n_modes;

            // 1. Compute u at quad node q on the fly
            double u_node_val = 0.0;
            for (int k = 0; k < n_modes; ++k) {
                 u_node_val += u[base_idx + k] * basis_at_quad[basis_row_start + k];
            }

            // 2. Compute flux/term contribution
            double flux_val = u_node_val * a;

            // 3. Accumulate to volume integrals
            // Access weighted_d_basis_at_quad sequentially
            for (int k = 0; k < n_modes; ++k) {
                 volume_ints_scratch[k] += flux_val * weighted_d_basis_at_quad[basis_row_start + k];
            }
        }

        for (int k = 0; k < n_modes; ++k) {
            double volume_int = volume_ints_scratch[k];
            
            // Surface terms
            // P_k(1) = 1, P_k(-1) = (-1)^k
            double surf_right = flux_surf_right * 1.0; 
            double surf_left  = flux_surf_left * ((k % 2 == 0) ? 1.0 : -1.0);
            
            double total_rhs = volume_int - (surf_right - surf_left);
            
            // Optimization: Use precomputed inverse mass matrix
            rhs[base_idx + k] = total_rhs * inv_mass_matrix[k];
        }
    }
}

void DiscontinuousGalerkinSolver::update_state(double dt) {
    // Because u and rhs are now flattened vectors of the same size,
    // we can use a single loop.
    // This allows the compiler to vectorize efficiently.
    size_t total_size = u.size();
    for (size_t idx = 0; idx < total_size; ++idx) {
        u[idx] += dt * rhs[idx];
    }
}

std::vector<std::pair<double, double>> DiscontinuousGalerkinSolver::get_solution() const {
    std::vector<std::pair<double, double>> sol;
    int points_per_elem = 5; // Resolution for plotting
    for (int i = 0; i < n_elements; ++i) {
        double x_center = x_start + (i + 0.5) * dx;
        for (int j = 0; j <= points_per_elem; ++j) {
            double xi = -1.0 + 2.0 * j / points_per_elem;
            double val = evaluate_element(i, xi);
            double x_phys = x_center + xi * dx / 2.0;
            sol.push_back({x_phys, val});
        }
    }
    return sol;
}

double DiscontinuousGalerkinSolver::get_left_boundary_value() const {
    return evaluate_element(0, -1.0);
}

double DiscontinuousGalerkinSolver::get_right_boundary_value() const {
    return evaluate_element(n_elements - 1, 1.0);
}

void DiscontinuousGalerkinSolver::set_left_neighbor_value(double val) {
    left_ghost = val;
}

void DiscontinuousGalerkinSolver::set_right_neighbor_value(double val) {
    right_ghost = val;
}
