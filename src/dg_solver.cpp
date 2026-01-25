#include "dg_solver.hpp"
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <limits>

#if defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

int DiscontinuousGalerkinSolver::validate_order(int p) {
    if (p < 0) {
        throw std::invalid_argument("Polynomial order must be non-negative.");
    }
    // Max supported order is 3 because n_modes = p + 1.
    // We need n_modes + 1 quadrature points.
    // Max quadrature points supported is 5.
    // p + 1 + 1 <= 5 => p + 2 <= 5 => p <= 3.
    if (p > 3) {
        throw std::invalid_argument("Polynomial order too high. Max supported order is 3 (5 quadrature points).");
    }
    return p;
}

DiscontinuousGalerkinSolver::DiscontinuousGalerkinSolver(int p_order)
    : poly_order(validate_order(p_order)), n_modes(p_order + 1) {

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

    // Precompute stiffness matrix K_km = Integral(P_m * P'_k) dxi
    // Optimization: Matrix is strictly lower triangular (0 for m >= k).
    // We pack it to store only m < k elements.
    // Size = sum(k=0 to n_modes-1) of k = n_modes * (n_modes - 1) / 2
    size_t packed_size = n_modes * (n_modes - 1) / 2;
    stiffness_matrix.resize(packed_size);

    size_t packed_idx = 0;
    // Iterate rows k
    for (int k = 0; k < n_modes; ++k) {
        // Iterate cols m < k
        for (int m = 0; m < k; ++m) {
            double sum = 0.0;
            for (size_t q = 0; q < quad_nodes.size(); ++q) {
                // sum += w_q * P_m(q) * P'_k(q)
                sum += basis_at_quad[q * n_modes + m] * weighted_d_basis_at_quad[q * n_modes + k];
            }
            stiffness_matrix[packed_idx++] = sum;
        }
    }

    is_initialized = false;
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

    // Security: Prevent Denial of Service (DoS) via excessive memory allocation.
    // Limit to 50 million elements.
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

    is_initialized = true;
}

void DiscontinuousGalerkinSolver::set_initial_condition(double (*func)(double)) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
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
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");
    if (element_idx < 0 || element_idx >= n_elements) {
        throw std::out_of_range("Element index out of bounds: " + std::to_string(element_idx));
    }
    if (std::abs(xi) > 1.0 + 1e-12) {
        throw std::domain_error("Local coordinate xi must be in [-1, 1]. Got: " + std::to_string(xi));
    }

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
        // Optimization: Use recurrence relation to evaluate all P_k(xi) in one pass O(N).
        // Avoiding repetitive calls to numerics::legendre reduces complexity from O(N^2) to O(N).

        double p_prev = 1.0; // P_0(xi)
        double p_curr = xi;  // P_1(xi)

        // Precomputed inverse integers to avoid division in recurrence
        // Supports up to k=5 (P=4, n_modes=5)
        // We add extra padding just in case.
        static constexpr double inv_k[] = {0.0, 1.0, 0.5, 0.3333333333333333, 0.25, 0.2, 0.1666666666666667, 0.14285714285714285};

        // n_modes >= 1 is guaranteed by constructor (P >= 0)
        val += u[base_idx] * p_prev;
        if (n_modes > 1) val += u[base_idx + 1] * p_curr;

        for (int k = 2; k < n_modes; ++k) {
            // Check bounds for inv_k table to be safe against future P increases
            double inv_k_val = (k < 8) ? inv_k[k] : 1.0 / static_cast<double>(k);

            // Recurrence: k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}
            // P_k = ((2k-1) x P_{k-1} - (k-1) P_{k-2}) / k
            double p_next = ((2.0 * k - 1.0) * xi * p_curr - (k - 1.0) * p_prev) * inv_k_val;

            val += u[base_idx + k] * p_next;
            p_prev = p_curr;
            p_curr = p_next;
        }
    }
    return val;
}

namespace {
    template <int N>
    void compute_rhs_optimized(int n_elements,
                               const double* RESTRICT u,
                               double* RESTRICT rhs,
                               const double* RESTRICT stiffness_matrix,
                               const double* RESTRICT inv_mass_matrix,
                               double left_ghost,
                               double a) {
        double prev_boundary_val = left_ghost;

        // Optimization: Precompute scaled inverse mass matrix
        // This avoids multiplying by 'a' inside the inner loop for every term.
        double scaled_inv_mass[N];
        for (int k = 0; k < N; ++k) {
            scaled_inv_mass[k] = inv_mass_matrix[k] * a;
        }

        // Optimization: Copy stiffness matrix to stack to ensure fast access and help compiler with aliasing.
        // Size is N*(N-1)/2. For N=4, size is 6.
        constexpr int K_size = (N * (N - 1) / 2);
        // Handle N=1 where K_size is 0.
        constexpr int local_K_size = (K_size > 0) ? K_size : 1;
        double local_K[local_K_size];
        for (int i = 0; i < K_size; ++i) local_K[i] = stiffness_matrix[i];

        for (int i = 0; i < n_elements; ++i) {
            double u_left_boundary_val = prev_boundary_val;

            // Optimization: Load u into local registers/stack to avoid repeated memory access.
            // Also compute u_right_boundary_val during the load loop.
            double u_local[N];
            double u_right_boundary_val = 0.0;
            const double* u_elem = &u[i * N];

            for (int k = 0; k < N; ++k) {
                 double val = u_elem[k];
                 u_local[k] = val;
                 u_right_boundary_val += val;
            }
            prev_boundary_val = u_right_boundary_val;

            // Optimization: Pre-calculate surface terms to avoid branches/mults in inner loop
            double val_surf_left = u_left_boundary_val;
            double val_surf_right = u_right_boundary_val;

            // For P_k(-1) = (-1)^k. Even k: +1, Odd k: -1.
            // Even k: surf_left term is val_surf_left * 1. Total = - (right - left) = -right + left
            // Odd k: surf_left term is val_surf_left * -1. Total = - (right - (-left)) = -right - left
            double surf_term_even = val_surf_right - val_surf_left;
            double surf_term_odd  = val_surf_right + val_surf_left;

            double* rhs_elem = &rhs[i * N];

            if constexpr (N <= 4) {
                // Optimization: Fully unroll volume integral calculation and RHS assembly using if constexpr.
                // This eliminates loop overhead, complex pointer arithmetic (K_ptr), and temporary arrays.

                // k=0 (Always present)
                {
                    double vol = 0.0;
                    double total_rhs = vol - surf_term_even;
                    rhs_elem[0] = total_rhs * scaled_inv_mass[0];
                }

                if constexpr (N >= 2) {
                     // k=1, m=0
                     double vol = local_K[0] * u_local[0];
                     double total_rhs = vol - surf_term_odd;
                     rhs_elem[1] = total_rhs * scaled_inv_mass[1];
                }
                if constexpr (N >= 3) {
                     // k=2, m=0,1
                     double vol = local_K[1] * u_local[0] + local_K[2] * u_local[1];
                     double total_rhs = vol - surf_term_even;
                     rhs_elem[2] = total_rhs * scaled_inv_mass[2];
                }
                if constexpr (N >= 4) {
                     // k=3, m=0,1,2
                     double vol = local_K[3] * u_local[0] + local_K[4] * u_local[1] + local_K[5] * u_local[2];
                     double total_rhs = vol - surf_term_odd;
                     rhs_elem[3] = total_rhs * scaled_inv_mass[3];
                }
            } else {
                // Generic fallback for N > 4
                const double* K_ptr = local_K;
                for (int k = 0; k < N; ++k) {
                    double volume_int = 0.0;
                    // Inner loop m < k matches the packed row length of k.
                    for (int m = 0; m < k; ++m) {
                         volume_int += K_ptr[m] * u_local[m];
                    }

                    double total_rhs = volume_int - ((k % 2 == 0) ? surf_term_even : surf_term_odd);
                    rhs_elem[k] = total_rhs * scaled_inv_mass[k];

                    // Advance pointer by row length (k) for packed storage
                    K_ptr += k;
                }
            }
        }
    }
}

void DiscontinuousGalerkinSolver::compute_rhs(double /*t*/, double a) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    // Optimization: Dispatch to fixed-size template for common polynomial orders (P=0 to P=3)
    // This allows the compiler to fully unroll the inner loops over modes (N)
    // n_modes = p_order + 1.
    switch (n_modes) {
        case 1: // P=0
            compute_rhs_optimized<1>(n_elements, u.data(), rhs.data(), stiffness_matrix.data(), inv_mass_matrix.data(), left_ghost, a);
            return;
        case 2: // P=1
            compute_rhs_optimized<2>(n_elements, u.data(), rhs.data(), stiffness_matrix.data(), inv_mass_matrix.data(), left_ghost, a);
            return;
        case 3: // P=2
            compute_rhs_optimized<3>(n_elements, u.data(), rhs.data(), stiffness_matrix.data(), inv_mass_matrix.data(), left_ghost, a);
            return;
        case 4: // P=3
            compute_rhs_optimized<4>(n_elements, u.data(), rhs.data(), stiffness_matrix.data(), inv_mass_matrix.data(), left_ghost, a);
            return;
    }

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
        const double* u_elem = &u[base_idx];

        for (int k = 0; k < n_modes; ++k) {
             u_right_boundary_val += u_elem[k];
        }
        prev_boundary_val = u_right_boundary_val;
        
        // Upwind fluxes
        double flux_surf_left = a * u_left_boundary_val;
        double flux_surf_right = a * u_right_boundary_val; 

        // Compute volume integrals for all modes using Precomputed Stiffness Matrix
        // Optimization: Replace quadrature loop with matrix-vector multiplication
        // VolInt_k = a * Integral(u * dP_k/dxi) = a * sum_m (u_m * K_km)

        const double* K_ptr = stiffness_matrix.data();
        double* rhs_elem = &rhs[base_idx];
        double sign = 1.0; // toggles 1.0, -1.0 for P_k(-1) = (-1)^k

        for (int k = 0; k < n_modes; ++k) {
            double volume_int = 0.0;
            // Matrix-vector multiplication: Row k of Stiffness * u vector
            // Optimization: Stiffness matrix is strictly lower triangular (m < k)
            // K_km = Integral(P_m * P'_k) dxi. Since deg(P'_k) = k-1 and deg(P_m) = m,
            // integral is 0 if m > k-1 (orthogonality). So only loop m < k.
            // Packed storage ensures we only store and access m < k elements.
            for (int m = 0; m < k; ++m) {
                 volume_int += K_ptr[m] * u_elem[m];
            }
            volume_int *= a;
            
            // Surface terms
            // P_k(1) = 1, P_k(-1) = (-1)^k
            double surf_right = flux_surf_right; // * 1.0
            double surf_left  = flux_surf_left * sign;
            
            double total_rhs = volume_int - (surf_right - surf_left);
            
            // Optimization: Use precomputed inverse mass matrix
            rhs_elem[k] = total_rhs * inv_mass_matrix[k];

            // Optimization: Advance pointer by row length (k) for packed storage
            K_ptr += k;
            sign = -sign;
        }
    }
}

void DiscontinuousGalerkinSolver::update_state(double dt) {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

    // Because u and rhs are now flattened vectors of the same size,
    // we can use a single loop.
    // Optimization: explicit restrict pointers to help compiler vectorization.
    double* RESTRICT u_ptr = u.data();
    const double* RESTRICT rhs_ptr = rhs.data();
    size_t total_size = u.size();

    for (size_t idx = 0; idx < total_size; ++idx) {
        u_ptr[idx] += dt * rhs_ptr[idx];
    }
}

std::vector<std::pair<double, double>> DiscontinuousGalerkinSolver::get_solution() const {
    if (!is_initialized) throw std::runtime_error("Solver not initialized. Call initialize() first.");

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
