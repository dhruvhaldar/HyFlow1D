#pragma once

#include "solver.hpp"
#include "numerics.hpp"
#include <vector>

class DiscontinuousGalerkinSolver : public Solver1D {
private:
    int n_elements;
    int poly_order; // P
    int n_modes;    // P + 1
    double dx;
    double x_start, x_end;
    
    // Storage: u[element_idx][mode_idx]
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> rhs;

    // Precomputed basis values at quadrature points
    // basis_at_quad[mode][quad_point]
    std::vector<std::vector<double>> basis_at_quad;
    std::vector<std::vector<double>> d_basis_at_quad;

    // Quadrature rules for integration
    std::vector<double> quad_nodes;
    std::vector<double> quad_weights;

    // Ghost values for numerical flux
    double left_ghost;
    double right_ghost;

public:
    DiscontinuousGalerkinSolver(int p_order);

    void initialize(double x_start, double x_end, int n_elements) override;
    void set_initial_condition(double (*func)(double)) override;
    void compute_rhs(double t, double advection_speed) override;
    void update_state(double dt) override;
    std::vector<std::pair<double, double>> get_solution() const override;

    double get_left_boundary_value() const override;
    double get_right_boundary_value() const override;
    void set_left_neighbor_value(double val) override;
    void set_right_neighbor_value(double val) override;

private:
    // Helper to evaluate solution at local coordinate xi in [-1, 1] for element i
    double evaluate_element(int element_idx, double xi) const;
};
