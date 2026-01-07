#pragma once

#include "solver.hpp"
#include <vector>

class FiniteVolumeSolver : public Solver1D {
private:
    int n_elements;
    double dx;
    double x_start, x_end;
    
    std::vector<double> u;     // Cell averages
    std::vector<double> rhs;   // d/dt(u)
    
    // Boundary conditions / Ghost values
    double left_ghost;
    double right_ghost;

    bool is_initialized = false;

public:
    void initialize(double x_start, double x_end, int n_elements) override;
    void set_initial_condition(double (*func)(double)) override;
    void compute_rhs(double t, double advection_speed) override;
    void update_state(double dt) override;
    std::vector<std::pair<double, double>> get_solution() const override;

    double get_left_boundary_value() const override;
    double get_right_boundary_value() const override;
    void set_left_neighbor_value(double val) override;
    void set_right_neighbor_value(double val) override;
};
