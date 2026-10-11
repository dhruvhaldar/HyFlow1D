#pragma once

#include "solver.hpp"
#include <vector>

class FiniteVolumeSolver : public Solver1D {
private:
    int n_elements;
    double dx;
    double x_start, x_end;
    
    std::vector<double> u;        // Cell averages
    std::vector<double> u_stage0; // Stage buffer for SSPRK3
    std::vector<double> rhs;      // d/dt(u)
    
    // Boundary conditions / Ghost values
    double left_ghost;
    double right_ghost;

    bool is_initialized = false;

public:
    void initialize(double x_start, double x_end, int n_elements) override;
    void set_initial_condition(double (*func)(double)) override;
    void compute_rhs(double t, double advection_speed) override;
    void update_state(double dt) override;
    void save_stage() override;
    void blend_stage(double alpha, double beta, double dt) override;
    void step(double dt, double advection_speed) override;
    std::vector<std::pair<double, double>> get_solution() const override;

    double get_domain_start() const override { return x_start; }
    double get_domain_end() const override { return x_end; }
    int get_num_elements() const override { return n_elements; }

    double compute_total_mass() const override;
    double compute_total_energy() const override;

    double get_left_boundary_value() const override;
    double get_right_boundary_value() const override;
    void set_left_neighbor_value(double val) override;
    void set_right_neighbor_value(double val) override;
};
