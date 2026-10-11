#pragma once

#include "solver.hpp"
#include <vector>
#include <memory>

class HybridDomain {
private:
    std::unique_ptr<Solver1D> left_domain;  // Usually FV
    std::unique_ptr<Solver1D> right_domain; // Usually DG
    
    // Interface location
    double interface_x;

public:
    enum class TimeIntegrator { ForwardEuler, SSPRK3 };

private:
    TimeIntegrator integrator = TimeIntegrator::SSPRK3;

public:
    HybridDomain(std::unique_ptr<Solver1D> left, std::unique_ptr<Solver1D> right);

    void set_time_integrator(TimeIntegrator ti) { integrator = ti; }
    TimeIntegrator get_time_integrator() const { return integrator; }

    void initialize(double t_start, double initial_dt);
    
    // Exchange boundary information between domains
    void exchange_boundaries();
    
    // Perform one time step using configured integrator
    void step(double dt, double advection_speed);

    // Explicit single-step methods
    void step_euler(double dt, double advection_speed);
    void step_ssprk3(double dt, double advection_speed);

    // Global diagnostic reductions
    double compute_total_mass() const;
    double compute_total_energy() const;

    // Equation type configuration
    void set_equation_type(EquationType eq) {
        if (left_domain) left_domain->set_equation_type(eq);
        if (right_domain) right_domain->set_equation_type(eq);
    }
    EquationType get_equation_type() const {
        return left_domain ? left_domain->get_equation_type() : EquationType::LinearAdvection;
    }

    // Modal shock sensor diagnostics
    std::vector<double> compute_modal_sensor() const;
    int get_troubled_cell_count(double threshold = -3.0) const;

    // Solver accessors
    const Solver1D* get_left_domain() const { return left_domain.get(); }
    const Solver1D* get_right_domain() const { return right_domain.get(); }

    // Get combined solution for plotting
    std::vector<std::pair<double, double>> get_solution() const;
};
