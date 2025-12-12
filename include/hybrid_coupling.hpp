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
    HybridDomain(std::unique_ptr<Solver1D> left, std::unique_ptr<Solver1D> right);

    void initialize(double t_start, double initial_dt);
    
    // Exchange boundary information between domains
    void exchange_boundaries();
    
    // Perform one time step
    void step(double dt, double advection_speed);

    // Get combined solution for plotting
    std::vector<std::pair<double, double>> get_solution() const;
};
