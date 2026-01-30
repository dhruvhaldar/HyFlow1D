#pragma once

#include <vector>
#include <string>
#include <memory>

// Abstract base class for a solver (or domain partition)
class Solver1D {
public:
    // Security: Global limit on elements to prevent DoS via memory exhaustion
    static constexpr int MAX_ELEMENTS = 50'000'000;

    virtual ~Solver1D() = default;

    // Initialize with domain bounds and number of elements
    virtual void initialize(double x_start, double x_end, int n_elements) = 0;

    // Set initial condition based on a function
    virtual void set_initial_condition(double (*func)(double)) = 0;

    // Compute rates of change (RHS) for RK4
    // t is current time (not used for linear advection but kept for interface)
    virtual void compute_rhs(double t, double advection_speed) = 0;

    // Update state: u = u + dt * rhs
    virtual void update_state(double dt) = 0;

    // Perform a full time step (compute_rhs + update_state)
    // Overridable for optimizations (e.g., fused update)
    virtual void step(double dt, double t, double advection_speed) {
        compute_rhs(t, advection_speed);
        update_state(dt);
    }

    // Get the solution at cell centers (for plotting)
    // Returns pairs of (x, u)
    virtual std::vector<std::pair<double, double>> get_solution() const = 0;

    // --- Coupling Interfaces ---
    
    // Get the value at the left boundary of the domain
    virtual double get_left_boundary_value() const = 0;
    
    // Get the value at the right boundary of the domain
    virtual double get_right_boundary_value() const = 0;

    // Set the ghost value (or flux info) from the left neighbor
    virtual void set_left_neighbor_value(double val) = 0;

    // Set the ghost value (or flux info) from the right neighbor
    virtual void set_right_neighbor_value(double val) = 0;
};
