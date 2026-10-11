#include "hybrid_coupling.hpp"
#include <stdexcept>
#include <cmath>

HybridDomain::HybridDomain(std::unique_ptr<Solver1D> left, std::unique_ptr<Solver1D> right)
    : left_domain(std::move(left)), right_domain(std::move(right)) {
    if (!left_domain) {
        throw std::invalid_argument("Left domain solver cannot be null.");
    }
    if (!right_domain) {
        throw std::invalid_argument("Right domain solver cannot be null.");
    }

    // Security Enhancement: Validate Domain Connectivity
    // Ensure the interface is continuous (within tolerance) to prevent
    // physically invalid configurations (teleportation/gaps).
    double left_end = left_domain->get_domain_end();
    double right_start = right_domain->get_domain_start();

    if (std::abs(left_end - right_start) > 1e-9) {
        throw std::invalid_argument("Domain mismatch: Gap or Overlap detected at interface ("
                                    + std::to_string(left_end) + " vs "
                                    + std::to_string(right_start) + ").");
    }
}

void HybridDomain::exchange_boundaries() {
    // Left domain's right ghost value comes from Right domain's left boundary
    double right_val = right_domain->get_left_boundary_value();
    left_domain->set_right_neighbor_value(right_val);

    // Right domain's left ghost value comes from Left domain's right boundary
    double left_val = left_domain->get_right_boundary_value();
    right_domain->set_left_neighbor_value(left_val);
}

void HybridDomain::step_euler(double dt, double advection_speed) {
    exchange_boundaries();
    left_domain->compute_rhs(0.0, advection_speed);
    right_domain->compute_rhs(0.0, advection_speed);
    left_domain->update_state(dt);
    right_domain->update_state(dt);
}

void HybridDomain::step_ssprk3(double dt, double advection_speed) {
    // Stage 1: u(1) = u(0) + dt * L(u(0))
    exchange_boundaries();
    left_domain->compute_rhs(0.0, advection_speed);
    right_domain->compute_rhs(0.0, advection_speed);

    left_domain->save_stage();
    right_domain->save_stage();

    left_domain->blend_stage(0.0, 1.0, dt);
    right_domain->blend_stage(0.0, 1.0, dt);

    // Stage 2: u(2) = 3/4 u(0) + 1/4 (u(1) + dt * L(u(1)))
    exchange_boundaries();
    left_domain->compute_rhs(0.0, advection_speed);
    right_domain->compute_rhs(0.0, advection_speed);

    left_domain->blend_stage(0.75, 0.25, dt);
    right_domain->blend_stage(0.75, 0.25, dt);

    // Stage 3: u(n+1) = 1/3 u(0) + 2/3 (u(2) + dt * L(u(2)))
    exchange_boundaries();
    left_domain->compute_rhs(0.0, advection_speed);
    right_domain->compute_rhs(0.0, advection_speed);

    left_domain->blend_stage(1.0 / 3.0, 2.0 / 3.0, dt);
    right_domain->blend_stage(1.0 / 3.0, 2.0 / 3.0, dt);
}

void HybridDomain::step(double dt, double advection_speed) {
    if (integrator == TimeIntegrator::SSPRK3) {
        step_ssprk3(dt, advection_speed);
    } else {
        step_euler(dt, advection_speed);
    }
}

double HybridDomain::compute_total_mass() const {
    return left_domain->compute_total_mass() + right_domain->compute_total_mass();
}

double HybridDomain::compute_total_energy() const {
    return left_domain->compute_total_energy() + right_domain->compute_total_energy();
}

std::vector<std::pair<double, double>> HybridDomain::get_solution() const {
    auto sol_left = left_domain->get_solution();
    auto sol_right = right_domain->get_solution();
    
    sol_left.insert(sol_left.end(), sol_right.begin(), sol_right.end());
    return sol_left;
}
