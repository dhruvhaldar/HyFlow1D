#include <iostream>
#include <cmath>
#include <memory>
#include <cassert>
#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"

double gaussian_pulse(double x) {
    double center = 0.25;
    double width = 0.05;
    return std::exp(-(x - center) * (x - center) / (2.0 * width * width));
}

int main() {
    std::cout << "Running test_ssprk3..." << std::endl;

    double x_start = 0.0;
    double x_interface = 0.5;
    double x_end = 1.0;
    int n_fv = 50;
    int n_dg = 10;
    int p_order = 3;

    auto fv = std::make_unique<FiniteVolumeSolver>();
    fv->initialize(x_start, x_interface, n_fv);
    fv->set_initial_condition(gaussian_pulse);

    auto dg = std::make_unique<DiscontinuousGalerkinSolver>(p_order);
    dg->initialize(x_interface, x_end, n_dg);
    dg->set_initial_condition(gaussian_pulse);

    HybridDomain hybrid(std::move(fv), std::move(dg));
    hybrid.set_time_integrator(HybridDomain::TimeIntegrator::SSPRK3);

    double initial_mass = hybrid.compute_total_mass();
    double initial_energy = hybrid.compute_total_energy();

    assert(initial_mass > 0.0);
    assert(initial_energy > 0.0);

    double advection_speed = 1.0;
    double dt = 0.0005;
    double t_final = 0.2; // pulse stays inside domain
    double t = 0.0;

    while (t < t_final) {
        hybrid.step(dt, advection_speed);
        t += dt;
    }

    double final_mass = hybrid.compute_total_mass();
    double final_energy = hybrid.compute_total_energy();

    // In linear advection with zero flux at boundaries before pulse exits, mass should be conserved
    double mass_error = std::abs(final_mass - initial_mass) / initial_mass;
    std::cout << "SSPRK3 Mass conservation relative error: " << mass_error << std::endl;
    assert(mass_error < 1e-4);

    // Energy should be stable (non-increasing due to upwind dissipation)
    std::cout << "Initial Energy: " << initial_energy << ", Final Energy: " << final_energy << std::endl;
    assert(final_energy <= initial_energy + 1e-12);

    std::cout << "test_ssprk3 passed successfully." << std::endl;
    return 0;
}
