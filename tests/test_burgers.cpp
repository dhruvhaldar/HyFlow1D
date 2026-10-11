#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"

// Smooth initial condition
double smooth_sine(double x) {
    return 0.5 + 0.5 * std::sin(2.0 * M_PI * x);
}

// Test FV solver with Burgers equation and Rusanov flux
void test_fv_burgers() {
    std::cout << "Testing FiniteVolumeSolver on Burgers' equation..." << std::endl;
    FiniteVolumeSolver fv;
    fv.initialize(0.0, 1.0, 100);
    fv.set_equation_type(EquationType::Burgers);
    fv.set_initial_condition(smooth_sine);

    double initial_mass = fv.compute_total_mass();
    assert(std::abs(initial_mass - 0.5) < 1e-2);

    double dt = 0.001;
    // Step forward 50 steps
    for (int step = 0; step < 50; ++step) {
        fv.set_left_neighbor_value(fv.get_right_boundary_value());
        fv.set_right_neighbor_value(fv.get_left_boundary_value());
        fv.step(dt, 1.0);
    }

    double final_mass = fv.compute_total_mass();
    std::cout << "  FV Burgers mass: initial=" << initial_mass << ", final=" << final_mass << std::endl;
    // Rusanov flux is strictly conservative
    assert(std::abs(final_mass - initial_mass) < 1e-4);
    std::cout << "  Passed." << std::endl;
}

// Test DG solver with Burgers equation, Gauss-Legendre quadrature, and Rusanov flux
void test_dg_burgers() {
    std::cout << "Testing DiscontinuousGalerkinSolver on Burgers' equation..." << std::endl;
    DiscontinuousGalerkinSolver dg(2); // P=2
    dg.initialize(0.0, 1.0, 20);
    dg.set_equation_type(EquationType::Burgers);
    dg.set_initial_condition(smooth_sine);

    double initial_mass = dg.compute_total_mass();
    std::cout << "  DG initial mass: " << initial_mass << std::endl;
    assert(std::abs(initial_mass - 0.5) < 1e-2);

    double dt = 0.0005;
    for (int step = 0; step < 50; ++step) {
        dg.set_left_neighbor_value(dg.get_right_boundary_value());
        dg.set_right_neighbor_value(dg.get_left_boundary_value());
        dg.compute_rhs(0.0, 1.0);
        dg.update_state(dt);
    }

    double final_mass = dg.compute_total_mass();
    std::cout << "  DG Burgers mass: initial=" << initial_mass << ", final=" << final_mass << std::endl;
    assert(std::abs(final_mass - initial_mass) < 1e-4);
    std::cout << "  Passed." << std::endl;
}

// Test HybridDomain coupling with Burgers equation and SSPRK3
void test_hybrid_burgers_ssprk3() {
    std::cout << "Testing HybridDomain coupling on Burgers' equation with SSPRK3..." << std::endl;
    auto fv = std::make_unique<FiniteVolumeSolver>();
    fv->initialize(0.0, 0.5, 50);
    fv->set_equation_type(EquationType::Burgers);
    fv->set_initial_condition(smooth_sine);

    auto dg = std::make_unique<DiscontinuousGalerkinSolver>(2);
    dg->initialize(0.5, 1.0, 20);
    dg->set_equation_type(EquationType::Burgers);
    dg->set_initial_condition(smooth_sine);

    HybridDomain hybrid(std::move(fv), std::move(dg));
    hybrid.set_equation_type(EquationType::Burgers);
    hybrid.set_time_integrator(HybridDomain::TimeIntegrator::SSPRK3);

    double initial_mass = hybrid.compute_total_mass();
    double dt = 0.0005;

    for (int step = 0; step < 50; ++step) {
        hybrid.step(dt, 1.0);
    }

    double final_mass = hybrid.compute_total_mass();
    std::cout << "  Hybrid Burgers SSPRK3 mass: initial=" << initial_mass << ", final=" << final_mass << std::endl;
    // With inflow/outflow, mass is bounded and finite
    assert(!std::isnan(final_mass));
    assert(!std::isinf(final_mass));
    std::cout << "  Passed." << std::endl;
}

int main() {
    try {
        test_fv_burgers();
        test_dg_burgers();
        test_hybrid_burgers_ssprk3();
        std::cout << "All Burgers equation tests passed successfully!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}
