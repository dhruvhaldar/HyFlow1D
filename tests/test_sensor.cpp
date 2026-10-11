#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include "dg_solver.hpp"
#include "fv_solver.hpp"
#include "hybrid_coupling.hpp"

// Smooth profile
double smooth_wave(double x) {
    return std::sin(2.0 * M_PI * x);
}

// Discontinuous profile (Heaviside step located inside an element)
double step_wave(double x) {
    return (x >= 0.52) ? 1.0 : 0.0;
}

void test_smooth_decay() {
    std::cout << "Testing Persson-Peraire sensor on smooth solution..." << std::endl;
    DiscontinuousGalerkinSolver dg(3); // P=3
    dg.initialize(0.0, 1.0, 16);
    dg.set_initial_condition(smooth_wave);

    auto sensor = dg.compute_modal_sensor();
    assert(sensor.size() == 16);

    int troubled = dg.get_troubled_cell_count(-3.0);
    std::cout << "  Smooth solution troubled elements (threshold=-3.0): " << troubled << std::endl;
    for (size_t i = 0; i < sensor.size(); ++i) {
        std::cout << "    Elem " << i << ": S_e = " << sensor[i] << std::endl;
        // On smooth wave, modal energy should decay significantly
        assert(sensor[i] < -3.0);
    }
    assert(troubled == 0);
    std::cout << "  Passed." << std::endl;
}

void test_discontinuous_detection() {
    std::cout << "Testing Persson-Peraire sensor on discontinuous step..." << std::endl;
    DiscontinuousGalerkinSolver dg(3); // P=3
    dg.initialize(0.0, 1.0, 16);
    dg.set_initial_condition(step_wave);

    auto sensor = dg.compute_modal_sensor();
    assert(sensor.size() == 16);

    int troubled = dg.get_troubled_cell_count(-3.0);
    std::cout << "  Discontinuous solution troubled elements: " << troubled << std::endl;
    // Discontinuity must trigger at least one troubled element
    assert(troubled >= 1);

    // Specifically element near 0.5 should be flagged
    // Element around index 8 spans [0.5, 0.5625]
    bool found_shock_cell = false;
    for (size_t i = 0; i < sensor.size(); ++i) {
        if (sensor[i] > -3.0) {
            found_shock_cell = true;
            std::cout << "    Detected troubled cell at elem " << i << " (S_e = " << sensor[i] << ")" << std::endl;
        }
    }
    assert(found_shock_cell);
    std::cout << "  Passed." << std::endl;
}

void test_hybrid_sensor() {
    std::cout << "Testing HybridDomain modal sensor aggregation..." << std::endl;
    auto fv = std::make_unique<FiniteVolumeSolver>();
    fv->initialize(0.0, 0.5, 10);
    fv->set_initial_condition(step_wave);

    auto dg = std::make_unique<DiscontinuousGalerkinSolver>(3);
    dg->initialize(0.5, 1.0, 10);
    dg->set_initial_condition(step_wave);

    HybridDomain hybrid(std::move(fv), std::move(dg));
    auto sensor = hybrid.compute_modal_sensor();
    assert(sensor.size() == 20); // 10 FV + 10 DG

    int troubled = hybrid.get_troubled_cell_count(-3.0);
    std::cout << "  Hybrid troubled elements: " << troubled << std::endl;
    assert(troubled >= 1);
    std::cout << "  Passed." << std::endl;
}

int main() {
    try {
        test_smooth_decay();
        test_discontinuous_detection();
        test_hybrid_sensor();
        std::cout << "All modal shock sensor tests passed successfully!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Sensor test failed: " << e.what() << std::endl;
        return 1;
    }
}
