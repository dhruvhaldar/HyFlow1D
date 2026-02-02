#include <iostream>
#include <vector>
#include <chrono>
#include "dg_solver.hpp"

// Benchmark for fused step (compute_rhs + update_state)
int main() {
    // 1,000,000 elements * 4 modes * 8 bytes = 32MB.
    // This ensures we exceed L3 cache (typically < 20MB) to stress memory bandwidth.
    int n_elements = 1000000;
    int p_order = 3;
    DiscontinuousGalerkinSolver dg(p_order);
    dg.initialize(0.0, 1.0, n_elements);

    // Set dummy data
    dg.set_initial_condition([](double x) { return x; });

    // We need to set neighbors to avoid uninitialized read if ghosts are used
    dg.set_left_neighbor_value(0.0);
    dg.set_right_neighbor_value(0.0);

    int iterations = 100;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        dg.step(0.0001, 1.0);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    size_t total_elements_updates = (size_t)n_elements * iterations;

    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Steps/s: " << iterations / elapsed.count() << std::endl;
    std::cout << "Element Updates: " << total_elements_updates / elapsed.count() / 1e6 << " M/s" << std::endl;

    return 0;
}
