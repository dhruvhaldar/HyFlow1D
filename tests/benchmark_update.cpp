#include <iostream>
#include <vector>
#include <chrono>
#include "dg_solver.hpp"

// Benchmark for update_state
int main() {
    int n_elements = 50000;
    int p_order = 3;
    DiscontinuousGalerkinSolver dg(p_order);
    dg.initialize(0.0, 1.0, n_elements);

    // Set dummy data
    dg.set_initial_condition([](double x) { return x; });

    // Populate RHS with something non-zero to force CPU work
    // Note: Since rhs and u are private, we can't write directly.
    // We can call compute_rhs once to populate it.
    dg.compute_rhs(0.0, 1.0);

    int iterations = 100000;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        dg.update_state(0.0001);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    size_t total_updates = (size_t)n_elements * (p_order + 1) * iterations;

    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Throughput: " << iterations / elapsed.count() << " calls/s" << std::endl;
    std::cout << "Element Updates: " << total_updates / elapsed.count() / 1e6 << " M/s" << std::endl;

    return 0;
}
