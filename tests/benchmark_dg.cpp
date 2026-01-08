#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include "dg_solver.hpp"

// Simple benchmark to measure compute_rhs performance
int main() {
    int n_elements = 200;
    int p_order = 3;
    // int n_modes = p_order + 1; // Unused
    double start = 0.0;
    double end = 1.0;

    DiscontinuousGalerkinSolver dg(p_order);
    dg.initialize(start, end, n_elements);

    // Initialize with some data
    dg.set_initial_condition([](double x) { return std::sin(10.0 * x); });

    // Warmup
    for (int i = 0; i < 100; ++i) {
        dg.compute_rhs(0.0, 1.0);
    }

    int iterations = 100000;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        dg.compute_rhs(0.0, 1.0);
        // We don't update state to avoid divergence, just benchmarking compute_rhs
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Throughput: " << iterations / elapsed.count() << " ops/s" << std::endl;

    return 0;
}
