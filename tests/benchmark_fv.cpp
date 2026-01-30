#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include "fv_solver.hpp"
#include <cmath>

// Simple benchmark to measure compute_rhs performance for FV
int main() {
    // FV is lighter per element, so use more elements to get significant timing
    int n_elements = 10000;
    double start = 0.0;
    double end = 1.0;

    FiniteVolumeSolver fv;
    fv.initialize(start, end, n_elements);

    // Initialize with some data
    fv.set_initial_condition([](double x) { return std::sin(10.0 * x); });

    int iterations = 50000;

    // --- Baseline: Split compute_rhs + update_state ---
    // Warmup
    for (int i = 0; i < 100; ++i) {
        fv.compute_rhs(0.0, 1.0);
        fv.update_state(0.001);
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        fv.compute_rhs(0.0, 1.0);
        fv.update_state(0.001);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_split = end_time - start_time;

    // --- Optimized: Fused step ---
    // Warmup
    for (int i = 0; i < 100; ++i) {
        fv.step(0.001, 0.0, 1.0);
    }

    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        fv.step(0.001, 0.0, 1.0);
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_fused = end_time - start_time;

    std::cout << "FV Benchmark (" << n_elements << " elements, " << iterations << " iterations)" << std::endl;

    std::cout << "Split (compute+update): " << elapsed_split.count() << " s ("
              << (iterations / elapsed_split.count()) << " steps/s)" << std::endl;

    std::cout << "Fused (step):           " << elapsed_fused.count() << " s ("
              << (iterations / elapsed_fused.count()) << " steps/s)" << std::endl;

    std::cout << "Speedup: " << (elapsed_split.count() / elapsed_fused.count()) << "x" << std::endl;

    return 0;
}
