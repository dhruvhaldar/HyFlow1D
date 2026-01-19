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

    // Warmup
    for (int i = 0; i < 100; ++i) {
        fv.compute_rhs(0.0, 1.0);
    }

    int iterations = 50000;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        fv.compute_rhs(0.0, 1.0);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    double ops_per_sec = iterations / elapsed.count();
    double ns_per_op = (elapsed.count() * 1e9) / iterations;

    std::cout << "FV Benchmark (" << n_elements << " elements, " << iterations << " iterations)" << std::endl;
    std::cout << "Total Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Throughput: " << ops_per_sec << " calls/s" << std::endl;
    std::cout << "Avg Latency: " << ns_per_op << " ns/call" << std::endl;

    return 0;
}
