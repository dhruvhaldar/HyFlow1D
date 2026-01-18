#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include <cmath>
#include "fv_solver.hpp"

// Simple benchmark to measure compute_rhs performance
int main() {
    int n_elements = 10000;
    double start = 0.0;
    double end = 1.0;

    FiniteVolumeSolver fv;
    fv.initialize(start, end, n_elements);

    // Initialize with some data
    fv.set_initial_condition([](double x) -> double { return std::sin(10.0 * x); });

    // Warmup
    for (int i = 0; i < 100; ++i) {
        fv.compute_rhs(0.0, 1.0);
    }

    int iterations = 100000;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        fv.compute_rhs(0.0, 1.0);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Throughput: " << iterations / elapsed.count() << " ops/s" << std::endl;

    return 0;
}
