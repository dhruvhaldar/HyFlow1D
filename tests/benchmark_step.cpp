#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include "fv_solver.hpp"
#include <cmath>

// Benchmark to measure full step performance (compute_rhs + update_state) for FV
int main() {
    int n_elements = 10000;
    double start = 0.0;
    double end = 1.0;

    FiniteVolumeSolver fv;
    fv.initialize(start, end, n_elements);
    fv.set_initial_condition([](double x) { return std::sin(10.0 * x); });

    double dt = 0.0001;
    double advection_speed = 1.0;

    // Warmup
    for (int i = 0; i < 100; ++i) {
        fv.step(dt, 0.0, advection_speed);
    }

    int iterations = 50000;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        fv.step(dt, 0.0, advection_speed);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    double steps_per_sec = iterations / elapsed.count();
    double ns_per_step = (elapsed.count() * 1e9) / iterations;

    std::cout << "FV Full Step Benchmark (" << n_elements << " elements, " << iterations << " iterations)" << std::endl;
    std::cout << "Total Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Throughput: " << steps_per_sec << " steps/s" << std::endl;
    std::cout << "Avg Latency: " << ns_per_step << " ns/step" << std::endl;

    return 0;
}
