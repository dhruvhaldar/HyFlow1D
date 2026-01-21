#include <iostream>
#include <vector>
#include <chrono>
#include "dg_solver.hpp"

// Benchmark for get_solution (evaluate_element)
int main() {
    int n_elements = 5000;
    int p_order = 3;
    DiscontinuousGalerkinSolver dg(p_order);
    dg.initialize(0.0, 1.0, n_elements);
    dg.set_initial_condition([](double x) { return x * x; });

    int iterations = 1000;
    auto start_time = std::chrono::high_resolution_clock::now();

    size_t total_points = 0;
    for (int i = 0; i < iterations; ++i) {
        auto sol = dg.get_solution();
        total_points += sol.size();
        // Prevent dead code elimination
        if (!sol.empty()) {
            volatile double dummy = sol[0].second;
            (void)dummy;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "Calls: " << iterations / elapsed.count() << " /s" << std::endl;
    std::cout << "Points: " << total_points / elapsed.count() / 1e6 << " M/s" << std::endl;

    return 0;
}
