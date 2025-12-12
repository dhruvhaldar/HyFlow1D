#include <iostream>
#include <cassert>
#include <cmath>
#include "fv_solver.hpp"
#include "dg_solver.hpp"

// Simple linear advection test
double initial_condition(double x) {
    return std::exp(-100.0 * (x - 0.5) * (x - 0.5));
}

int main() {
    std::cout << "Running Solver Tests..." << std::endl;

    // Test FV Solver
    {
        FiniteVolumeSolver fv;
        fv.initialize(0.0, 1.0, 100);
        fv.set_initial_condition(initial_condition);
        
        // Check conservation (sum of u should be approx constant if periodic)
        // But here we have 0 boundary conditions by default ghost 0
        // Just check it initializes
        auto sol = fv.get_solution();
        assert(sol.size() == 100);
        std::cout << "FV Initialized correctly." << std::endl;
    }

    // Test DG Solver
    {
        DiscontinuousGalerkinSolver dg(2); // P=2
        dg.initialize(0.0, 1.0, 20);
        dg.set_initial_condition(initial_condition);
        
        auto sol = dg.get_solution();
        // 20 elements * 6 points each (0 to 5)
        assert(sol.size() == 20 * 6); 
        std::cout << "DG Initialized correctly." << std::endl;
        
        // Check value at peak approx 1.0
        double max_val = 0.0;
        for(auto p : sol) if(p.second > max_val) max_val = p.second;
        assert(std::abs(max_val - 1.0) < 0.1); // Projection error
        std::cout << "DG Peak value check passed: " << max_val << std::endl;
    }

    std::cout << "All Solver Tests Passed!" << std::endl;
    return 0;
}
