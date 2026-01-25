#include <iostream>
#include <cassert>
#include <stdexcept>
#include <limits>
#include "fv_solver.hpp"
#include "dg_solver.hpp"

// Security Test Suite
// Verifies that solvers enforce limits on resource usage to prevent DoS.

int main() {
    std::cout << "Running Security Tests..." << std::endl;

    constexpr int OVER_LIMIT = Solver1D::MAX_ELEMENTS + 1;

    // Test FV Solver Limits
    {
        std::cout << "  Testing FV Solver DoS protection... ";
        FiniteVolumeSolver fv;
        bool caught = false;
        try {
            fv.initialize(0.0, 1.0, OVER_LIMIT);
        } catch (const std::length_error& e) {
            caught = true;
        } catch (...) {
            std::cout << "FAILED (Wrong exception type)" << std::endl;
            return 1;
        }

        if (caught) std::cout << "PASSED" << std::endl;
        else {
            std::cout << "FAILED (No exception thrown)" << std::endl;
            return 1;
        }
    }

    // Test DG Solver Constructor Validation (Integer Overflow/UB prevention)
    {
        std::cout << "  Testing DG Solver Constructor Validation... ";
        bool caught = false;
        try {
            // Attempt to create solver with INT_MAX to trigger potential overflow
            DiscontinuousGalerkinSolver dg(std::numeric_limits<int>::max());
        } catch (const std::invalid_argument& e) {
            caught = true;
        } catch (...) {
            std::cout << "FAILED (Wrong exception type)" << std::endl;
            return 1;
        }

        if (caught) std::cout << "PASSED" << std::endl;
        else {
            std::cout << "FAILED (No exception thrown)" << std::endl;
            return 1;
        }
    }

    // Test DG Solver Limits
    {
        std::cout << "  Testing DG Solver DoS protection... ";
        DiscontinuousGalerkinSolver dg(1); // P=1
        bool caught = false;
        try {
            dg.initialize(0.0, 1.0, OVER_LIMIT);
        } catch (const std::length_error& e) {
            caught = true;
        } catch (...) {
            std::cout << "FAILED (Wrong exception type)" << std::endl;
            return 1;
        }

        if (caught) std::cout << "PASSED" << std::endl;
        else {
            std::cout << "FAILED (No exception thrown)" << std::endl;
            return 1;
        }
    }

    std::cout << "All Security Tests Passed!" << std::endl;
    return 0;
}
