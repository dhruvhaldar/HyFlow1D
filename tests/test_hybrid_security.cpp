#include <iostream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "hybrid_coupling.hpp"
#include "fv_solver.hpp"
#include "dg_solver.hpp"

// Security Test for HybridDomain
// Ensures that HybridDomain rejects null pointers to prevent Segfaults (DoS).

int main() {
    std::cout << "Running Hybrid Domain Security Tests..." << std::endl;

    // Helper to check if HybridDomain throws invalid_argument on null pointers
    auto check_null_rejection = [](std::unique_ptr<Solver1D> left, std::unique_ptr<Solver1D> right, const std::string& case_name) {
        std::cout << "  Testing " << case_name << "... ";
        try {
            HybridDomain hybrid(std::move(left), std::move(right));
            // If we get here, no exception was thrown
            std::cout << "FAILED (No exception thrown)" << std::endl;
            return false;
        } catch (const std::invalid_argument& e) {
            std::cout << "PASSED (Caught expected exception: " << e.what() << ")" << std::endl;
            return true;
        } catch (...) {
            std::cout << "FAILED (Caught wrong exception type)" << std::endl;
            return false;
        }
    };

    bool all_passed = true;

    // Case 1: Left solver is null
    {
        auto right = std::make_unique<DiscontinuousGalerkinSolver>(2);
        if (!check_null_rejection(nullptr, std::move(right), "Left=nullptr")) {
            all_passed = false;
        }
    }

    // Case 2: Right solver is null
    {
        auto left = std::make_unique<FiniteVolumeSolver>();
        if (!check_null_rejection(std::move(left), nullptr, "Right=nullptr")) {
            all_passed = false;
        }
    }

    // Case 3: Both solvers are null
    {
        if (!check_null_rejection(nullptr, nullptr, "Both=nullptr")) {
            all_passed = false;
        }
    }

    if (all_passed) {
        std::cout << "All Hybrid Security Tests Passed!" << std::endl;
        return 0;
    } else {
        std::cerr << "Some tests failed." << std::endl;
        return 1;
    }
}
