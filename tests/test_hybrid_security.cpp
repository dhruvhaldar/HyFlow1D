#include <iostream>
#include <memory>
#include <cassert>
#include <stdexcept>
#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"

// Security Test for Hybrid Domain Connectivity
// Ensures that gaps or overlaps in the simulation domain are rejected.

int main() {
    std::cout << "Running Hybrid Domain Security Tests..." << std::endl;

    // Test Case 1: Valid Configuration
    {
        std::cout << "  Testing Valid Connectivity... ";
        auto fv = std::make_unique<FiniteVolumeSolver>();
        fv->initialize(0.0, 0.5, 10);

        auto dg = std::make_unique<DiscontinuousGalerkinSolver>(1);
        dg->initialize(0.5, 1.0, 5);

        try {
            HybridDomain hybrid(std::move(fv), std::move(dg));
            std::cout << "PASSED" << std::endl;
        } catch (const std::exception& e) {
            std::cout << "FAILED (Unexpected exception: " << e.what() << ")" << std::endl;
            return 1;
        }
    }

    // Test Case 2: Gap Detection
    {
        std::cout << "  Testing Gap Detection... ";
        auto fv = std::make_unique<FiniteVolumeSolver>();
        fv->initialize(0.0, 0.4, 10); // Ends at 0.4

        auto dg = std::make_unique<DiscontinuousGalerkinSolver>(1);
        dg->initialize(0.6, 1.0, 5); // Starts at 0.6

        bool caught = false;
        try {
            HybridDomain hybrid(std::move(fv), std::move(dg));
        } catch (const std::invalid_argument& e) {
            caught = true;
            // Verify message content
            std::string msg = e.what();
            if (msg.find("Domain mismatch") != std::string::npos) {
                // Good
            } else {
                std::cout << "FAILED (Wrong message: " << msg << ")" << std::endl;
                return 1;
            }
        } catch (...) {
            std::cout << "FAILED (Wrong exception type)" << std::endl;
            return 1;
        }

        if (caught) std::cout << "PASSED" << std::endl;
        else {
            std::cout << "FAILED (No exception thrown for gap)" << std::endl;
            return 1;
        }
    }

    // Test Case 3: Overlap Detection
    {
        std::cout << "  Testing Overlap Detection... ";
        auto fv = std::make_unique<FiniteVolumeSolver>();
        fv->initialize(0.0, 0.6, 10); // Ends at 0.6

        auto dg = std::make_unique<DiscontinuousGalerkinSolver>(1);
        dg->initialize(0.5, 1.0, 5); // Starts at 0.5

        bool caught = false;
        try {
            HybridDomain hybrid(std::move(fv), std::move(dg));
        } catch (const std::invalid_argument& e) {
            caught = true;
        }

        if (caught) std::cout << "PASSED" << std::endl;
        else {
            std::cout << "FAILED (No exception thrown for overlap)" << std::endl;
            return 1;
        }
    }

    std::cout << "All Hybrid Security Tests Passed!" << std::endl;
    return 0;
}
