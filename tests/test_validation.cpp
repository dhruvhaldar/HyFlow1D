#include <cassert>
#include <stdexcept>
#include <iostream>
#include "fv_solver.hpp"
#include "dg_solver.hpp"

int main() {
    FiniteVolumeSolver fv;
    DiscontinuousGalerkinSolver dg(2);

    // Helper to check if a specific exception type is thrown
    auto check_fail = [&](auto&& f) {
        try {
            f();
            return false;
        }
        catch (const std::invalid_argument&) { return true; }
        catch (const std::length_error&) { return true; } // For excessive memory check
        catch (const std::overflow_error&) { return true; } // For integer overflow
        catch (...) { return false; } // Unexpected exception type
    };

    assert(check_fail([&]{ fv.initialize(0, 1, 0); }));
    assert(check_fail([&]{ fv.initialize(1, 0, 10); }));

    // Security: Test excessive memory limit
    // 50,000,001 elements > 50,000,000 limit
    assert(check_fail([&]{ fv.initialize(0, 1, 50000001); }));

    assert(check_fail([&]{ dg.initialize(0, 1, 0); }));
    assert(check_fail([&]{ DiscontinuousGalerkinSolver bad_dg(-1); }));

    (void)check_fail; // Suppress unused variable warning in Release builds with NDEBUG

    std::cout << "All validation tests passed." << std::endl;

    return 0;
}
