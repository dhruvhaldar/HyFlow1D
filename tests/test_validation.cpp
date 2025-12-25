#include <cassert>
#include <stdexcept>
#include "fv_solver.hpp"
#include "dg_solver.hpp"

int main() {
    FiniteVolumeSolver fv;
    DiscontinuousGalerkinSolver dg(2);

    auto check_fail = [&](auto&& f) {
        try { f(); return false; } catch (const std::invalid_argument&) { return true; }
    };

    assert(check_fail([&]{ fv.initialize(0, 1, 0); }));
    assert(check_fail([&]{ fv.initialize(1, 0, 10); }));
    assert(check_fail([&]{ dg.initialize(0, 1, 0); }));
    assert(check_fail([&]{ DiscontinuousGalerkinSolver bad_dg(-1); }));

    return 0;
}
