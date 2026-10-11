#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "fv_solver.hpp"
#include "dg_solver.hpp"

#ifdef HYFLOW_ENABLE_OPENMP
#include <omp.h>
#endif

double test_func(double x) {
    return std::sin(2.0 * M_PI * x) + 0.5 * std::cos(4.0 * M_PI * x);
}

int main() {
    std::cout << "Running test_openmp..." << std::endl;

#ifdef HYFLOW_ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    std::cout << "Max OpenMP threads detected: " << max_threads << std::endl;

    // Run FV with 1 thread
    omp_set_num_threads(1);
    FiniteVolumeSolver fv1;
    fv1.initialize(0.0, 1.0, 1000);
    fv1.set_initial_condition(test_func);
    fv1.compute_rhs(0.0, 1.0);
    double mass1 = fv1.compute_total_mass();
    double energy1 = fv1.compute_total_energy();

    // Run FV with multiple threads
    int test_threads = std::min(4, max_threads);
    omp_set_num_threads(test_threads);
    FiniteVolumeSolver fvN;
    fvN.initialize(0.0, 1.0, 1000);
    fvN.set_initial_condition(test_func);
    fvN.compute_rhs(0.0, 1.0);
    double massN = fvN.compute_total_mass();
    double energyN = fvN.compute_total_energy();

    assert(std::abs(mass1 - massN) < 1e-12);
    assert(std::abs(energy1 - energyN) < 1e-12);

    // Run DG with 1 thread vs N threads
    omp_set_num_threads(1);
    DiscontinuousGalerkinSolver dg1(3);
    dg1.initialize(0.0, 1.0, 200);
    dg1.set_initial_condition(test_func);
    dg1.compute_rhs(0.0, 1.0);
    double dg_mass1 = dg1.compute_total_mass();
    double dg_energy1 = dg1.compute_total_energy();

    omp_set_num_threads(test_threads);
    DiscontinuousGalerkinSolver dgN(3);
    dgN.initialize(0.0, 1.0, 200);
    dgN.set_initial_condition(test_func);
    dgN.compute_rhs(0.0, 1.0);
    double dg_massN = dgN.compute_total_mass();
    double dg_energyN = dgN.compute_total_energy();

    assert(std::abs(dg_mass1 - dg_massN) < 1e-12);
    assert(std::abs(dg_energy1 - dg_energyN) < 1e-12);

    std::cout << "OpenMP thread equivalence confirmed across 1 and " << test_threads << " threads." << std::endl;
#else
    std::cout << "OpenMP not enabled at compile time. Skipping thread check." << std::endl;
#endif

    std::cout << "test_openmp passed successfully." << std::endl;
    return 0;
}
