#include <iostream>
#include <cmath>
#include <memory>
#include <cassert>

#ifdef HYFLOW_ENABLE_MPI
#include <mpi.h>
#include "mpi_hybrid_coupling.hpp"
#include "fv_solver.hpp"
#include "dg_solver.hpp"

double initial_pulse(double x) {
    double center = 0.25;
    double width = 0.05;
    return std::exp(-(x - center) * (x - center) / (2.0 * width * width));
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Running test_mpi with " << size << " ranks..." << std::endl;
    }

    if (size != 2) {
        if (rank == 0) {
            std::cout << "test_mpi requires exactly 2 ranks. Skipping." << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    std::unique_ptr<Solver1D> local_domain;
    if (rank == 0) {
        // Rank 0 owns FV [0.0, 0.5]
        auto fv = std::make_unique<FiniteVolumeSolver>();
        fv->initialize(0.0, 0.5, 50);
        fv->set_initial_condition(initial_pulse);
        local_domain = std::move(fv);
    } else {
        // Rank 1 owns DG [0.5, 1.0]
        auto dg = std::make_unique<DiscontinuousGalerkinSolver>(3);
        dg->initialize(0.5, 1.0, 10);
        dg->set_initial_condition(initial_pulse);
        local_domain = std::move(dg);
    }

    MpiHybridDomain mpi_hybrid(std::move(local_domain), MPI_COMM_WORLD);
    mpi_hybrid.set_time_integrator(HybridDomain::TimeIntegrator::SSPRK3);

    // Initial global mass and energy
    double m0 = mpi_hybrid.compute_total_mass();
    double e0 = mpi_hybrid.compute_total_energy();

    if (rank == 0) {
        std::cout << "MPI Initial Mass: " << m0 << ", Initial Energy: " << e0 << std::endl;
        assert(m0 > 0.0);
        assert(e0 > 0.0);
    }

    // Step a few iterations
    double dt = 0.0005;
    double a = 1.0;
    for (int step = 0; step < 50; ++step) {
        mpi_hybrid.step(dt, a);
    }

    double m_final = mpi_hybrid.compute_total_mass();
    double e_final = mpi_hybrid.compute_total_energy();

    if (rank == 0) {
        std::cout << "MPI Final Mass: " << m_final << ", Final Energy: " << e_final << std::endl;
        assert(std::abs(m_final - m0) / m0 < 1e-4);
        assert(e_final <= e0 + 1e-12);
    }

    auto full_sol = mpi_hybrid.gather_solution();
    if (rank == 0) {
        assert(!full_sol.empty());
        std::cout << "Gathered " << full_sol.size() << " coordinate points from all ranks." << std::endl;
        std::cout << "test_mpi passed successfully on Rank 0." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
#else
int main() {
    std::cout << "HYFLOW_ENABLE_MPI not defined. Skipping test_mpi." << std::endl;
    return 0;
}
#endif
