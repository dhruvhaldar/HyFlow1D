#pragma once

#ifdef HYFLOW_ENABLE_MPI

#include <mpi.h>
#include <memory>
#include <vector>
#include <utility>
#include "solver.hpp"
#include "hybrid_coupling.hpp"

// Distributed Hybrid Domain using MPI for inter-process domain decomposition
// Rank 0 typically owns Subdomain 0 (FV), Rank 1 typically owns Subdomain 1 (DG).
class MpiHybridDomain {
private:
    std::unique_ptr<Solver1D> local_domain;
    int mpi_rank;
    int mpi_size;
    MPI_Comm comm;

    HybridDomain::TimeIntegrator integrator = HybridDomain::TimeIntegrator::SSPRK3;

public:
    MpiHybridDomain(std::unique_ptr<Solver1D> domain, MPI_Comm communicator = MPI_COMM_WORLD);
    ~MpiHybridDomain() = default;

    int get_rank() const { return mpi_rank; }
    int get_size() const { return mpi_size; }
    Solver1D* get_local_domain() { return local_domain.get(); }
    const Solver1D* get_local_domain() const { return local_domain.get(); }

    void set_time_integrator(HybridDomain::TimeIntegrator ti) { integrator = ti; }
    HybridDomain::TimeIntegrator get_time_integrator() const { return integrator; }

    // Non-blocking MPI halo exchange between neighboring subdomains
    void exchange_boundaries();

    // Perform one time step
    void step(double dt, double advection_speed);
    void step_euler(double dt, double advection_speed);
    void step_ssprk3(double dt, double advection_speed);

    // Distributed reductions across all MPI ranks
    double compute_total_mass() const;
    double compute_total_energy() const;

    // Gather solution from all ranks onto Rank 0 (for unified I/O)
    std::vector<std::pair<double, double>> gather_solution() const;
};

#endif // HYFLOW_ENABLE_MPI
