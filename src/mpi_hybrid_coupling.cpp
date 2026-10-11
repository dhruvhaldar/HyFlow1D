#ifdef HYFLOW_ENABLE_MPI

#include "mpi_hybrid_coupling.hpp"
#include <stdexcept>
#include <iostream>

MpiHybridDomain::MpiHybridDomain(std::unique_ptr<Solver1D> domain, MPI_Comm communicator)
    : local_domain(std::move(domain)), comm(communicator) {
    if (!local_domain) {
        throw std::invalid_argument("Local domain solver cannot be null.");
    }
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);
}

void MpiHybridDomain::exchange_boundaries() {
    MPI_Request reqs[4];
    int num_reqs = 0;

    double send_right = 0.0, recv_right = 0.0;
    double send_left = 0.0, recv_left = 0.0;

    // Send to and receive from right neighbor (rank + 1)
    if (mpi_rank < mpi_size - 1) {
        send_right = local_domain->get_right_boundary_value();
        MPI_Isend(&send_right, 1, MPI_DOUBLE, mpi_rank + 1, 100, comm, &reqs[num_reqs++]);
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, mpi_rank + 1, 101, comm, &reqs[num_reqs++]);
    }

    // Send to and receive from left neighbor (rank - 1)
    if (mpi_rank > 0) {
        send_left = local_domain->get_left_boundary_value();
        MPI_Isend(&send_left, 1, MPI_DOUBLE, mpi_rank - 1, 101, comm, &reqs[num_reqs++]);
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, mpi_rank - 1, 100, comm, &reqs[num_reqs++]);
    }

    if (num_reqs > 0) {
        MPI_Waitall(num_reqs, reqs, MPI_STATUSES_IGNORE);
    }

    if (mpi_rank < mpi_size - 1) {
        local_domain->set_right_neighbor_value(recv_right);
    }
    if (mpi_rank > 0) {
        local_domain->set_left_neighbor_value(recv_left);
    }
}

void MpiHybridDomain::step_euler(double dt, double advection_speed) {
    exchange_boundaries();
    local_domain->compute_rhs(0.0, advection_speed);
    local_domain->update_state(dt);
}

void MpiHybridDomain::step_ssprk3(double dt, double advection_speed) {
    // Stage 1: u(1) = u(0) + dt * L(u(0))
    exchange_boundaries();
    local_domain->compute_rhs(0.0, advection_speed);
    local_domain->save_stage();
    local_domain->blend_stage(0.0, 1.0, dt);

    // Stage 2: u(2) = 3/4 u(0) + 1/4 (u(1) + dt * L(u(1)))
    exchange_boundaries();
    local_domain->compute_rhs(0.0, advection_speed);
    local_domain->blend_stage(0.75, 0.25, dt);

    // Stage 3: u(n+1) = 1/3 u(0) + 2/3 (u(2) + dt * L(u(2)))
    exchange_boundaries();
    local_domain->compute_rhs(0.0, advection_speed);
    local_domain->blend_stage(1.0 / 3.0, 2.0 / 3.0, dt);
}

void MpiHybridDomain::step(double dt, double advection_speed) {
    if (integrator == HybridDomain::TimeIntegrator::SSPRK3) {
        step_ssprk3(dt, advection_speed);
    } else {
        step_euler(dt, advection_speed);
    }
}

double MpiHybridDomain::compute_total_mass() const {
    double local_mass = local_domain->compute_total_mass();
    double global_mass = 0.0;
    MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global_mass;
}

double MpiHybridDomain::compute_total_energy() const {
    double local_energy = local_domain->compute_total_energy();
    double global_energy = 0.0;
    MPI_Allreduce(&local_energy, &global_energy, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global_energy;
}

std::vector<std::pair<double, double>> MpiHybridDomain::gather_solution() const {
    auto local_sol = local_domain->get_solution();

    std::vector<double> local_buf;
    local_buf.reserve(local_sol.size() * 2);
    for (const auto& pt : local_sol) {
        local_buf.push_back(pt.first);
        local_buf.push_back(pt.second);
    }

    if (mpi_rank != 0) {
        int count = static_cast<int>(local_buf.size());
        MPI_Send(&count, 1, MPI_INT, 0, 200, comm);
        MPI_Send(local_buf.data(), count, MPI_DOUBLE, 0, 201, comm);
        return {};
    }

    // Rank 0 collects from all other ranks
    std::vector<std::pair<double, double>> global_sol = local_sol;

    for (int r = 1; r < mpi_size; ++r) {
        int count = 0;
        MPI_Recv(&count, 1, MPI_INT, r, 200, comm, MPI_STATUS_IGNORE);
        std::vector<double> r_buf(count);
        MPI_Recv(r_buf.data(), count, MPI_DOUBLE, r, 201, comm, MPI_STATUS_IGNORE);

        for (int i = 0; i < count; i += 2) {
            global_sol.push_back({r_buf[i], r_buf[i + 1]});
        }
    }

    return global_sol;
}

int MpiHybridDomain::get_troubled_cell_count(double threshold) const {
    int local_count = local_domain ? local_domain->get_troubled_cell_count(threshold) : 0;
    int global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, comm);
    return global_count;
}

std::vector<double> MpiHybridDomain::gather_modal_sensor() const {
    std::vector<double> local_sensor = local_domain ? local_domain->compute_modal_sensor() : std::vector<double>{};
    if (local_sensor.empty() && local_domain) {
        local_sensor.assign(local_domain->get_num_elements(), -16.0);
    }

    if (mpi_rank != 0) {
        int count = static_cast<int>(local_sensor.size());
        MPI_Send(&count, 1, MPI_INT, 0, 300, comm);
        if (count > 0) {
            MPI_Send(local_sensor.data(), count, MPI_DOUBLE, 0, 301, comm);
        }
        return {};
    }

    std::vector<double> global_sensor = local_sensor;
    for (int r = 1; r < mpi_size; ++r) {
        int count = 0;
        MPI_Recv(&count, 1, MPI_INT, r, 300, comm, MPI_STATUS_IGNORE);
        if (count > 0) {
            std::vector<double> r_buf(count);
            MPI_Recv(r_buf.data(), count, MPI_DOUBLE, r, 301, comm, MPI_STATUS_IGNORE);
            global_sensor.insert(global_sensor.end(), r_buf.begin(), r_buf.end());
        }
    }
    return global_sensor;
}

#endif // HYFLOW_ENABLE_MPI
