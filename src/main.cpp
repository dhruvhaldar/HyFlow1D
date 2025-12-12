#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <string>

#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"

// Initial Condition: Gaussian Pulse
double initial_condition(double x) {
    double center = 0.25;
    double width = 0.05;
    return std::exp(-(x - center) * (x - center) / (2 * width * width));
}

void write_solution(const std::string& filename, const std::vector<std::pair<double, double>>& solution) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outfile << "x,u\n";
    for (const auto& p : solution) {
        outfile << p.first << "," << p.second << "\n";
    }
    outfile.close();
}

int main() {
    std::cout << "Starting Hybrid FV-DG Simulation..." << std::endl;

    // Simulation Parameters
    double x_start = 0.0;
    double x_interface = 0.5; // Interface at 0.5
    double x_end = 1.0;
    
    // FV Region: [0, 0.5], 50 elements -> dx = 0.01
    int n_fv = 50; 
    // DG Region: [0.5, 1.0], 10 elements -> dx = 0.05. P=3.
    // Degrees of freedom roughly: 10 * 4 = 40. Comparable resolution? 
    // FV is 50 dofs. DG is 40 dofs.
    int n_dg = 10;
    int p_order = 3; 

    double advection_speed = 1.0;
    double cfl = 0.1; // Low CFL for stability with Forward Euler
    double dx_min = std::min( (x_interface - x_start)/n_fv, (x_end - x_interface)/n_dg / (2*p_order+1) );
    double dt = cfl * dx_min / advection_speed;
    double t_final = 0.6; // Pulse should move from 0.25 to 0.85 (well into DG region)

    // Setup Solvers
    auto fv = std::make_unique<FiniteVolumeSolver>();
    fv->initialize(x_start, x_interface, n_fv);
    fv->set_initial_condition(initial_condition);

    auto dg = std::make_unique<DiscontinuousGalerkinSolver>(p_order);
    dg->initialize(x_interface, x_end, n_dg);
    dg->set_initial_condition(initial_condition); // Should be approx 0

    // Setup Hybrid Domain
    HybridDomain hybrid(std::move(fv), std::move(dg));

    // Time Loop
    double t = 0.0;
    int step = 0;
    int output_interval = 100;
    
    std::cout << "dt = " << dt << ", Total Steps = " << int(t_final/dt) << std::endl;

    while (t < t_final) {
        if (step % output_interval == 0) {
            std::string fname = "solution_" + std::to_string(step) + ".csv";
            write_solution(fname, hybrid.get_solution());
            std::cout << "Step " << step << ", t = " << t << std::endl;
        }

        hybrid.step(dt, advection_speed);
        t += dt;
        step++;
    }

    // Final output
    write_solution("solution_final.csv", hybrid.get_solution());
    std::cout << "Simulation Complete." << std::endl;

    return 0;
}
