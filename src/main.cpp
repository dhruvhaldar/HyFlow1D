#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <cstring>
#include <filesystem>

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

void show_usage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [OPTIONS]\n"
              << "Options:\n"
              << "  -h, --help      Show this help message\n"
              << "  -v, --verbose   Show detailed step output (disables progress bar)\n"
              << "  -o, --output    Specify output directory (default: output)\n"
              << "\n"
              << "Description:\n"
              << "  Runs a Hybrid FV-DG Simulation for 1D Advection.\n";
}

void draw_progress_bar(int step, int total_steps) {
    float progress = (float)step / total_steps;
    if (progress > 1.0f) progress = 1.0f;
    int barWidth = 50;

    std::cout << "\r[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %";
    std::cout.flush();
}

int main(int argc, char* argv[]) {
    try {
        bool verbose = false;
        std::string output_dir = "output";

        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                show_usage(argv[0]);
                return 0;
            } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
                verbose = true;
            } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
                if (i + 1 < argc) {
                    output_dir = argv[++i];
                } else {
                    std::cerr << "Error: Output directory not specified after " << argv[i] << std::endl;
                    return 1;
                }
            }
        }

        std::cout << "Starting Hybrid FV-DG Simulation..." << std::endl;

        // Create output directory
        namespace fs = std::filesystem;

        // Security Check: Verify output_dir is not a file
        if (fs::exists(output_dir) && !fs::is_directory(output_dir)) {
            throw std::runtime_error("Output path '" + output_dir + "' exists and is not a directory.");
        }

        if (!fs::exists(output_dir)) {
            fs::create_directory(output_dir);
        }

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

        int total_steps = int(t_final/dt);
        std::cout << "dt = " << dt << ", Total Steps = " << total_steps << std::endl;

        while (t < t_final) {
            if (step % output_interval == 0) {
                std::string fname = output_dir + "/solution_" + std::to_string(step) + ".csv";
                write_solution(fname, hybrid.get_solution());
                if (verbose) {
                    std::cout << "Step " << step << ", t = " << t << std::endl;
                }
            }

            if (!verbose) {
                draw_progress_bar(step, total_steps);
            }

            hybrid.step(dt, advection_speed);
            t += dt;
            step++;
        }

        if (!verbose) {
            draw_progress_bar(total_steps, total_steps);
            std::cout << std::endl;
        }

        // Final output
        write_solution(output_dir + "/solution_final.csv", hybrid.get_solution());
        std::cout << "Simulation Complete. Results saved in " << output_dir << "/" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\nFATAL ERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "\nFATAL ERROR: Unknown exception occurred." << std::endl;
        return 1;
    }

    return 0;
}
