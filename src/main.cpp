#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <cstring>
#include <filesystem>
#include <chrono>
#include <iomanip>

#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"
#include "term_colors.hpp"

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
    std::cout << "Usage: " << Color::Bold << prog_name << " [OPTIONS]" << Color::Reset << "\n"
              << "Options:\n"
              << "  " << Color::Yellow << "-h, --help" << Color::Reset << "      Show this help message\n"
              << "  " << Color::Yellow << "-v, --verbose" << Color::Reset << "   Show detailed step output (disables progress bar)\n"
              << "  " << Color::Yellow << "-o, --output" << Color::Reset << "    Specify output directory (default: output)\n"
              << "\n"
              << "Description:\n"
              << "  Runs a Hybrid FV-DG Simulation for 1D Advection.\n";
}

void draw_progress_bar(int step, int total_steps, double elapsed_seconds, double t_current, double t_final) {
    float progress = (float)step / total_steps;
    if (progress > 1.0f) progress = 1.0f;
    int barWidth = 50;

    // Spinner animation
    const char spinner[] = {'|', '/', '-', '\\'};
    char spin_char = (progress >= 1.0f) ? ' ' : spinner[step % 4];

    std::cout << "\r" << Color::Blue << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "█";
        else if (i == pos) std::cout << "█";
        else std::cout << " ";
    }
    std::cout << "] " << Color::Reset << Color::Bold << int(progress * 100.0) << " % " << Color::Reset;

    // Simulation time info
    std::cout << Color::Cyan << " (t=" << std::fixed << std::setprecision(2) << t_current << "/" << t_final << ")" << Color::Reset << " ";

    // Spinner & ETA
    if (progress < 1.0f) {
        std::cout << spin_char << " ";
        if (progress > 0.0f) {
            double total_estimated = elapsed_seconds / progress;
            double remaining = total_estimated - elapsed_seconds;
            int r_min = (int)remaining / 60;
            int r_sec = (int)remaining % 60;
            std::cout << "ETA: " << std::setfill('0') << std::setw(2) << r_min << ":"
                      << std::setw(2) << r_sec;
        }
    }

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
            } else {
                // Check if it looks like a flag
                if (argv[i][0] == '-') {
                    std::cerr << Color::BoldRed << "Error: Unknown option '" << argv[i] << "'" << Color::Reset << std::endl;
                    show_usage(argv[0]);
                    return 1;
                } else {
                     // Could be a positional argument, but we don't support any yet.
                    std::cerr << Color::BoldRed << "Error: Unexpected argument '" << argv[i] << "'" << Color::Reset << std::endl;
                    show_usage(argv[0]);
                    return 1;
                }
            }
        }

        std::cout << Color::BoldCyan << "Starting Hybrid FV-DG Simulation..." << Color::Reset << std::endl;

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
        std::cout << Color::Yellow << "dt = " << dt << ", Total Steps = " << total_steps << Color::Reset << std::endl;

        auto start_time = std::chrono::steady_clock::now();

        while (t < t_final) {
            if (step % output_interval == 0) {
                std::string fname = output_dir + "/solution_" + std::to_string(step) + ".csv";
                write_solution(fname, hybrid.get_solution());
                if (verbose) {
                    std::cout << "Step " << step << ", t = " << t << std::endl;
                }
            }

            if (!verbose) {
                auto now = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed = now - start_time;
                draw_progress_bar(step, total_steps, elapsed.count(), t, t_final);
            }

            hybrid.step(dt, advection_speed);
            t += dt;
            step++;
        }

        if (!verbose) {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - start_time;
            draw_progress_bar(total_steps, total_steps, elapsed.count(), t_final, t_final);
            std::cout << std::endl;
        }

        // Final output
        write_solution(output_dir + "/solution_final.csv", hybrid.get_solution());
        std::cout << Color::BoldGreen << "Simulation Complete. Results saved in " << output_dir << "/" << Color::Reset << std::endl;

    } catch (const std::exception& e) {
        std::cerr << Color::BoldRed << "\nFATAL ERROR: " << e.what() << Color::Reset << std::endl;
        return 1;
    } catch (...) {
        std::cerr << Color::BoldRed << "\nFATAL ERROR: Unknown exception occurred." << Color::Reset << std::endl;
        return 1;
    }

    return 0;
}
