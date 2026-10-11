#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <cstring>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <locale>
#include <csignal>
#include <algorithm>
#include <numeric>
#include <cerrno>
#include <cctype>

#ifdef HYFLOW_ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef HYFLOW_ENABLE_MPI
#include <mpi.h>
#include "mpi_hybrid_coupling.hpp"
#endif

// Custom numpunct to add thousands separator (comma)
struct comma_numpunct : std::numpunct<char> {
    char do_thousands_sep() const override { return ','; }
    std::string do_grouping() const override { return "\3"; }
};

// Check for unistd.h availability for isatty
#if defined(__unix__) || defined(__APPLE__) || defined(__linux__)
    #include <unistd.h>
    #include <fcntl.h>
    #define HAS_UNISTD 1
#elif defined(_WIN32) || defined(_WIN64)
    #include <io.h>
    #define HAS_IO_H 1
#endif

#include <sys/stat.h> // For umask/ _umask
#include "fv_solver.hpp"
#include "dg_solver.hpp"
#include "hybrid_coupling.hpp"
#include "term_colors.hpp"

// Signal Handler for Graceful Shutdown
volatile std::sig_atomic_t g_signal_status = 0;

void signal_handler(int signum) {
    if (signum == SIGINT) {
        g_signal_status = signum;
    }
}

// Initial Condition: Gaussian Pulse
double initial_condition(double x) {
    double center = 0.25;
    double width = 0.05;
    return std::exp(-(x - center) * (x - center) / (2.0 * width * width));
}

// Security: Validate path to prevent directory traversal
bool is_safe_path(const std::string& path_str) {
    if (path_str.empty()) return false;
    std::filesystem::path path(path_str);
    for (const auto& part : path) {
        if (part == "..") return false;
    }
    return true;
}

void write_solution(const std::string& filename, const std::vector<std::pair<double, double>>& solution, double time = -1.0, const std::string& extra_header = "") {
#if defined(HAS_UNISTD)
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_NOFOLLOW | O_CLOEXEC, 0600);

    if (fd == -1) {
        if (errno == ELOOP) {
            std::cerr << Color::BoldRed << "Error: Output file '" << filename << "' is a symbolic link." << Color::Reset << std::endl;
            std::cerr << "Refusing to overwrite symbolic links to prevent security risks." << std::endl;
        } else {
            std::cerr << Color::BoldRed << "Error opening file '" << filename << "': " << std::strerror(errno) << Color::Reset << std::endl;
        }
        return;
    }

    FILE* outfile = fdopen(fd, "w");
    if (!outfile) {
        std::cerr << "Error: fdopen failed." << std::endl;
        close(fd);
        return;
    }

    struct FileDeleter {
        void operator()(FILE* f) const { std::fclose(f); }
    };
    std::unique_ptr<FILE, FileDeleter> file_guard(outfile);

    if (time >= 0.0) {
        std::fprintf(outfile, "# t=%g\n", time);
    }
    if (!extra_header.empty()) {
        std::fprintf(outfile, "# %s\n", extra_header.c_str());
    }

    std::fprintf(outfile, "x,u\n");
    for (const auto& pt : solution) {
        std::fprintf(outfile, "%.17g,%.17g\n", pt.first, pt.second);
    }

#else
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    if (time >= 0.0) {
        outfile << "# t=" << time << "\n";
    }
    if (!extra_header.empty()) {
        outfile << "# " << extra_header << "\n";
    }

    outfile << "x,u\n";
    for (const auto& pt : solution) {
        outfile << std::setprecision(16) << pt.first << "," << pt.second << "\n";
    }
    outfile.close();

    try {
        std::filesystem::permissions(filename,
            std::filesystem::perms::owner_read | std::filesystem::perms::owner_write,
            std::filesystem::perm_options::replace);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << Color::Yellow << "Warning: Failed to set secure permissions on '" << filename << "': " << e.what() << Color::Reset << std::endl;
    }
#endif
}

void show_usage(const char* prog_name) {
    std::cout << "\n" << Color::BoldCyan << "  HyFlow1D  " << Color::Reset << " - High-Performance Hybrid FV-DG Aerodynamic Solver\n\n"
              << "  " << Color::Bold << "Usage:" << Color::Reset << " " << prog_name << " [OPTIONS]\n\n"
              << "  " << Color::Bold << "Options:" << Color::Reset << "\n"
              << "    " << Color::Yellow << "-h, --help" << Color::Reset << "            Show this help message\n"
              << "    " << Color::Yellow << "-v, --verbose" << Color::Reset << "         Show detailed step output (disables progress bar)\n"
              << "    " << Color::Yellow << "-o, --output <dir>" << Color::Reset << "   Specify output directory (default: output)\n"
              << "    " << Color::Yellow << "-c, --clean" << Color::Reset << "           Clean output directory before running\n"
              << "    " << Color::Yellow << "-n, --dry-run" << Color::Reset << "         Print configuration and exit without running\n"
              << "    " << Color::Yellow << "-t, --threads <N>" << Color::Reset << "     Set number of OpenMP threads\n"
              << "    " << Color::Yellow << "-i, --integrator <name>" << Color::Reset << " Select time integrator ('ssprk3' or 'euler', default: ssprk3)\n\n"
              << "  " << Color::Bold << "Simulation Details:" << Color::Reset << "\n"
              << "    Domain: [0, 1] (FV: [0, 0.5], DG: [0.5, 1])\n"
              << "    Config: 50 FV cells, 10 DG elements (P=3)\n"
              << "    Physics: Linear Advection (a=1.0), Gaussian Pulse\n"
              << "    Parallelism: MPI Distributed Domain Decomposition + OpenMP Multithreading\n\n"
              << "  " << Color::Bold << "Examples:" << Color::Reset << "\n"
              << "    " << prog_name << "\n"
              << "    " << prog_name << " -t 4 --integrator ssprk3\n"
              << "    mpirun -n 2 " << prog_name << " -o output_mpi\n";
}

int levenshtein_distance(std::string_view s1, std::string_view s2) {
    const size_t m = s1.length();
    const size_t n = s2.length();
    if (m == 0) return n;
    if (n == 0) return m;

    std::vector<int> costs(n + 1);
    std::iota(costs.begin(), costs.end(), 0);

    for (size_t i = 0; i < m; ++i) {
        int prev = costs[0];
        costs[0] = i + 1;
        for (size_t j = 0; j < n; ++j) {
            int temp = costs[j + 1];
            costs[j + 1] = std::min({ prev + (s1[i] == s2[j] ? 0 : 1), costs[j] + 1, costs[j + 1] + 1 });
            prev = temp;
        }
    }
    return costs[n];
}

void suggest_flag(const std::string& invalid_flag) {
    const std::vector<std::string> valid_flags = {
        "--help", "-h", "--verbose", "-v", "--output", "-o",
        "--clean", "-c", "--dry-run", "-n", "--threads", "-t",
        "--integrator", "-i"
    };
    std::string best_match;
    int min_dist = 100;

    for (const auto& flag : valid_flags) {
        int dist = levenshtein_distance(invalid_flag, flag);
        if (dist < min_dist) {
            min_dist = dist;
            best_match = flag;
        }
    }

    if (min_dist <= 3) {
         std::cerr << "       Did you mean '" << Color::Yellow << best_match << Color::Reset << "'?" << std::endl;
    }
}

std::string shell_escape(const std::string& arg) {
    if (arg.empty()) return "''";
    bool safe = true;
    for (char c : arg) {
        if (!isalnum(c) && c != '-' && c != '_' && c != '.' && c != '/') {
            safe = false;
            break;
        }
    }
    if (safe) return arg;

    std::string result = "'";
    for (char c : arg) {
        if (c == '\'') result += "'\\''";
        else result += c;
    }
    result += "'";
    return result;
}

namespace {
    struct CursorGuard {
        bool active;
        CursorGuard(bool enable) : active(enable) {
            if (active) std::cout << Color::CursorHide << std::flush;
        }
        ~CursorGuard() {
            if (active) std::cout << Color::CursorShow << std::flush;
        }
    };
}

void draw_progress_bar(int step, int total_steps, double elapsed_seconds, double t_current, double t_final) {
    float progress = (float)step / total_steps;
    if (progress > 1.0f) progress = 1.0f;

    const int barWidth = 30;
    const char spinner[] = {'|', '/', '-', '\\'};
    char spin_char = (progress >= 1.0f) ? ' ' : spinner[step % 4];

    std::cout << "\r" << Color::Blue << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i <= pos) std::cout << "█";
        else std::cout << " ";
    }
    std::cout << "] " << Color::Reset << Color::Bold << int(progress * 100.0) << " % " << Color::Reset;
    std::cout << Color::Cyan << " (t=" << std::fixed << std::setprecision(2) << t_current << "/" << t_final << ")" << Color::Reset << " ";

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
#ifdef HYFLOW_ENABLE_MPI
    MPI_Init(&argc, &argv);
    int mpi_rank = 0, mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
    int mpi_rank = 0, mpi_size = 1;
#endif

#if defined(_WIN32) || defined(_WIN64)
    _umask(0077);
#else
    umask(0077);
#endif

    std::signal(SIGINT, signal_handler);

    try {
        namespace fs = std::filesystem;
        bool verbose = false;
        bool clean_output = false;
        bool dry_run = false;
        int num_threads = 0;
        std::string integrator_str = "ssprk3";
        std::string output_dir = "output";

        bool is_tty = (mpi_rank == 0);
#if defined(HAS_UNISTD)
        if (!isatty(fileno(stdout))) is_tty = false;
#elif defined(HAS_IO_H)
        if (!_isatty(_fileno(stdout))) is_tty = false;
#endif

        if (!is_tty || mpi_rank != 0) {
            Color::enabled = false;
        } else {
            std::locale comma_locale(std::locale::classic(), new comma_numpunct());
            std::cout.imbue(comma_locale);
        }

        if (std::getenv("NO_COLOR") != nullptr) {
            Color::enabled = false;
        }

        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                if (mpi_rank == 0) show_usage(argv[0]);
#ifdef HYFLOW_ENABLE_MPI
                MPI_Finalize();
#endif
                return 0;
            } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
                verbose = true;
            } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
                if (i + 1 < argc) {
                    output_dir = argv[++i];
                    if (!is_safe_path(output_dir)) {
                        if (mpi_rank == 0) std::cerr << Color::BoldRed << "Error: Invalid output directory. Path traversal ('..') is not allowed." << Color::Reset << std::endl;
#ifdef HYFLOW_ENABLE_MPI
                        MPI_Finalize();
#endif
                        return 1;
                    }
                } else {
                    if (mpi_rank == 0) std::cerr << "Error: Output directory not specified after " << argv[i] << std::endl;
#ifdef HYFLOW_ENABLE_MPI
                    MPI_Finalize();
#endif
                    return 1;
                }
            } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--clean") == 0) {
                clean_output = true;
            } else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--dry-run") == 0) {
                dry_run = true;
            } else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--threads") == 0) {
                if (i + 1 < argc) {
                    num_threads = std::atoi(argv[++i]);
                }
            } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--integrator") == 0) {
                if (i + 1 < argc) {
                    integrator_str = argv[++i];
                }
            } else {
                if (argv[i][0] == '-') {
                    if (mpi_rank == 0) {
                        std::cerr << Color::BoldRed << "Error: Unknown option '" << argv[i] << "'" << Color::Reset << std::endl;
                        suggest_flag(argv[i]);
                    }
                } else {
                    if (mpi_rank == 0) {
                        std::cerr << Color::BoldRed << "Error: Unexpected argument '" << argv[i] << "'" << Color::Reset << std::endl;
                    }
                }
#ifdef HYFLOW_ENABLE_MPI
                MPI_Finalize();
#endif
                return 1;
            }
        }

#ifdef HYFLOW_ENABLE_OPENMP
        if (num_threads > 0) {
            omp_set_num_threads(num_threads);
        }
        int active_threads = omp_get_max_threads();
#else
        int active_threads = 1;
#endif

        HybridDomain::TimeIntegrator chosen_integrator = HybridDomain::TimeIntegrator::SSPRK3;
        if (integrator_str == "euler" || integrator_str == "forward_euler") {
            chosen_integrator = HybridDomain::TimeIntegrator::ForwardEuler;
        }

        // Only rank 0 prepares the output directory
        if (mpi_rank == 0) {
            if (fs::exists(output_dir) && !fs::is_directory(output_dir)) {
                throw std::runtime_error("Output path '" + output_dir + "' exists and is not a directory.");
            }

            if (!fs::exists(output_dir)) {
                if (fs::create_directory(output_dir)) {
                    fs::permissions(output_dir, fs::perms::owner_all, fs::perm_options::replace | fs::perm_options::nofollow);
                }
            } else {
                if (!fs::is_symlink(output_dir)) {
                    auto perms = fs::status(output_dir).permissions();
                    if ((perms & (fs::perms::group_all | fs::perms::others_all)) != fs::perms::none) {
                        std::cerr << Color::Yellow << "Warning: Output directory '" << output_dir << "' has insecure permissions." << Color::Reset << std::endl;
                        try {
                            fs::permissions(output_dir, fs::perms::owner_all, fs::perm_options::replace | fs::perm_options::nofollow);
                        } catch (...) {}
                    }
                }

                bool has_solution_files = false;
                std::vector<fs::path> files_to_remove;
                for (const auto& entry : fs::directory_iterator(output_dir)) {
                    if (entry.is_regular_file()) {
                        std::string fname = entry.path().filename().string();
                        if (fname.rfind("solution_", 0) == 0 && entry.path().extension() == ".csv") {
                            has_solution_files = true;
                            if (clean_output) files_to_remove.push_back(entry.path());
                        }
                    }
                }

                if (clean_output && has_solution_files) {
                    for (const auto& p : files_to_remove) fs::remove(p);
                    std::cout << Color::Blue << "Cleaned " << files_to_remove.size() << " file(s) from '" << output_dir << "'." << Color::Reset << std::endl;
                } else if (has_solution_files) {
                    std::cout << Color::Yellow << "Warning: Output directory '" << output_dir << "' contains existing solution files. Use --clean to remove them." << Color::Reset << "\n" << std::endl;
                }
            }
        }

        double x_start = 0.0;
        double x_interface = 0.5;
        double x_end = 1.0;
        int n_fv = 50;
        int n_dg = 10;
        int p_order = 3;

        std::ostringstream meta_ss;
        meta_ss << "x_domain=" << x_start << ":" << x_interface << ":" << x_end;
        std::string domain_meta = meta_ss.str();

        double advection_speed = 1.0;
        double cfl = 0.1;
        double dx_min = std::min((x_interface - x_start)/n_fv, (x_end - x_interface)/n_dg / (2*p_order+1));
        double dt = cfl * dx_min / advection_speed;
        double t_final = 0.6;
        int total_steps = int(t_final / dt);

        if (mpi_rank == 0) {
            std::cout << Color::BoldCyan << "Starting Hybrid FV-DG Simulation..." << Color::Reset << "\n\n"
                      << Color::Bold << "Simulation Configuration:" << Color::Reset << "\n"
                      << "  " << Color::Blue << "Domain:       " << Color::Reset << "[" << std::fixed << std::setprecision(2) << x_start << "] |---FV---| [" << x_interface << "] |---DG---| [" << x_end << "]\n"
                      << "  " << Color::Blue << "FV Cells:     " << Color::Reset << Color::Bold << n_fv << Color::Reset << "\n"
                      << "  " << Color::Blue << "DG Elements:  " << Color::Reset << Color::Bold << n_dg << Color::Reset << " (Order P=" << Color::Bold << p_order << Color::Reset << ")\n"
                      << "  " << Color::Blue << "Time:         " << Color::Reset << Color::Bold << "0.0 -> " << t_final << Color::Reset << "\n"
                      << "  " << Color::Blue << "Time Step:    " << Color::Reset << Color::Bold << std::defaultfloat << std::setprecision(6) << dt << Color::Reset << " (" << total_steps << " steps)\n"
                      << "  " << Color::Blue << "Integrator:   " << Color::Reset << Color::Bold << (chosen_integrator == HybridDomain::TimeIntegrator::SSPRK3 ? "SSPRK3 (3rd-Order SSP Runge-Kutta)" : "Forward Euler (1st-Order)") << Color::Reset << "\n"
                      << "  " << Color::Blue << "Parallelism:  " << Color::Reset << "MPI Ranks: " << Color::Bold << mpi_size << Color::Reset << " | OpenMP Threads: " << Color::Bold << active_threads << Color::Reset << "\n"
                      << std::endl;
        }

        if (dry_run) {
            if (mpi_rank == 0) {
                std::cout << Color::Green << "Dry run complete." << Color::Reset << " Configuration valid.\n"
                          << "  (No simulation performed, no files written)" << std::endl;
            }
#ifdef HYFLOW_ENABLE_MPI
            MPI_Finalize();
#endif
            return 0;
        }

        auto start_time = std::chrono::steady_clock::now();
        CursorGuard cursor_guard(is_tty && !verbose);

        double initial_mass = 0.0, initial_energy = 0.0;
        double final_mass = 0.0, final_energy = 0.0;

#ifdef HYFLOW_ENABLE_MPI
        if (mpi_size > 1) {
            // Distributed MPI Hybrid Execution
            std::unique_ptr<Solver1D> local_solver;
            if (mpi_rank == 0) {
                auto fv = std::make_unique<FiniteVolumeSolver>();
                fv->initialize(x_start, x_interface, n_fv);
                fv->set_initial_condition(initial_condition);
                local_solver = std::move(fv);
            } else {
                auto dg = std::make_unique<DiscontinuousGalerkinSolver>(p_order);
                dg->initialize(x_interface, x_end, n_dg);
                dg->set_initial_condition(initial_condition);
                local_solver = std::move(dg);
            }

            MpiHybridDomain mpi_hybrid(std::move(local_solver), MPI_COMM_WORLD);
            mpi_hybrid.set_time_integrator(chosen_integrator);

            initial_mass = mpi_hybrid.compute_total_mass();
            initial_energy = mpi_hybrid.compute_total_energy();

            double t = 0.0;
            int step = 0;
            int output_interval = 100;

            while (t < t_final) {
                if (g_signal_status == SIGINT) break;

                if (step % output_interval == 0) {
                    auto full_sol = mpi_hybrid.gather_solution();
                    if (mpi_rank == 0) {
                        std::string fname = output_dir + "/solution_" + std::to_string(step) + ".csv";
                        write_solution(fname, full_sol, t, domain_meta);
                        if (verbose) std::cout << "Step " << step << ", t = " << t << std::endl;
                    }
                }

                if (is_tty && !verbose) {
                    auto now = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed = now - start_time;
                    draw_progress_bar(step, total_steps, elapsed.count(), t, t_final);
                }

                mpi_hybrid.step(dt, advection_speed);
                t += dt;
                step++;
            }

            auto final_sol = mpi_hybrid.gather_solution();
            if (mpi_rank == 0) {
                write_solution(output_dir + "/solution_final.csv", final_sol, t_final, domain_meta);
            }

            final_mass = mpi_hybrid.compute_total_mass();
            final_energy = mpi_hybrid.compute_total_energy();
        } else
#endif
        {
            // Single-process Hybrid Execution with OpenMP multithreading
            auto fv = std::make_unique<FiniteVolumeSolver>();
            fv->initialize(x_start, x_interface, n_fv);
            fv->set_initial_condition(initial_condition);

            auto dg = std::make_unique<DiscontinuousGalerkinSolver>(p_order);
            dg->initialize(x_interface, x_end, n_dg);
            dg->set_initial_condition(initial_condition);

            HybridDomain hybrid(std::move(fv), std::move(dg));
            hybrid.set_time_integrator(chosen_integrator);

            initial_mass = hybrid.compute_total_mass();
            initial_energy = hybrid.compute_total_energy();

            double t = 0.0;
            int step = 0;
            int output_interval = 100;

            while (t < t_final) {
                if (g_signal_status == SIGINT) break;

                if (step % output_interval == 0) {
                    std::string fname = output_dir + "/solution_" + std::to_string(step) + ".csv";
                    write_solution(fname, hybrid.get_solution(), t, domain_meta);
                    if (verbose) std::cout << "Step " << step << ", t = " << t << std::endl;
                }

                if (is_tty && !verbose) {
                    auto now = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed = now - start_time;
                    draw_progress_bar(step, total_steps, elapsed.count(), t, t_final);
                }

                hybrid.step(dt, advection_speed);
                t += dt;
                step++;
            }

            write_solution(output_dir + "/solution_final.csv", hybrid.get_solution(), t_final, domain_meta);
            final_mass = hybrid.compute_total_mass();
            final_energy = hybrid.compute_total_energy();
        }

        if (is_tty && !verbose) {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - start_time;
            draw_progress_bar(total_steps, total_steps, elapsed.count(), t_final, t_final);
            std::cout << std::endl;
        }

        if (mpi_rank == 0) {
            auto end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> total_elapsed = end_time - start_time;
            double steps_per_sec = (total_elapsed.count() > 0) ? (total_steps / total_elapsed.count()) : 0.0;

            std::cout << Color::BoldGreen << "Simulation complete in ";
            if (total_elapsed.count() < 1.0) {
                 std::cout << std::setprecision(0) << (total_elapsed.count() * 1000.0) << "ms";
            } else {
                 std::cout << std::fixed << std::setprecision(2) << total_elapsed.count() << "s";
            }
            std::cout << " (" << std::setprecision(0) << steps_per_sec << " steps/s)." << Color::Reset << "\n\n";

            // Physics conservation diagnostics
            double mass_loss = (initial_mass > 0.0) ? std::abs(final_mass - initial_mass) / initial_mass * 100.0 : 0.0;
            double energy_retention = (initial_energy > 0.0) ? (final_energy / initial_energy) * 100.0 : 0.0;

            std::cout << Color::Bold << "Conservation Diagnostics:" << Color::Reset << "\n"
                      << "  " << Color::Blue << "Initial Mass:     " << Color::Reset << std::setprecision(8) << initial_mass << "\n"
                      << "  " << Color::Blue << "Final Mass:       " << Color::Reset << std::setprecision(8) << final_mass
                      << " (" << (mass_loss < 0.1 ? Color::Green : Color::Yellow) << "Rel Error: " << std::setprecision(4) << mass_loss << "%" << Color::Reset << ")\n"
                      << "  " << Color::Blue << "Kinetic Energy:   " << Color::Reset << std::setprecision(8) << initial_energy << " -> " << final_energy
                      << " (" << Color::Cyan << "Retention: " << std::setprecision(2) << energy_retention << "%" << Color::Reset << ")\n\n"
                      << "Results saved in " << Color::Bold << output_dir << "/" << Color::Reset << std::endl;

            std::string script_path = "scripts/plot_results.py";
            if (!fs::exists(script_path) && fs::exists("../scripts/plot_results.py")) {
                script_path = "../scripts/plot_results.py";
            }
            std::cout << "Visualize with: " << Color::Yellow << "python3 " + script_path << Color::Reset << std::endl;
        }

    } catch (const std::exception& e) {
        if (mpi_rank == 0) std::cerr << Color::BoldRed << "\nFATAL ERROR: " << e.what() << Color::Reset << std::endl;
#ifdef HYFLOW_ENABLE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

#ifdef HYFLOW_ENABLE_MPI
    MPI_Finalize();
#endif
    return 0;
}
