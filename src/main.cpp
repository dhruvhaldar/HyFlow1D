#include <iostream>
#include <cstdio>
#include <fstream>
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
// Note: We use volatile sig_atomic_t for async-signal-safety.
// We avoid calling unsafe functions (like IO) in the handler directly.
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
    return std::exp(-(x - center) * (x - center) / (2 * width * width));
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

void write_solution(const std::string& filename, const std::vector<std::pair<double, double>>& solution, double time = -1.0) {
#if defined(HAS_UNISTD)
    // Security: Prevent Symlink Attack (Arbitrary File Overwrite) via O_NOFOLLOW
    // This atomic check prevents TOCTOU race conditions.
    // 0600 permissions are set atomically at creation.
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

    // Convert fd to FILE* for convenient formatting
    FILE* outfile = fdopen(fd, "w");
    if (!outfile) {
        std::cerr << "Error: fdopen failed." << std::endl;
        close(fd);
        return;
    }

    // Ensure strict closing via RAII
    struct FileDeleter {
        void operator()(FILE* f) const { std::fclose(f); }
    };
    std::unique_ptr<FILE, FileDeleter> file_guard(outfile);

    if (time >= 0.0) {
        fprintf(outfile, "# t=%.17g\n", time);
    }
    fprintf(outfile, "x,u\n");
    for (const auto& p : solution) {
        fprintf(outfile, "%.17g,%.17g\n", p.first, p.second);
    }
    // fclose handled by unique_ptr
#else
    // Windows/Other: Fallback to std::ofstream with basic checks (Best Effort)
    std::filesystem::path path(filename);
    if (std::filesystem::is_symlink(path)) {
        std::cerr << Color::BoldRed << "Error: Output file '" << filename << "' is a symbolic link." << Color::Reset << std::endl;
        std::cerr << "Refusing to overwrite symbolic links to prevent security risks." << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    if (time >= 0.0) {
        outfile << "# t=" << time << "\n";
    }
    outfile << "x,u\n";
    for (const auto& p : solution) {
        outfile << p.first << "," << p.second << "\n";
    }
    outfile.close();

    // Security: Restrict file permissions to owner only
    try {
        std::filesystem::permissions(filename,
            std::filesystem::perms::owner_read | std::filesystem::perms::owner_write,
            std::filesystem::perm_options::replace | std::filesystem::perm_options::nofollow);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << Color::Yellow << "Warning: Failed to set secure permissions on '" << filename << "': " << e.what() << Color::Reset << std::endl;
    }
#endif
}

void show_usage(const char* prog_name) {
    std::cout << "\n" << Color::BoldCyan << "  HyFlow1D  " << Color::Reset << " - Hybrid FV-DG Advection Solver\n\n"
              << "  " << Color::Bold << "Usage:" << Color::Reset << " " << prog_name << " [OPTIONS]\n\n"
              << "  " << Color::Bold << "Options:" << Color::Reset << "\n"
              << "    " << Color::Yellow << "-h, --help" << Color::Reset << "      Show this help message\n"
              << "    " << Color::Yellow << "-v, --verbose" << Color::Reset << "   Show detailed step output (disables progress bar)\n"
              << "    " << Color::Yellow << "-o, --output" << Color::Reset << "    Specify output directory (default: output)\n"
              << "    " << Color::Yellow << "-c, --clean" << Color::Reset << "     Clean output directory before running\n\n"
              << "  " << Color::Bold << "Simulation Details (Hardcoded):" << Color::Reset << "\n"
              << "    Domain: [0, 1] (FV: [0, 0.5], DG: [0.5, 1])\n"
              << "    Config: 50 FV cells, 10 DG elements (P=3)\n"
              << "    Physics: Linear Advection (a=1.0), Gaussian Pulse\n\n"
              << "  " << Color::Bold << "Examples:" << Color::Reset << "\n"
              << "    " << prog_name << "\n"
              << "    " << prog_name << " -o my_results -v\n";
}

// UX Helper: Calculate Levenshtein distance for "Did you mean?" suggestions
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
    const std::vector<std::string> valid_flags = {"--help", "-h", "--verbose", "-v", "--output", "-o", "--clean", "-c"};
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
    // Security: Set file creation mask to 0077 (rw-------)
    // This ensures all files created by this process are private by default,
    // eliminating the race condition between file creation and permission setting.
#if defined(_WIN32) || defined(_WIN64)
    _umask(0077);
#else
    umask(0077);
#endif

    // Register signal handler
    std::signal(SIGINT, signal_handler);

    try {
        namespace fs = std::filesystem;
        bool verbose = false;
        bool clean_output = false;
        std::string output_dir = "output";

        // Detect TTY to auto-configure UX. Default to true to be safe (fail-open).
        bool is_tty = true;
#if defined(HAS_UNISTD)
        if (!isatty(fileno(stdout))) {
            is_tty = false;
        }
#elif defined(HAS_IO_H)
        if (!_isatty(_fileno(stdout))) {
            is_tty = false;
        }
#endif
        // Auto-disable colors if not TTY
        if (!is_tty) {
            Color::enabled = false;
        } else {
            // Palette UX: Enable thousands separator for human-readable output
            // Only apply if is_tty is true (human is watching), otherwise keep raw numbers
            std::locale comma_locale(std::locale::classic(), new comma_numpunct());
            std::cout.imbue(comma_locale);
        }

        // Accessibility: Respect NO_COLOR standard (https://no-color.org/)
        if (std::getenv("NO_COLOR") != nullptr) {
            Color::enabled = false;
        }

        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                show_usage(argv[0]);
                return 0;
            } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
                verbose = true;
            } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
                if (i + 1 < argc) {
                    output_dir = argv[++i];
                    // Security: Validate output directory input
                    if (!is_safe_path(output_dir)) {
                        std::cerr << Color::BoldRed << "Error: Invalid output directory. Path traversal ('..') is not allowed." << Color::Reset << std::endl;
                        return 1;
                    }
                } else {
                    std::cerr << "Error: Output directory not specified after " << argv[i] << std::endl;
                    return 1;
                }
            } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--clean") == 0) {
                clean_output = true;
            } else {
                // Check if it looks like a flag
                if (argv[i][0] == '-') {
                    std::cerr << Color::BoldRed << "Error: Unknown option '" << argv[i] << "'" << Color::Reset << std::endl;
                    suggest_flag(argv[i]);
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

        // Create output directory
        // Security Check: Verify output_dir is not a file
        if (fs::exists(output_dir) && !fs::is_directory(output_dir)) {
            throw std::runtime_error("Output path '" + output_dir + "' exists and is not a directory.");
        }

        if (!fs::exists(output_dir)) {
            // Fix TOCTOU vulnerability: Only set permissions if we actually created the directory.
            // If create_directory returns false, it means the directory (or a symlink to it)
            // already existed, so we should NOT modify its permissions.
            if (fs::create_directory(output_dir)) {
                // Security: Restrict output directory permissions to owner only (0700)
                // to prevent unauthorized access to simulation results.
                // Added nofollow to prevent following symlinks if they are somehow introduced.
                fs::permissions(output_dir, fs::perms::owner_all, fs::perm_options::replace | fs::perm_options::nofollow);
            }
        } else {
            // Palette UX: Warn if output directory contains previous results
            // This prevents confusion when new results are mixed with old ones (e.g. fewer steps)
            bool has_solution_files = false;
            std::vector<fs::path> files_to_remove;

            for (const auto& entry : fs::directory_iterator(output_dir)) {
                if (entry.is_regular_file()) {
                    std::string fname = entry.path().filename().string();
                    if (fname.rfind("solution_", 0) == 0 && entry.path().extension() == ".csv") {
                        has_solution_files = true;
                        if (clean_output) {
                            files_to_remove.push_back(entry.path());
                        }
                    }
                }
            }

            if (clean_output && has_solution_files) {
                // Perform cleanup
                for (const auto& p : files_to_remove) {
                    fs::remove(p);
                }
                std::cout << Color::Blue << "ℹ️  Cleaned " << files_to_remove.size() << " file(s) from '" << output_dir << "'." << Color::Reset << std::endl;
            } else if (has_solution_files) {
                 std::cout << Color::Yellow << "⚠️  Warning: Output directory '" << output_dir << "' contains existing solution files.\n"
                           << "    New results may mix with old ones. Use " << Color::Bold << "--clean" << Color::Reset << Color::Yellow << " to remove them." << Color::Reset << "\n" << std::endl;
            }
        }

        // Print Start Message
        std::cout << Color::BoldCyan << "Starting Hybrid FV-DG Simulation..." << Color::Reset << std::endl;

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

        // Print Configuration Summary
        std::cout << "\n" << Color::Bold << "Simulation Configuration:" << Color::Reset << "\n"
                  << "  " << Color::Blue << "Domain:      " << Color::Reset
                  << "[" << std::fixed << std::setprecision(2) << x_start << "] |---FV---| ["
                  << x_interface << "] |---DG---| [" << x_end << "]\n"
                  << "  " << Color::Blue << "FV Cells:    " << Color::Reset << Color::Bold << n_fv << Color::Reset << "\n"
                  << "  " << Color::Blue << "DG Elements: " << Color::Reset << Color::Bold << n_dg << Color::Reset << " (Order P=" << Color::Bold << p_order << Color::Reset << ")\n"
                  << "  " << Color::Blue << "Time:        " << Color::Reset << Color::Bold << "0.0 -> " << t_final << Color::Reset << "\n"
                  << "  " << Color::Blue << "CFL:         " << Color::Reset << Color::Bold << cfl << Color::Reset << "\n"
                  << "  " << Color::Yellow << "dt:          " << Color::Reset << Color::Bold << std::defaultfloat << std::setprecision(6) << dt << Color::Reset << "\n"
                  << "  " << Color::Yellow << "Total Steps: " << Color::Reset << Color::Bold << total_steps << Color::Reset << "\n"
                  << std::endl;

        auto start_time = std::chrono::steady_clock::now();

        while (t < t_final) {
            if (g_signal_status == SIGINT) {
                std::cout << Color::Reset << "\n";
                std::cerr << Color::Yellow << "\n[!] Simulation paused by user." << Color::Reset << std::endl;
                std::cerr << "    Stopped at step " << step << " (t=" << t << ")" << std::endl;
                std::cerr << "    Partial results saved in " << Color::Bold << output_dir << "/" << Color::Reset << std::endl;

                // Palette UX: Suggest visualization for partial results
                std::string script_path = "scripts/plot_results.py";
                if (!fs::exists(script_path) && fs::exists("../scripts/plot_results.py")) {
                    script_path = "../scripts/plot_results.py";
                }
                std::string viz_cmd = "python3 " + script_path;
                if (output_dir != "output") viz_cmd += " " + output_dir;

                std::cerr << "Visualize with: " << Color::Yellow << viz_cmd << Color::Reset << std::endl;

                // Break loop to allow RAII cleanup (e.g. file buffers, solvers)
                return 130;
            }

            if (step % output_interval == 0) {
                std::string fname = output_dir + "/solution_" + std::to_string(step) + ".csv";
                write_solution(fname, hybrid.get_solution(), t);
                if (verbose) {
                    std::cout << "Step " << step << ", t = " << t << std::endl;
                }
            }

            if (!verbose && is_tty) {
                auto now = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed = now - start_time;
                draw_progress_bar(step, total_steps, elapsed.count(), t, t_final);
            } else if (!verbose && !is_tty && step % std::max(1, total_steps / 10) == 0) {
                 // Simple progress for non-TTY (every 10%)
                 std::cout << "Progress: " << int((double)step/total_steps * 100)
                           << "% (t=" << std::fixed << std::setprecision(2) << t << "/" << t_final << ")" << std::endl;
            }

            hybrid.step(dt, advection_speed);
            t += dt;
            step++;
        }

        if (!verbose && is_tty) {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - start_time;
            draw_progress_bar(total_steps, total_steps, elapsed.count(), t_final, t_final);
            std::cout << std::endl;
        }

        // Final output
        write_solution(output_dir + "/solution_final.csv", hybrid.get_solution(), t_final);

        auto end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> total_elapsed = end_time - start_time;
        double steps_per_sec = (total_elapsed.count() > 0) ? (total_steps / total_elapsed.count()) : 0.0;

        std::cout << Color::BoldGreen << "✔ Simulation Complete in ";
        if (total_elapsed.count() < 1.0) {
             std::cout << std::setprecision(0) << (total_elapsed.count() * 1000.0) << "ms";
        } else {
             std::cout << std::fixed << std::setprecision(2) << total_elapsed.count() << "s";
        }
        std::cout << " (" << std::setprecision(0) << steps_per_sec << " steps/s)." << Color::Reset << "\n"
                  << "Results saved in " << Color::Bold << output_dir << "/" << Color::Reset << std::endl;

        // Palette UX: Suggest next step with smart path detection
        std::string script_path = "scripts/plot_results.py";
        if (!fs::exists(script_path)) {
            if (fs::exists("../scripts/plot_results.py")) {
                script_path = "../scripts/plot_results.py";
            }
        }

        std::string viz_cmd = "python3 " + script_path;
        if (output_dir != "output") {
            viz_cmd += " " + output_dir;
        }
        std::cout << "Visualize with: " << Color::Yellow << viz_cmd << Color::Reset << std::endl;

    } catch (const std::exception& e) {
        std::cerr << Color::BoldRed << "\nFATAL ERROR: " << e.what() << Color::Reset << std::endl;
        return 1;
    } catch (...) {
        std::cerr << Color::BoldRed << "\nFATAL ERROR: Unknown exception occurred." << Color::Reset << std::endl;
        return 1;
    }

    return 0;
}
