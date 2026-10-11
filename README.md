# HyFlow1D

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Build: CMake](https://img.shields.io/badge/Build-CMake-brightgreen.svg)](CMakeLists.txt)
[![Tests: CTest](https://img.shields.io/badge/Tests-5%2F5%20Passing-brightgreen.svg)](tests/)

**HyFlow1D** is a lightweight, high-performance research code implementing a coupled **Hybrid Finite Volume (FV) / Discontinuous Galerkin (DG)** solver for 1D aerodynamic flow problems and hyperbolic conservation laws.

The framework provides an extensible testbed for domain-decomposition coupling: applying robust, lower-order industrial finite volume methods near complex geometries, shocks, or boundary layers, while utilizing high-order discontinuous Galerkin schemes in smooth or far-field regions to preserve waveforms with minimal numerical dispersion and dissipation.

---

## Features

- **Hybrid Domain Coupling**: Bidirectional interface coupling linking a 1st-order Finite Volume scheme with a high-order Discontinuous Galerkin formulation.
- **Arbitrary Order DG**: Discontinuous Galerkin implementation supporting arbitrary polynomial degrees (default $P=3$).
- **Rigorous Numerics**: Gauss-Legendre quadrature integration with recurrence-optimized Legendre basis evaluations and analytical stiffness matrix computations.
- **Dependency-Free Core**: Written in clean, modern C++17 using only the standard library for core solvers.
- **Security-Hardened I/O**: Atomic file creation with `O_NOFOLLOW` symlink protections, path traversal sanitization, and restricted directory permissions (`0700`/`0600`).
- **Dual Visualization Ecosystem**:
  - Python CLI tool with color-coded terminal telemetry, inline sparklines, and matplotlib figures.
  - Interactive full-stack web dashboard built with Vite, React 19, Tailwind CSS, Recharts, and Three.js backed by Bun.
- **Comprehensive Test Suite**: Automated unit and security tests covering numerical precision, interface flux consistency, solver stability, and system security.

---

## Mathematical and Numerical Formulation

### 1. Governing Equation

HyFlow1D solves the 1D linear advection equation for a scalar conserved quantity $u(x, t)$:

$$\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} = 0, \quad x \in [0, 1], \; t \in [0, t_{\text{final}}]$$

with constant advection velocity $a = 1.0$ and periodic or transmissive boundary conditions. The default benchmark propagates a smooth Gaussian pulse:

$$u(x, 0) = \exp\left(-\frac{(x - x_0)^2}{2 \sigma^2}\right)$$

where $x_0 = 0.25$ and $\sigma = 0.05$.

### 2. Domain Decomposition

The computational domain $[0, 1]$ is partitioned at interface $x = 0.5$:
- **Near-field Subdomain** $\Omega_{\text{FV}} = [0, 0.5]$: Discretized with $N_{\text{FV}} = 50$ uniform finite volume cells ($\Delta x_{\text{FV}} = 0.01$).
- **Far-field Subdomain** $\Omega_{\text{DG}} = [0.5, 1.0]$: Discretized with $N_{\text{DG}} = 10$ elements ($\Delta x_{\text{DG}} = 0.05$) with polynomial degree $P = 3$.

### 3. Finite Volume Method (FV)

For cell $i \in [1, N_{\text{FV}}]$ with spatial extent $[x_{i-1/2}, x_{i+1/2}]$, the cell-averaged state is defined as:

$$\bar{u}_i(t) = \frac{1}{\Delta x_i} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x, t) \, dx$$

Integrating the conservation law yields the semi-discrete system:

$$\frac{d \bar{u}_i}{dt} = -\frac{1}{\Delta x_i} \left( F_{i+1/2} - F_{i-1/2} \right)$$

Using an upwind numerical flux ($a > 0$):

$$F_{i+1/2} = a \, \bar{u}_i$$

### 4. Discontinuous Galerkin Method (DG)

Within each element $e = [x_{e-1/2}, x_{e+1/2}]$, the solution is approximated by a degree-$P$ polynomial expansion mapped to the reference interval $\xi \in [-1, 1]$:

$$u_h^e(\xi, t) = \sum_{k=0}^P \hat{u}_k^e(t) P_k(\xi)$$

where $P_k(\xi)$ are orthogonal Legendre polynomials satisfying:

$$\int_{-1}^1 P_k(\xi) P_m(\xi) \, d\xi = \frac{2}{2k+1} \delta_{km}$$

Multiplying by test functions $P_k(\xi)$, integrating over the element, and applying integration by parts yields the weak formulation:

$$\frac{\Delta x_e}{2} \frac{2}{2k+1} \frac{d \hat{u}_k^e}{dt} - a \int_{-1}^1 u_h^e(\xi) \frac{d P_k}{d\xi} \, d\xi + \left[ \hat{F}(u_h^e) P_k(\xi) \right]_{-1}^{1} = 0$$

- **Mass Matrix**: Diagonal with entries $M_{kk} = \frac{\Delta x_e}{2k+1}$.
- **Stiffness Matrix**: $K_{km} = \int_{-1}^1 P_m(\xi) P_k'(\xi) \, d\xi$, which is strictly lower-triangular due to polynomial parity and orthogonality.
- **Interface Numerical Flux**: Upwind Riemann flux across element boundaries.

### 5. Hybrid Interface Coupling

At the coupling boundary $x = 0.5$:
1. The right boundary of the FV domain passes its cell average $\bar{u}_{N_{\text{FV}}}$ to the inflow face of the first DG element as an exterior Dirichlet state.
2. The left boundary of the first DG element evaluates its polynomial reconstruction $u_h^1(\xi = -1)$ to supply the numerical flux $F_{\text{interface}} = a \, \bar{u}_{N_{\text{FV}}}$ entering the DG domain, guaranteeing strict mass conservation across the hybrid boundary.
3. Both subdomains are advanced simultaneously in time using Forward Euler with a global time step restricted by the DG CFL stability criterion:

$$\Delta t = \text{CFL} \cdot \frac{\min\left(\Delta x_{\text{FV}}, \, \frac{\Delta x_{\text{DG}}}{2P + 1}\right)}{a}$$

---

## Directory Structure

```
HyFlow1D/
├── CMakeLists.txt              # CMake build configuration and hardening flags
├── include/                    # Header files
│   ├── dg_solver.hpp           # Discontinuous Galerkin solver interface
│   ├── fv_solver.hpp           # Finite Volume solver interface
│   ├── hybrid_coupling.hpp     # Hybrid domain coupling manager
│   ├── numerics.hpp            # Legendre polynomials & Gauss-Legendre quadrature
│   ├── solver.hpp              # Abstract solver base interface
│   └── term_colors.hpp         # Terminal ANSI color formatting utilities
├── src/                        # Implementation files
│   ├── dg_solver.cpp           # DG matrix assembly, quadrature, and kernels
│   ├── fv_solver.cpp           # FV flux calculations and time advancement
│   ├── hybrid_coupling.cpp     # Hybrid boundary exchange implementation
│   ├── main.cpp                # Simulation CLI driver and configuration
│   └── numerics.cpp            # Quadrature weights and polynomial evaluation
├── scripts/                    # Automation and visualization tools
│   └── plot_results.py         # Matplotlib visualization script with CLI telemetry
├── tests/                      # Verification and benchmarks
│   ├── benchmark_*.cpp         # Computational performance benchmarks
│   ├── test_numerics.cpp       # Quadrature and basis accuracy tests
│   ├── test_solvers.cpp        # Solver unit tests
│   ├── test_validation.cpp     # Interface conservation validation
│   ├── test_security.cpp       # File system and path sanitization tests
│   └── test_hybrid_security.cpp# Hybrid coupling exception and security tests
└── webapp/                     # Interactive web visualization dashboard
    ├── package.json            # Web app dependencies (Bun, React, Three.js)
    ├── server.ts               # Bun backend API server
    └── src/                    # Frontend React dashboard components
```

---

## Build and Installation

### Prerequisites

- **C++ Compiler**: GCC 7+, Clang 5+, or MSVC supporting C++17.
- **Build System**: CMake 3.10 or later.
- **Python** (optional): Python 3.7+ with `matplotlib` for CLI plots.
- **Bun** (optional): Bun runtime for running the web dashboard.

### Compilation

```bash
# Clone the repository
git clone https://github.com/dhruvhaldar/HyFlow1D.git
cd HyFlow1D

# Create and enter build directory
mkdir -p build && cd build

# Configure and compile
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The build process produces the `hyflow1d` executable and all test/benchmark targets within `build/`.

---

## Usage

### Command Line Interface

Run the solver directly from the command line:

```bash
./build/hyflow1d [OPTIONS]
```

#### Available Options

| Option | Long Option | Description |
| :--- | :--- | :--- |
| `-h` | `--help` | Display command usage and simulation details |
| `-v` | `--verbose` | Output step-by-step time log (disables animated progress bar) |
| `-o <dir>` | `--output <dir>` | Set output directory for CSV files (default: `output`) |
| `-c` | `--clean` | Remove previous `solution_*.csv` files before running |
| `-n` | `--dry-run` | Print simulation configuration and exit without running |

#### Execution Examples

```bash
# Run simulation with default parameters and clean previous results
./build/hyflow1d --clean

# Run simulation with custom output directory and verbose logging
./build/hyflow1d -o results_run1 -v

# Perform a dry run to verify configuration
./build/hyflow1d --dry-run
```

Simulation snapshots (`solution_*.csv`) and the final solution (`solution_final.csv`) will be saved in the designated output directory.

---

## Visualization

### 1. Python CLI Plotter

Use the included visualization script to plot snapshot profiles and inspect peak decay across the interface:

```bash
# Plot results from the default output directory
python3 scripts/plot_results.py

# Specify an alternate results directory and output filename
python3 scripts/plot_results.py results_run1 -o comparison.png

# Automatically launch the default system image viewer upon saving
python3 scripts/plot_results.py -p
```

The script displays terminal diagnostics, including:
- ASCII sparklines previewing spatial profiles.
- Min/max value tracking and peak decay metrics.
- Domain segmentation highlighting the FV and DG regions.

### 2. Full-Stack Web Dashboard

HyFlow1D includes a modern browser dashboard featuring interactive Recharts line visualizations and a Three.js 3D background.

```bash
cd webapp

# Install dependencies
bun install

# Start the Bun backend API server (listens on port 3001)
bun run server.ts &

# Start the Vite development server
bun run dev
```

Navigate to `http://localhost:5173` to trigger simulations and inspect results in real time. For more information, refer to [webapp/README.md](webapp/README.md).

---

## Verification and Testing

The repository includes a comprehensive CTest test suite covering numerical correctness, solver convergence, conservation properties, and security invariants.

To execute the test suite:

```bash
cd build
ctest --output-on-failure
```

### Test Targets

1. **`numerics_test`** (`tests/test_numerics.cpp`): Verifies Legendre polynomial orthogonality, recurrence relation evaluations, and Gauss-Legendre quadrature integration exactness.
2. **`solvers_test`** (`tests/test_solvers.cpp`): Validates individual time stepping and spatial discretization in FV and DG solvers.
3. **`validation_test`** (`tests/test_validation.cpp`): Confirms mass conservation across the hybrid coupling interface.
4. **`security_test`** (`tests/test_security.cpp`): Asserts path traversal prevention (`..`), symbolic link overwrite defense, and atomic permission masking.
5. **`hybrid_security_test`** (`tests/test_hybrid_security.cpp`): Checks input validation, null pointer handling, and parameter bounds checking in the hybrid coupling framework.

### Performance Benchmarks

Performance micro-benchmarks are built alongside the tests to measure computational throughput:

```bash
./build/benchmark_dg
./build/benchmark_fv
./build/benchmark_step
./build/benchmark_update
./build/benchmark_plot
```

---

## Security and Hardening

HyFlow1D implements security hardening suited for computational research environments:

- **Path Sanitization**: User-supplied paths are validated to block directory traversal attacks (`..`).
- **Atomic File Creation**: File writers use `open()` with `O_CREAT | O_EXCL | O_NOFOLLOW` on supported platforms to prevent Time-of-Check to Time-of-Use (TOCTOU) symlink vulnerabilities.
- **Least Privilege Permissions**: Output files and directories are restricted to owner-only permissions (`0600` for files, `0700` for directories).
- **Compiler Hardening**: Release builds incorporate `-fstack-protector-strong`, `_FORTIFY_SOURCE=2`, Position Independent Executable (`PIE`), and Full RELRO (`-Wl,-z,relro,-z,now`).

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
