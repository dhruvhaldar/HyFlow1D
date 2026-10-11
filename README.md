# HyFlow1D

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Build: CMake](https://img.shields.io/badge/Build-CMake-brightgreen.svg)](CMakeLists.txt)
[![Parallelism: MPI + OpenMP](https://img.shields.io/badge/Parallelism-MPI%20%2B%20OpenMP-orange.svg)](src/main.cpp)
[![Tests: CTest](https://img.shields.io/badge/Tests-10%2F10%20Passing-brightgreen.svg)](tests/)

**HyFlow1D** is a lightweight, high-performance research Proof of Concept (PoC) code implementing a coupled **Hybrid Finite Volume (FV) / Discontinuous Galerkin (DG)** solver for 1D aerodynamic flows and hyperbolic conservation laws.

As a research testbed, HyFlow1D demonstrates domain-decomposition coupling: deploying robust, lower-order industrial finite volume methods near complex geometries, shocks, or boundary layers, while utilizing high-order discontinuous Galerkin schemes in smooth or far-field regions to preserve waveforms with minimal numerical dispersion and dissipation. The code demonstrates hybrid distributed/shared-memory parallelism combining **MPI** and **OpenMP**.

---

## Features

- **Hybrid Domain Coupling**: Bidirectional interface coupling linking a 1st-order Finite Volume scheme with a high-order Discontinuous Galerkin formulation.
- **Governing Equations & Numerical Fluxes**:
  - Linear advection ($\partial_t u + a \partial_x u = 0$) and nonlinear Burgers' equation ($\partial_t u + \partial_x (\frac{1}{2} u^2) = 0$).
  - Monotone Local Lax-Friedrichs (**Rusanov**) numerical interface flux handling shocks, rarefactions, and bidirectional waves.
- **Persson & Peraire Modal Shock Sensor**:
  - Modal energy decay sensor $S_e = \log_{10} \left( \frac{\hat{u}_P^2 / (2P+1)}{\sum \hat{u}_k^2 / (2k+1)} + \epsilon \right)$ to detect shocks and steep gradient elements.
- **Hybrid HPC Parallelism (MPI + OpenMP)**:
  - Distributed-memory domain decomposition using non-blocking MPI halo communication (`MPI_Isend`, `MPI_Irecv`, `MPI_Waitall`).
  - Shared-memory multi-threading within subdomains using OpenMP parallel loops and reductions.
- **Higher-Order Time Integration**: Multi-stage 3rd-order Strong Stability Preserving Runge-Kutta (**SSPRK3**) integrator alongside 1st-order Forward Euler.
- **Arbitrary Order DG**: Discontinuous Galerkin implementation supporting arbitrary polynomial degrees (default $P = 3$).
- **Conservation & Invariant Diagnostics**: Real-time analytical tracking of total conserved mass ($\int u \, dx$) and kinetic energy ($\frac{1}{2}\int u^2 \, dx$) across the hybrid boundary.
- **ParaView VTK Export**: Built-in ASCII `.vtu` (VTK XML UnstructuredGrid) export for 2D/3D ParaView post-processing.
- **Rigorous Numerics**: Gauss-Legendre quadrature integration with recurrence-optimized Legendre basis evaluations and analytical stiffness matrix computations.
- **Security-Hardened I/O**: Atomic file creation with `O_NOFOLLOW` symlink protections, path traversal sanitization, and restricted directory permissions (`0700`/`0600`).
- **Dual Visualization Ecosystem**:
  - Python CLI tool with color-coded terminal telemetry, inline sparklines, and matplotlib figures.
  - Interactive full-stack web dashboard built with Vite, React 19, Tailwind CSS, Recharts, and Three.js backed by Node.js or Bun.
- **Automated Test Suite**: 10 unit, integration, physics, and security test targets covering numerical integration, solver convergence, Burgers equation, shock sensor decay, OpenMP thread consistency, and distributed MPI communication.

---

## Architecture Details

HyFlow1D is structured around clean object-oriented abstraction, modular domain decomposition, and multi-tier parallelism suited for scientific computing research.

```
                              hyflow1d Application Driver
                                      (main.cpp)
                                           │
             ┌─────────────────────────────┴─────────────────────────────┐
             ▼                                                           ▼
     HybridDomain (Single-Node)                            MpiHybridDomain (Distributed)
     ├── Left Subdomain (FV)                               ├── Rank 0: Subdomain 0 (FV)
     └── Right Subdomain (DG)                              └── Rank 1: Subdomain 1 (DG)
             │                                                           │
             │   Local Pointer Exchange                                  │   Non-Blocking MPI Halo
             │   (Shared Memory)                                         │   (MPI_Isend / MPI_Irecv)
             ▼                                                           ▼
┌─────────────────────────┐                                 ┌─────────────────────────┐
│   FiniteVolumeSolver    │                                 │ DiscontinuousGalerkin   │
│   (P=0 Piecewise Const) │                                 │ (P=3 Modal Legendre)    │
│   #pragma omp parallel  │                                 │ #pragma omp parallel    │
└─────────────────────────┘                                 └─────────────────────────┘
             │                                                           │
             └─────────────────────────────┬─────────────────────────────┘
                                           ▼
                                 Numerics Engine (numerics.cpp)
                                 ├── Gauss-Legendre Quadrature
                                 └── Legendre Polynomial Recurrence
```

### 1. Class Hierarchy and Core Abstractions

- **`Solver1D` (`include/solver.hpp`)**: Abstract base class defining the contract for any 1D subdomain solver.
  - Lifecycle: `initialize(start, end, n_elements)`, `set_initial_condition(func)`.
  - Discretization & Time-stepping: `compute_rhs(t, a)`, `update_state(dt)`, `save_stage()`, `blend_stage(alpha, beta, dt)`.
  - Coupling Interface: `get_left_boundary_value()`, `get_right_boundary_value()`, `set_left_neighbor_value(val)`, `set_right_neighbor_value(val)`.
  - Diagnostic Reductions: `compute_total_mass()`, `compute_total_energy()`.
- **`FiniteVolumeSolver` (`include/fv_solver.hpp`)**: Standard cell-averaged finite volume discretization with 1st-order upwind numerical fluxes and OpenMP loop-level parallelism.
- **`DiscontinuousGalerkinSolver` (`include/dg_solver.hpp`)**: Modal Discontinuous Galerkin implementation utilizing Legendre polynomial basis functions mapped to reference element $\xi \in [-1, 1]$.
- **`HybridDomain` (`include/hybrid_coupling.hpp`)**: Coordinates boundary flux matching and multi-stage time stepping between two subdomains in a single process.
- **`MpiHybridDomain` (`include/mpi_hybrid_coupling.hpp`)**: Manages domain decomposition across distributed MPI ranks with non-blocking halo exchanges and distributed reduction collectives.

### 2. Multi-Level Parallelism Model

```
   Cluster / Node Level (Distributed Memory - MPI):
   Rank 0 [0.0, 0.5] (FV)                       Rank 1 [0.5, 1.0] (DG)
   ┌────────────────────────────────┐           ┌────────────────────────────────┐
   │ Subdomain 0 Elements           │  MPI Halo │ Subdomain 1 Elements           │
   │ [Cell 0 ... Cell 49]           │◄─────────►│ [Elem 0 ... Elem 9]            │
   └────────────────────────────────┘  Exchange └────────────────────────────────┘
                 │                                             │
                 ▼                                             ▼
   Core / Socket Level (Shared Memory - OpenMP):
   Thread 0  Thread 1  ... Thread T             Thread 0  Thread 1  ... Thread T
   ┌─────────┬─────────┬────────────┐           ┌─────────┬─────────┬────────────┐
   │ Cells   │ Cells   │ Cells      │           │ Elems   │ Elems   │ Elems      │
   │ 0..15   │ 16..31  │ 32..49     │           │ 0..2    │ 3..5    │ 6..9       │
   └─────────┴─────────┴────────────┘           └─────────┴─────────┴────────────┘
```

1. **Distributed Memory (MPI)**:
   - Rank 0 owns the Finite Volume subdomain ($[0.0, 0.5]$).
   - Rank 1 owns the Discontinuous Galerkin subdomain ($[0.5, 1.0]$).
   - Boundary exchange operates through non-blocking point-to-point primitives (`MPI_Isend` and `MPI_Irecv` with distinct communication tags for left-to-right and right-to-left flux propagation), finalized via `MPI_Waitall`.
   - Global diagnostic reductions (`compute_total_mass`, `compute_total_energy`) utilize `MPI_Allreduce` with `MPI_SUM`.
   - Snapshot gathering uses coordinated MPI messages onto Rank 0 for consolidated CSV generation.
2. **Shared Memory (OpenMP)**:
   - Independent element loops in `compute_rhs`, `blend_stage`, and `update_state` are parallelized across available hardware threads with static scheduling.
   - Diagnostic metrics employ OpenMP reduction clauses (`reduction(+:total)`), preventing race conditions without lock contention.

### 3. Multi-Stage SSPRK3 Time Integration

To preserve stability and accuracy across the high-order DG subdomain, HyFlow1D implements the 3rd-order Strong Stability Preserving Runge-Kutta scheme (Shu-Osher):

$$\mathbf{u}^{(1)} = \mathbf{u}^n + \Delta t \, \mathcal{L}(\mathbf{u}^n)$$

$$\mathbf{u}^{(2)} = \frac{3}{4}\mathbf{u}^n + \frac{1}{4}\left(\mathbf{u}^{(1)} + \Delta t \, \mathcal{L}(\mathbf{u}^{(1)})\right)$$

$$\mathbf{u}^{n+1} = \frac{1}{3}\mathbf{u}^n + \frac{2}{3}\left(\mathbf{u}^{(2)} + \Delta t \, \mathcal{L}(\mathbf{u}^{(2)})\right)$$

Between every intermediate stage, `exchange_boundaries()` is evaluated so interface flux boundary conditions remain fully synchronized at all fractional time levels.

---

## Mathematical and Numerical Formulation

### 1. Governing Equations & Numerical Fluxes

HyFlow1D supports both linear advection and nonlinear conservation laws for scalar quantity $u(x, t)$:

$$\frac{\partial u}{\partial t} + \frac{\partial f(u)}{\partial x} = 0, \quad x \in [0, 1], \; t \in [0, t_{\text{final}}]$$

1. **Linear Advection**: $f(u) = a u$ with constant velocity $a = 1.0$. The initial condition is a smooth Gaussian pulse:

$$u(x, 0) = \exp\left(-\frac{(x - x_0)^2}{2 \sigma^2}\right), \quad x_0 = 0.25, \; \sigma = 0.05$$

2. **Nonlinear Burgers' Equation**: $f(u) = \frac{1}{2} u^2$. Represents nonlinear convective steepening and shock formation. The benchmark profile is a compression wave:

$$u(x, 0) = 0.5 + 0.5 \sin(2\pi x)$$

Interfaces between elements (and subdomains) employ the monotone **Local Lax-Friedrichs (Rusanov)** numerical flux:

$$F(u_L, u_R) = \frac{1}{2}\left(f(u_L) + f(u_R)\right) - \frac{1}{2}\lambda_{\max}(u_R - u_L)$$

where $\lambda_{\max} = \max\left(|f'(u_L)|, \, |f'(u_R)|\right)$ is the maximum local characteristic wave speed. For linear advection with $a > 0$, Rusanov reduces identically to the standard upwind flux $a u_L$.

### 2. Domain Decomposition

The computational domain $[0, 1]$ is partitioned at interface $x = 0.5$:
- **Near-field Subdomain** $\Omega_{\text{FV}} = [0, 0.5]$: Discretized with $N_{\text{FV}} = 50$ uniform finite volume cells ($\Delta x_{\text{FV}} = 0.01$).
- **Far-field Subdomain** $\Omega_{\text{DG}} = [0.5, 1.0]$: Discretized with $N_{\text{DG}} = 10$ elements ($\Delta x_{\text{DG}} = 0.05$) with polynomial degree $P = 3$.

### 3. Finite Volume Method (FV)

For cell $i \in [1, N_{\text{FV}}]$ with spatial extent $[x_{i-1/2}, x_{i+1/2}]$, the cell-averaged state is defined as:

$$\bar{u}_i(t) = \frac{1}{\Delta x_i} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x, t) \, dx$$

Integrating the conservation law yields the conservative semi-discrete system:

$$\frac{d \bar{u}_i}{dt} = -\frac{1}{\Delta x_i} \left( F_{i+1/2} - F_{i-1/2} \right)$$

### 4. Discontinuous Galerkin Method (DG)

Within each element $e = [x_{e-1/2}, x_{e+1/2}]$, the solution is approximated by a polynomial expansion of degree $P$ mapped to the reference interval $\xi \in [-1, 1]$:

$$u_h^e(\xi, t) = \sum_{k=0}^P \hat{u}_k^e(t) P_k(\xi)$$

where $P_k(\xi)$ are orthogonal Legendre polynomials satisfying:

$$\int_{-1}^1 P_k(\xi) P_m(\xi) \, d\xi = \frac{2}{2k+1} \delta_{km}$$

Multiplying by test functions $P_k(\xi)$, integrating over the element, and applying integration by parts yields the weak formulation:

$$\frac{\Delta x_e}{2k+1} \frac{d \hat{u}_k^e}{dt} = \int_{-1}^1 f(u_h^e(\xi)) \frac{d P_k}{d\xi} \, d\xi - \left( \hat{F}_{e+1/2} - (-1)^k \hat{F}_{e-1/2} \right)$$

- **Mass Matrix**: Diagonal with entries $M_{kk} = \frac{\Delta x_e}{2k+1}$.
- **Volume Integrals**: For linear advection, evaluated using precomputed analytical stiffness matrices. For nonlinear Burgers' equation, evaluated using Gauss-Legendre quadrature ($N = P + 2$ points) exact for polynomials up to degree $2N-1$.
- **Mode $k=0$ Conservation**: Since $P_0'(\xi) = 0$, the volume integral vanishes for the zeroth mode, ensuring exact cell-average conservation identically mirroring Finite Volume schemes.

### 5. Persson & Peraire Modal Shock Sensor

To detect shocks, steep fronts, or under-resolved gradient regions without empirical tuning, HyFlow1D implements the Persson & Peraire (2002) modal decay sensor. In element $e$, the relative energy in the highest Legendre polynomial mode $P$ is:

$$S_e = \log_{10} \left( \frac{\frac{\hat{u}_P^2}{2P+1}}{\sum_{k=0}^P \frac{\hat{u}_k^2}{2k+1}} + 10^{-16} \right)$$

- **Smooth Waveforms**: Legendre coefficients decay exponentially ($\hat{u}_k \sim \mathcal{O}(e^{-\sigma k})$), yielding $S_e \ll -4.0$.
- **Shocks / Discontinuities**: Decay slows to $\mathcal{O}(1/k)$, causing $S_e$ to jump to $\approx -1.5$ to $-3.0$.
- Elements with $S_e > S_0$ (default threshold $S_0 = -3.0$) are flagged as troubled elements.

### 6. Hybrid Interface Coupling

At the coupling boundary $x = 0.5$:
1. The right boundary of the FV domain passes its cell average $\bar{u}_{N_{\text{FV}}}$ to the inflow face of the first DG element.
2. The left boundary of the first DG element evaluates its polynomial reconstruction $u_h^1(\xi = -1)$ to supply boundary state to the Rusanov interface flux $F_{\text{interface}} = F(\bar{u}_{N_{\text{FV}}}, u_h^1(-1))$, guaranteeing strict mass conservation across the hybrid boundary.
3. Both subdomains are advanced simultaneously in time using SSPRK3 with a global time step restricted by the DG CFL stability criterion:

$$\Delta t = \text{CFL} \cdot \frac{\min\left(\Delta x_{\text{FV}}, \, \frac{\Delta x_{\text{DG}}}{2P + 1}\right)}{\lambda_{\max}}$$

---

## Directory Structure

```
HyFlow1D/
├── CMakeLists.txt              # CMake build configuration and hardening flags
├── include/                    # Header files
│   ├── dg_solver.hpp           # Discontinuous Galerkin solver interface
│   ├── fv_solver.hpp           # Finite Volume solver interface
│   ├── hybrid_coupling.hpp     # Hybrid domain coupling manager (shared-memory)
│   ├── mpi_hybrid_coupling.hpp # Distributed MPI hybrid domain manager
│   ├── numerics.hpp            # Legendre polynomials & Gauss-Legendre quadrature
│   ├── solver.hpp              # Abstract solver base interface
│   └── term_colors.hpp         # Terminal ANSI color formatting utilities
├── src/                        # Implementation files
│   ├── dg_solver.cpp           # DG matrix assembly, quadrature, and kernels
│   ├── fv_solver.cpp           # FV flux calculations and time advancement
│   ├── hybrid_coupling.cpp     # Shared-memory hybrid boundary exchange
│   ├── mpi_hybrid_coupling.cpp # Distributed MPI non-blocking halo exchanges
│   ├── main.cpp                # Simulation CLI driver, MPI/OpenMP runtime
│   └── numerics.cpp            # Quadrature weights and polynomial evaluation
├── scripts/                    # Automation and visualization tools
│   └── plot_results.py         # Matplotlib visualization script with CLI telemetry
├── tests/                      # Verification and benchmarks
│   ├── benchmark_*.cpp         # Computational performance benchmarks
│   ├── test_numerics.cpp       # Quadrature and basis accuracy tests
│   ├── test_solvers.cpp        # Solver unit tests
│   ├── test_validation.cpp     # Interface conservation validation
│   ├── test_security.cpp       # File system and path sanitization tests
│   ├── test_hybrid_security.cpp# Hybrid coupling exception and security tests
│   ├── test_ssprk3.cpp         # Multi-stage SSPRK3 time integration test
│   ├── test_openmp.cpp         # OpenMP thread equivalence test
│   ├── test_mpi.cpp            # Distributed MPI halo communication test
│   ├── test_burgers.cpp        # Nonlinear Burgers' equation & Rusanov flux test
│   └── test_sensor.cpp         # Persson & Peraire modal shock sensor test
└── webapp/                     # Interactive web visualization dashboard
    ├── package.json            # Web app dependencies (React, Three.js, Recharts)
    ├── server.ts               # Universal Node.js / Bun backend API server
    └── src/                    # Frontend React dashboard components
```

---

## Build and Installation

### Prerequisites

- **C++ Compiler**: GCC 7+, Clang 5+, or MSVC supporting C++17.
- **Build System**: CMake 3.10 or later.
- **OpenMP** (optional): OpenMP 4.5+ for multi-threading.
- **MPI** (optional): OpenMPI or MPICH for distributed-memory execution.
- **Python** (optional): Python 3.7+ with `matplotlib` for CLI plots.
- **Node.js / Bun** (optional): Node.js (>= 20) or Bun for the web dashboard.

### Compilation

```bash
# Clone the repository
git clone https://github.com/dhruvhaldar/HyFlow1D.git
cd HyFlow1D

# Create and enter build directory
mkdir -p build && cd build

# Configure and compile with OpenMP and MPI enabled
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=ON -DENABLE_MPI=ON ..
make -j$(nproc)
```

The build process produces the `hyflow1d` executable and all test targets within `build/`.

---

## Usage

### Command Line Interface

```bash
./build/hyflow1d [OPTIONS]
```

#### Available Options

| Option | Long Option | Description |
| :--- | :--- | :--- |
| `-h` | `--help` | Display command usage and simulation details |
| `-v` | `--verbose` | Output step-by-step time log (disables animated progress bar) |
| `-o <dir>` | `--output <dir>` | Set output directory for CSV / VTK files (default: `output`) |
| `-c` | `--clean` | Remove previous `solution_*.csv` files before running |
| `-n` | `--dry-run` | Print simulation configuration and exit without running |
| `-t <N>` | `--threads <N>` | Set number of OpenMP worker threads |
| `-i <name>`| `--integrator <name>`| Select time integrator (`ssprk3` or `euler`, default: `ssprk3`) |
| `-e <name>`| `--equation <name>`  | Select equation (`advection` or `burgers`, default: `advection`) |
| | `--vtk` | Export solution to VTK XML (`.vtu`) for ParaView visualization |
| | `--sensor-threshold <T>` | Set modal shock sensor threshold (default: `-3.0`) |

#### Execution Examples

```bash
# Linear advection with SSPRK3 and 4 OpenMP threads
./build/hyflow1d --clean -t 4 --integrator ssprk3

# Nonlinear Burgers' equation with shock capturing and ParaView VTK export
./build/hyflow1d --clean -e burgers -i ssprk3 --vtk

# Distributed execution across 2 MPI ranks (Rank 0: FV, Rank 1: DG)
mpirun -n 2 ./build/hyflow1d --clean -e burgers --vtk

# Hybrid MPI + OpenMP: 2 MPI ranks with 4 threads each
OMP_NUM_THREADS=4 mpirun -n 2 ./build/hyflow1d --clean -e burgers

# Perform a dry run to verify configuration
./build/hyflow1d --dry-run -e burgers
```

---

## Visualization

### 1. Python CLI Plotter

```bash
# Plot results from the default output directory
python3 scripts/plot_results.py

# Specify custom output image filename and open preview
python3 scripts/plot_results.py -o comparison.png -p
```

The script displays terminal diagnostics, including:
- ASCII sparklines previewing spatial wave profiles.
- Min/max value tracking and wave peak decay analysis.
- Domain segmentation highlighting the FV and DG subdomains.

### 2. ParaView VTK Visualization

Passing the `--vtk` flag generates `output/solution_final.vtu` (VTK XML UnstructuredGrid):

1. Open **ParaView** (`paraview output/solution_final.vtu`).
2. Click **Apply**.
3. Choose the `u` scalar field to visualize 1D wave profiles, shock fronts, and element reconstructions.

### 3. Full-Stack Web Dashboard

```bash
cd webapp
pnpm install # or: npm install / bun install
node server.ts & # or: bun run server.ts &
pnpm run dev # or: npm run dev / bun run dev
```

Navigate to `http://localhost:5173` to configure equation types, adjust thread counts, run simulations, and inspect results with real-time shock sensor metrics. Refer to [webapp/README.md](webapp/README.md) for full instructions.

---

## Verification and Testing

Execute the test suite using CTest:

```bash
cd build
ctest --output-on-failure
```

### Test Suite Overview (10 Test Targets)

1. **`numerics_test`** (`tests/test_numerics.cpp`): Legendre orthogonality and Gauss-Legendre quadrature exactness.
2. **`solvers_test`** (`tests/test_solvers.cpp`): Individual FV and DG solver discretization verification.
3. **`validation_test`** (`tests/test_validation.cpp`): Mass conservation across the hybrid coupling boundary.
4. **`security_test`** (`tests/test_security.cpp`): Path traversal (`..`) and symlink defenses.
5. **`hybrid_security_test`** (`tests/test_hybrid_security.cpp`): Null pointer and domain gap/overlap validation.
6. **`ssprk3_test`** (`tests/test_ssprk3.cpp`): 3rd-order SSP Runge-Kutta time accuracy and conservation.
7. **`openmp_test`** (`tests/test_openmp.cpp`): Verification of identical numerical results between 1 and $N$ OpenMP threads.
8. **`mpi_test`** (`tests/test_mpi.cpp`): Distributed execution verifying non-blocking halo exchanges and distributed reductions across 2 MPI ranks.
9. **`burgers_test`** (`tests/test_burgers.cpp`): Nonlinear Burgers' equation flux conservation and multi-stage time integration.
10. **`sensor_test`** (`tests/test_sensor.cpp`): Persson & Peraire modal shock sensor verification distinguishing smooth waveforms from steep fronts.

### Performance Benchmarks

```bash
./build/benchmark_dg
./build/benchmark_fv
./build/benchmark_step
./build/benchmark_update
./build/benchmark_plot
```

---

## Security and Hardening

- **Path Sanitization**: User-supplied paths are validated to block directory traversal attacks (`..`).
- **Atomic File Creation**: File writers use `open()` with `O_CREAT | O_EXCL | O_NOFOLLOW` on supported platforms to prevent Time-of-Check to Time-of-Use (TOCTOU) symlink vulnerabilities.
- **Least Privilege Permissions**: Output files and directories are restricted to owner-only permissions (`0600` for files, `0700` for directories).
- **Compiler Hardening**: Release builds incorporate `-fstack-protector-strong`, `_FORTIFY_SOURCE=2`, Position Independent Executable (`PIE`), and Full RELRO (`-Wl,-z,relro,-z,now`).

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
