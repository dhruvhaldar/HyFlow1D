# HyFlow1D

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Language: C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build: CMake](https://img.shields.io/badge/Build-CMake-green.svg)

**HyFlow1D** is a lightweight, high-performance research code implementing a **Hybrid Finite Volume (FV) / Discontinuous Galerkin (DG)** solver for 1D aerodynamic flows.

It serves as a Proof of Concept (PoC) for coupling robust low-order industrial schemes with high-order accurate methods to achieve superior fidelity in transient simulations.

## 🚀 Features

*   **Hybrid Coupling**: Seamlessly couples 1st-order FV (robust near-field) with High-Order DG (accurate far-field).
*   **High-Order Accuracy**: Arbitrary polynomial order DG implementation (default $P=3$).
*   **Efficient Numerics**: Uses Gauss-Legendre quadrature and optimized Legendre basis evaluation.
*   **Simple & Modular**: Written in clean C++17 with no external dependencies (standard library only).
*   **Visualization**: Python scripts included for immediate result analysis.

## 🛠️ Build & Install

### Prerequisites
*   **C++ Compiler**: GCC 7+, Clang 5+, or MSVC (C++17 support required).
*   **CMake**: Version 3.10 or higher.
*   **Python**: For plotting (optional, requires `matplotlib` and `pandas`).

### Compilation

```bash
# Clone the repository
git clone https://github.com/username/HyFlow1D.git
cd HyFlow1D

# Create build directory
mkdir build && cd build

# Configure and Build
cmake ..
make
```

## 🏃 Usage

### Running the Solver
After building, run the executable from the `build` directory:

```bash
./hyflow1d
```

The simulation propagates a Gaussian pulse through the domain. The domain $[0, 0.5]$ uses Finite Volume, and $[0.5, 1.0]$ uses Discontinuous Galerkin.

Output files (`solution_*.csv`) will be generated in the current directory.

### Visualizing Results
To see the wave move through the hybrid interface:

```bash
python3 ../scripts/plot_results.py
```
*(Ensure `matplotlib` and `pandas` are installed: `pip install matplotlib pandas`)*

## 🧪 Testing

The project includes unit tests for numerical integration and solver logic.

```bash
cd build
make test
# Output:
#     Start 1: numerics_test
# 1/2 Test #1: numerics_test ....................   Passed    0.00 sec
#     Start 2: solvers_test
# 2/2 Test #2: solvers_test .....................   Passed    0.00 sec
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Theory

**The Hybrid Strategy:**
1.  **Near-Field (FV)**: Uses standard Finite Volume methods ($P=0$). This handles complex geometries and discontinuities robustly but is dissipative.
2.  **Far-Field (DG)**: Uses Discontinuous Galerkin methods ($P=k$). This preserves wave shape and energy over long distances (low dispersion/dissipation).
3.  **Interface**: Fluxes are exchanged by treating the neighbor's state as a boundary condition. Reconstruction is performed to pass high-order information to the FV cells and vice-versa.

---
*Developed for research on Next-Generation Aerodynamic Simulation.*
