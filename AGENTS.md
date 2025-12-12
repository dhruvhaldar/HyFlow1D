# Agent Instructions

This project implements a hybrid Finite Volume (FV) / Discontinuous Galerkin (DG) solver for the 1D Linear Advection equation.

## Coding Standards

1.  **Language**: C++17.
2.  **Style**:
    *   Use `snake_case` for variables and functions.
    *   Use `PascalCase` for classes.
    *   Use `const` correctness wherever possible.
    *   Use `std::vector` for dynamic arrays.
3.  **Documentation**:
    *   All headers (`.hpp`) should have include guards (`#pragma once` or standard `#ifndef`).
    *   Classes and complex functions should have brief comments explaining their purpose.
4.  **Verification**:
    *   Always compile out-of-source: `mkdir build && cd build && cmake .. && make`.
    *   Run the Python plotting script to verify numerical results visually.

## Directory Structure

*   `src/`: Source files (`.cpp`).
*   `include/`: Header files (`.hpp`).
*   `scripts/`: Python visualization scripts.
*   `tests/`: Unit tests (if any).

## Compilation

```bash
mkdir build
cd build
cmake ..
make
```

## Running

```bash
./hybrid_solver
python3 ../scripts/plot_results.py
```
