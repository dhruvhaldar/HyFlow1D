## 2024-05-23 - Pre-computation in DG Solver
**Learning:** In high-order numerical methods like Discontinuous Galerkin, small per-element operations accumulate massively. Pre-multiplying quadrature weights into basis functions (`weighted_d_basis`) and pre-calculating inverse mass matrix diagonals (`inv_mass_matrix`) significantly reduces FLOPs in the hot `compute_rhs` loop without increasing memory complexity significantly.
**Action:** Always look for constant terms inside hot loops (integrals, geometric factors) that can be pre-calculated during initialization, especially if they are element-independent or only depend on the mesh geometry.

## 2024-05-24 - Flattened Vector Storage
**Learning:** In C++, `std::vector<std::vector<T>>` incurs significant overhead due to memory fragmentation and double indirection. Flattening 2D data into a single `std::vector<T>` with index arithmetic (`idx = row * width + col`) improved performance by ~16.5% in the DG solver by enhancing cache locality and reducing allocation overhead.
**Action:** Prefer flat 1D arrays for dense matrix-like structures in performance-critical code, especially when the dimensions are known or fixed at runtime.
