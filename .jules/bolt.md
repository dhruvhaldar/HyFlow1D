## 2024-05-23 - Pre-computation in DG Solver
**Learning:** In high-order numerical methods like Discontinuous Galerkin, small per-element operations accumulate massively. Pre-multiplying quadrature weights into basis functions (`weighted_d_basis`) and pre-calculating inverse mass matrix diagonals (`inv_mass_matrix`) significantly reduces FLOPs in the hot `compute_rhs` loop without increasing memory complexity significantly.
**Action:** Always look for constant terms inside hot loops (integrals, geometric factors) that can be pre-calculated during initialization, especially if they are element-independent or only depend on the mesh geometry.

## 2024-05-24 - Flattened Vector Storage
**Learning:** In C++, `std::vector<std::vector<T>>` incurs significant overhead due to memory fragmentation and double indirection. Flattening 2D data into a single `std::vector<T>` with index arithmetic (`idx = row * width + col`) improved performance by ~1.7% in the DG solver (in a larger benchmark) by enhancing cache locality and reducing allocation overhead.
**Action:** Prefer flat 1D arrays for dense matrix-like structures in performance-critical code, especially when the dimensions are known or fixed at runtime.

## 2024-05-25 - Loop Fusion in DG RHS
**Learning:** In the `compute_rhs` function, fusing the loop that computes solution values at quadrature points with the loop that integrates the flux term eliminated the need for an intermediate storage vector (`u_at_quad_scratch`). This reduced memory traffic (writes/reads to L1 cache) and improved instruction level parallelism, yielding a ~5% speedup.
**Action:** Look for opportunities to compute intermediate values "just-in-time" within the consumer loop to avoid temporary storage, especially when the intermediate data is used immediately and discarded.
