## 2024-05-23 - Pre-computation in DG Solver
**Learning:** In high-order numerical methods like Discontinuous Galerkin, small per-element operations accumulate massively. Pre-multiplying quadrature weights into basis functions (`weighted_d_basis`) and pre-calculating inverse mass matrix diagonals (`inv_mass_matrix`) significantly reduces FLOPs in the hot `compute_rhs` loop without increasing memory complexity significantly.
**Action:** Always look for constant terms inside hot loops (integrals, geometric factors) that can be pre-calculated during initialization, especially if they are element-independent or only depend on the mesh geometry.

## 2024-05-24 - Flattened Vector Storage
**Learning:** In C++, `std::vector<std::vector<T>>` incurs significant overhead due to memory fragmentation and double indirection. Flattening 2D data into a single `std::vector<T>` with index arithmetic (`idx = row * width + col`) improved performance by ~1.7% in the DG solver (in a larger benchmark) by enhancing cache locality and reducing allocation overhead.
**Action:** Prefer flat 1D arrays for dense matrix-like structures in performance-critical code, especially when the dimensions are known or fixed at runtime.

## 2024-05-25 - Loop Fusion in DG RHS
**Learning:** In the `compute_rhs` function, fusing the loop that computes solution values at quadrature points with the loop that integrates the flux term eliminated the need for an intermediate storage vector (`u_at_quad_scratch`). This reduced memory traffic (writes/reads to L1 cache) and improved instruction level parallelism, yielding a ~5% speedup.
**Action:** Look for opportunities to compute intermediate values "just-in-time" within the consumer loop to avoid temporary storage, especially when the intermediate data is used immediately and discarded.

## 2025-01-06 - Legendere Derivative Inlining (Rejected)
**Learning:** Inlining the `legendre_derivative` calculation to compute $P_n$ and $P_{n-1}$ in a single pass instead of calling `legendre` twice seemed like a good optimization, but in this specific benchmark, it resulted in a slowdown (25s -> 29.5s).
**Action:** The overhead of `legendre` might be negligible compared to other parts, or the compiler was already optimizing the repeated calls effectively.

## 2025-01-06 - Inlining Boundary Evaluation in Hot Loop
**Learning:** In the `compute_rhs` loop, replacing the generic `evaluate_element(i, 1.0)` call with an inline summation of coefficients (since $P_k(1) = 1$) reduced overhead significantly. This is because `evaluate_element` involves branching for `xi = 1.0` and potential function call overhead, which adds up inside the hot loop over elements and time steps.
**Action:** For boundary values in hot loops, specialize the evaluation manually if the basis function properties (like $P_k(1)=1$) allow for a simple sum, avoiding the cost of a general-purpose evaluation function.

## 2025-01-14 - Precomputed Stiffness Matrix
**Learning:** Precomputing the stiffness matrix ($K_{km} = \int P_m P'_k d\xi$) allows replacing the complex quadrature loop in `compute_rhs` with a simple matrix-vector multiplication. This reduced the complexity from $O(N_{modes}^2)$ ops per element (with large constant factor from quadrature) to a tight $O(N_{modes}^2)$ mat-vec, resulting in a ~43% speedup.
**Action:** When solving PDEs with fixed basis functions, always try to precompute operator matrices (Stiffness, Mass, Advection) instead of performing quadrature integration on the fly, as the geometric factors and basis integrals are often constant.
