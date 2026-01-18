## 2024-05-23 - DG Stiffness Matrix Sparsity
**Learning:** The stiffness matrix $K_{km} = \int P_m P'_k$ for Legendre polynomials is strictly lower triangular ($K_{km} = 0$ for $m \ge k$). This is because $P'_k$ is a polynomial of degree $k-1$, and due to orthogonality, its inner product with any $P_m$ where $m > k-1$ is zero.
**Action:** Always check mathematical properties of matrices in numerical methods. Loop bounds can be optimized from `0..N` to `0..k` (or similar) to skip zero blocks, providing significant speedups (e.g., ~50% reduction in MACs for small N).

## 2024-05-24 - Stack Allocation and Restrict Pointers
**Learning:** Copying small invariant matrices to the stack and using `__restrict__` qualifiers significantly helps the compiler optimize inner loops, likely by reducing aliasing assumptions and improving register allocation. In this case, it yielded a ~17% speedup in the DG solver's RHS computation.
**Action:** For critical loops involving small, repeated matrix ops, consider explicit stack copies and strict aliasing hints.

## 2024-05-24 - Fully Unrolling Small Matrix Ops
**Learning:** For small fixed sizes (e.g., N=1..4), fully unrolling loops using `if constexpr` and scalar variables (avoiding temporary arrays) eliminates loop overhead and allows better instruction scheduling. This yielded a ~27% speedup (904k -> 1.15M ops/s) by removing complex pointer arithmetic and dependency chains in the hot path.
**Action:** When N is small and known at compile time, prefer full unrolling with `if constexpr` over loops, even if the loop count is small.

## 2024-05-25 - FV Solver Branch Elimination
**Learning:** In the Finite Volume solver, the `if (i==0)` check inside the loop for boundary conditions inhibited vectorization and added branching overhead. Peeling the first iteration and precalculating the coefficient (`-a/dx`) improved throughput by ~25% (7.9k -> 9.9k ops/s).
**Action:** When handling boundary conditions in tight loops, peel the boundary iterations to allow the main loop to run without conditionals, enabling better compiler optimization.
