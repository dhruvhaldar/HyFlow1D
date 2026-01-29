## 2024-05-23 - DG Stiffness Matrix Sparsity
**Learning:** The stiffness matrix $K_{km} = \int P_m P'_k$ for Legendre polynomials is strictly lower triangular ($K_{km} = 0$ for $m \ge k$). This is because $P'_k$ is a polynomial of degree $k-1$, and due to orthogonality, its inner product with any $P_m$ where $m > k-1$ is zero.
**Action:** Always check mathematical properties of matrices in numerical methods. Loop bounds can be optimized from `0..N` to `0..k` (or similar) to skip zero blocks, providing significant speedups (e.g., ~50% reduction in MACs for small N).

## 2024-05-24 - Stack Allocation and Restrict Pointers
**Learning:** Copying small invariant matrices to the stack and using `__restrict__` qualifiers significantly helps the compiler optimize inner loops, likely by reducing aliasing assumptions and improving register allocation. In this case, it yielded a ~17% speedup in the DG solver's RHS computation.
**Action:** For critical loops involving small, repeated matrix ops, consider explicit stack copies and strict aliasing hints.

## 2024-05-24 - Fully Unrolling Small Matrix Ops
**Learning:** For small fixed sizes (e.g., N=1..4), fully unrolling loops using `if constexpr` and scalar variables (avoiding temporary arrays) eliminates loop overhead and allows better instruction scheduling. This yielded a ~27% speedup (904k -> 1.15M ops/s) by removing complex pointer arithmetic and dependency chains in the hot path.
**Action:** When N is small and known at compile time, prefer full unrolling with `if constexpr` over loops, even if the loop count is small.

## 2024-05-24 - Polynomial Evaluation Complexity
**Learning:** Evaluating a polynomial expansion $\sum u_k P_k(x)$ by calling `legendre(k, x)` for each term creates an $O(N^2)$ bottleneck because each call re-computes $P_0 \dots P_k$. Using the inline recurrence relation reduces this to $O(N)$ and yields a ~27% speedup in the evaluation kernel while preserving code clarity.
**Action:** Always prefer inline recurrence (e.g., Clenshaw-like) for evaluating orthogonal polynomial series over independent function calls.

## 2024-10-24 - FV Solver Stencil Optimization
**Learning:** For 1D stencil operations (like Finite Volume `u[i] - u[i-1]`), manually caching the previous element in a register variable (`u_prev`) combined with `__restrict__` pointers yielded a 3x speedup. The compiler failed to optimize the redundant memory load of `u[i-1]` (which was `u[i]` in the previous iteration) automatically, likely due to potential aliasing concerns or conservative analysis.
**Action:** When implementing stencil loops, explicitly use local variables to carry over dependencies and use `__restrict__` pointers to hint the compiler, even for simple 1st-order schemes.

## 2024-10-25 - Auto-Vectorization vs Scalar Replacement
**Learning:** Contrary to the 2024-10-24 finding, applying manual scalar replacement (`u_prev`) in the FV solver loop caused a 60% performance regression (210k -> 83k ops/s) with GCC 13.3. The compiler successfully auto-vectorizes the simple array access `u[i] - u[i-1]`, while the manual scalar dependency inhibits vectorization.
**Action:** Prioritize array indexing that enables auto-vectorization over manual scalar replacement for simple stencils. Verification via benchmarking is essential before applying such "optimizations".

## 2024-10-26 - DG Stiffness Matrix Sparsity Unrolling
**Learning:** The stiffness matrix $K_{km}$ has specific zero entries due to parity (k+m even -> 0). Exploiting this in the fully unrolled kernel for N=3 and N=4 saved ~3% overhead by removing 2 multiplies per element.
**Action:** When unrolling matrix operations, always check for mathematical zeros (sparsity) that might not be obvious in the generic loop form.
