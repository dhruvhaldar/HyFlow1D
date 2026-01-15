## 2024-05-23 - DG Stiffness Matrix Sparsity
**Learning:** The stiffness matrix $K_{km} = \int P_m P'_k$ for Legendre polynomials is strictly lower triangular ($K_{km} = 0$ for $m \ge k$). This is because $P'_k$ is a polynomial of degree $k-1$, and due to orthogonality, its inner product with any $P_m$ where $m > k-1$ is zero.
**Action:** Always check mathematical properties of matrices in numerical methods. Loop bounds can be optimized from `0..N` to `0..k` (or similar) to skip zero blocks, providing significant speedups (e.g., ~50% reduction in MACs for small N).

## 2024-05-24 - Stack Allocation and Restrict Pointers
**Learning:** Copying small invariant matrices to the stack and using `__restrict__` qualifiers significantly helps the compiler optimize inner loops, likely by reducing aliasing assumptions and improving register allocation. In this case, it yielded a ~17% speedup in the DG solver's RHS computation.
**Action:** For critical loops involving small, repeated matrix ops, consider explicit stack copies and strict aliasing hints.
