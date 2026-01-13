## 2024-05-23 - DG Stiffness Matrix Sparsity
**Learning:** The stiffness matrix $K_{km} = \int P_m P'_k$ for Legendre polynomials is strictly lower triangular ($K_{km} = 0$ for $m \ge k$). This is because $P'_k$ is a polynomial of degree $k-1$, and due to orthogonality, its inner product with any $P_m$ where $m > k-1$ is zero.
**Action:** Always check mathematical properties of matrices in numerical methods. Loop bounds can be optimized from `0..N` to `0..k` (or similar) to skip zero blocks, providing significant speedups (e.g., ~50% reduction in MACs for small N).

## 2024-05-24 - Aliasing and Recurrence Optimization
**Learning:** Even with `const double*`, compilers may assume aliasing with output pointers (`double* rhs`). Hoisting small invariant data (like stiffness matrix or mass matrix diagonals) into local stack arrays forces values into registers/L1 cache and eliminates aliasing checks, improving vectorization opportunities. Also, standard Legendre polynomial evaluation is $O(N^2)$ if implemented naively; using the recurrence relation inline makes it $O(N)$.
**Action:** Use `__restrict__` for tight numerical loops and copy small constant lookup tables to local stack variables. Always check if polynomial evaluations can be done via recurrence in a single pass.
