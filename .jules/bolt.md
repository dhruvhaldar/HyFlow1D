## 2024-05-23 - DG Stiffness Matrix Sparsity
**Learning:** The stiffness matrix $K_{km} = \int P_m P'_k$ for Legendre polynomials is strictly lower triangular ($K_{km} = 0$ for $m \ge k$). This is because $P'_k$ is a polynomial of degree $k-1$, and due to orthogonality, its inner product with any $P_m$ where $m > k-1$ is zero.
**Action:** Always check mathematical properties of matrices in numerical methods. Loop bounds can be optimized from `0..N` to `0..k` (or similar) to skip zero blocks, providing significant speedups (e.g., ~50% reduction in MACs for small N).
