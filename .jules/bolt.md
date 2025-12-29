## 2024-05-23 - [Transposed Basis for DG Solver]
**Learning:** In numerical solvers (DG/FEM), accessing basis functions in `[quad_point][mode]` layout is intuitive but causes cache misses during volume integration, which typically iterates `mode` then `quad_point`. Strided access into `std::vector<std::vector<>>` exacerbates this with pointer chasing.
**Action:** Always verify memory layout matches the innermost loop. Transposing to `[mode][quad_point]` ensures contiguous memory access. Pre-multiplying weights into the basis further reduces FLOPS.
