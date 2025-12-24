## 2024-05-23 - Flattening Vector<Vector>
**Learning:** Flattening `std::vector<std::vector<double>>` to a single `std::vector<double>` resulted in a ~25% slowdown (1.12s -> 1.40s). This might be due to the overhead of index arithmetic `i * n_modes + k` not being amortized effectively, or compiler auto-vectorization being hindered by the strided access pattern compared to the small, contiguous inner vectors of the original implementation.
**Action:** Do not assume flattening is always faster. Profile or verify with benchmarks.
