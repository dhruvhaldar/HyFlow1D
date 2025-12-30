## 2024-05-21 - Vector of Vectors Anti-Pattern
**Learning:** Using `std::vector<std::vector<double>>` for the state vector `u` (1000 elements x 4 modes) creates excessive heap allocations and pointer indirections. Each element `u[i]` is a separate heap allocation, leading to poor cache locality and memory fragmentation.
**Action:** Use a flat `std::vector<double>` with stride arithmetic (`u[i * n_modes + k]`) for dense field data. This ensures contiguous memory, reduces allocator overhead, and improves prefetching.
