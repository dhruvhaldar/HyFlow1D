## 2024-10-24 - CLI Progress Indicators
**Learning:** Adding a simple ASCII progress bar to a C++ CLI tool significantly improves the "delight" factor compared to scrolling text logs.
**Action:** For CLI tools with long-running loops, default to a progress bar and hide detailed logs behind a `-v` flag.
## 2024-05-23 - [CLI Colorization]
**Learning:** Adding ANSI color codes to C++ CLI tools dramatically improves user feedback (success vs. error) with minimal code overhead.
**Action:** Use a simple header-only `Color` namespace with `constexpr std::string_view` for zero-cost abstractions in future C++ CLI tasks.
