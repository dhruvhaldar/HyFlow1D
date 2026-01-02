## 2024-10-24 - CLI Progress Indicators
**Learning:** Adding a simple ASCII progress bar to a C++ CLI tool significantly improves the "delight" factor compared to scrolling text logs.
**Action:** For CLI tools with long-running loops, default to a progress bar and hide detailed logs behind a `-v` flag.
## 2024-05-23 - [CLI Colorization]
**Learning:** Adding ANSI color codes to C++ CLI tools dramatically improves user feedback (success vs. error) with minimal code overhead.
**Action:** Use a simple header-only `Color` namespace with `constexpr std::string_view` for zero-cost abstractions in future C++ CLI tasks.
## 2024-05-23 - [CLI Progress Bar Upgrade]
**Learning:** Adding visual block characters and an ETA timer transforms a basic CLI tool into a professional-feeling application.
**Action:** Always check if a long-running CLI process has an ETA indicator. If not, add one.
## 2024-10-24 - [CLI Argument Validation]
**Learning:** Silent failure on unknown flags is a common CLI frustration. Catching "unknown arguments" early prevents user confusion about why their flag (e.g., typo) isn't working.
**Action:** Always implement a strict `else` block in argument parsing loops to reject unknown flags/arguments.
## 2024-10-24 - [CLI Configuration Summary]
**Learning:** Displaying a configuration summary before execution (even if parameters are hardcoded) builds user trust and clarity.
**Action:** Always print a formatted summary of simulation parameters at startup.
