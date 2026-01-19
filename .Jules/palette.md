# Palette's UX Journal

This journal records critical UX and accessibility learnings for the HyFlow1D project.

## 2024-05-22 - CLI Accessibility and Help
**Learning:** Hardcoded simulation parameters in CLI tools hide critical context from the user, making "Help" documentation essential for understanding the tool's behavior without reading source code.
**Action:** Always include "Configuration Details" in CLI help messages when parameters are not user-configurable.

**Learning:** The `NO_COLOR` standard is a simple, high-impact accessibility feature for CLI tools that allows users to disable ANSI colors via environment variables, supporting screen readers and logs.
**Action:** Implement `NO_COLOR` check in CLI entry points.

## 2024-05-23 - Relative Paths in CLI Hints
**Learning:** CLI tools often suggest follow-up commands (e.g., visualization scripts). However, standard build workflows (like `mkdir build && cd build`) change the working directory, breaking relative path assumptions if the suggestion is hardcoded.
**Action:** Use runtime path detection (e.g., `std::filesystem`) to verify existence of helper scripts and dynamically adjust the suggested command path.

## 2024-05-24 - Typo Tolerance in CLI
**Learning:** Users frequently mistype short flags (e.g., `--outupt` vs `--output`). Implementing a simple "Did you mean?" suggestion system using Levenshtein distance significantly reduces friction and makes the tool feel more helpful rather than just rigid.
**Action:** Implement fuzzy matching for unknown arguments in CLI entry points.

## 2024-05-25 - Cross-Language CLI Consistency
**Learning:** When a project uses multiple languages for its toolchain (e.g., C++ solver + Python visualizer), users expect a consistent visual language (colors, emoji use) and behavior (auto-detection of build artifacts) across all tools.
**Action:** Ensure helper scripts implement the same `NO_COLOR` standards and success/error styling as the main application.

## 2024-05-26 - Metadata in Simulation Output
**Learning:** Simulation users think in terms of physics time, not iteration steps. Abstract filenames (e.g., `step_100.csv`) forcing users to manually calculate time creates unnecessary cognitive load.
**Action:** Embed critical metadata (like simulation time) directly into output file headers (e.g., as comments in CSV) to allow visualization tools to display meaningful labels automatically.

## 2024-10-27 - Actionable Error States
**Learning:** Error messages that only state the problem ("No files found") leave the user stranded. Providing the immediate next step ("Run ./hyflow1d") transforms a dead-end into a guided workflow.
**Action:** Always pair error states with a "Tip" or "Suggested Action" that solves the problem.

## 2024-10-27 - Process Transparency
**Learning:** When a tool performs a batch operation (like plotting multiple files), summarizing the selection *before* execution gives the user confidence that the correct data is being processed, especially when sampling is involved.
**Action:** Print a summary of selected inputs (filenames, timestamps) before starting long-running or batch operations.

## 2024-10-27 - Graceful Interrupts
**Learning:** Users often interrupt long-running simulations (Ctrl+C) to check partial results, not just to abort. Treating `SIGINT` as a "Pause" event with helpful context (step number, saved files) and next steps (visualization command) transforms a crash into a useful workflow.
**Action:** Catch `SIGINT`, print a "Paused" status to `stderr`, and provide the visualization command for partial results.
