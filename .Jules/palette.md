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
