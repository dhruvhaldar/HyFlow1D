# Palette's UX Journal

This journal records critical UX and accessibility learnings for the HyFlow1D project.

## 2024-05-22 - CLI Accessibility and Help
**Learning:** Hardcoded simulation parameters in CLI tools hide critical context from the user, making "Help" documentation essential for understanding the tool's behavior without reading source code.
**Action:** Always include "Configuration Details" in CLI help messages when parameters are not user-configurable.

**Learning:** The `NO_COLOR` standard is a simple, high-impact accessibility feature for CLI tools that allows users to disable ANSI colors via environment variables, supporting screen readers and logs.
**Action:** Implement `NO_COLOR` check in CLI entry points.
