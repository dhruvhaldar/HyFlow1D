## 2024-12-25 - [Input Validation for Mathematical Limits]
**Vulnerability:** Unchecked polynomial order in DG solver caused a runtime crash (`std::runtime_error`) when it exceeded the hardcoded limit of the quadrature library.
**Learning:** Mathematical utility functions often have implementation limits (e.g., number of precomputed quadrature points). These limits must be propagated as input validation checks in the calling layers to ensure "Fail Securely" and avoid unexpected runtime termination.
**Prevention:** Always validate configuration parameters against the capabilities of dependent libraries/modules at initialization time.

## 2024-05-23 - [C++ Compiler Hardening]
**Vulnerability:** Lack of explicit compiler hardening flags (`-fstack-protector-strong`, `_FORTIFY_SOURCE=2`, RELRO, PIE).
**Learning:** While some modern toolchains enable basic protections (like basic stack protector and PIE) by default, they don't always enable the stronger variants or full RELRO. Explicitly setting them in `CMakeLists.txt` ensures consistent security posture across different build environments.
**Prevention:** Always include a standard set of hardening flags in the project's build configuration to enforce "Defense in Depth".

## 2024-05-23 - [C++ Build Configuration]
**Vulnerability:** Misconfiguration of `_FORTIFY_SOURCE=2` in Debug builds caused build failures due to `-Werror`.
**Learning:** `_FORTIFY_SOURCE` requires optimization (typically -O1 or higher) to work. Defining it in Debug builds (which usually use -O0) can cause compiler warnings (e.g., `#warning _FORTIFY_SOURCE requires compiling with optimization`). When `-Werror` is active, this breaks the build.
**Prevention:** Use generator expressions (e.g., `$<$<CONFIG:Release>:_FORTIFY_SOURCE=2>`) to conditionally apply flags based on the build type.

## 2025-01-16 - [Filesystem Race Conditions (TOCTOU)]
**Vulnerability:** A Time-of-Check Time-of-Use (TOCTOU) vulnerability existed where `fs::exists` was checked before creating a directory and setting permissions. An attacker could introduce a symlink between the check and the usage, causing `fs::permissions` to modify an arbitrary target file/directory.
**Learning:** `fs::exists` is not atomic with subsequent operations. Checks on filesystem state are inherently racy.
**Prevention:** Rely on the return value of atomic operations like `fs::create_directory` (which returns `true` only if it created the directory). Additionally, use `fs::perm_options::nofollow` when setting permissions to prevent traversing unexpected symlinks.

## 2025-01-20 - [Uninitialized State & Out-of-Bounds Access]
**Vulnerability:** The `DiscontinuousGalerkinSolver` methods could be called before `initialize()`, leading to usage of uninitialized indices and empty vectors (Segfault). Also, internal helper `evaluate_element` lacked bounds checking, allowing potential out-of-bounds reads.
**Learning:** C++ classes with separate initialization methods (two-phase init) are prone to "Use Before Init" bugs. Internal helpers often assume valid state, but can become security liabilities if exposed or misused.
**Prevention:** Explicitly track initialization state (`is_initialized`) and enforce checks in all public entry points. Apply "Defense in Depth" by adding bounds checking even in internal helpers.
## 2025-01-20 - [Input Validation for File Paths]
**Vulnerability:** The application accepted arbitrary file paths for the output directory (`-o`), allowing Path Traversal attacks (e.g., `-o ../etc`) which could lead to arbitrary file writes outside the intended directory.
**Learning:** CLI tools often trust user input for file paths implicitly. Relying on the OS to handle permissions is insufficient defense-in-depth, especially if the tool might run with elevated privileges or in shared environments.
**Prevention:** Explicitly validate all user-provided file paths. Reject paths containing directory traversal components (`..`) unless specifically required and authorized. Use `std::filesystem` to parse and inspect path components reliably.

## 2025-01-22 - [Unsecured Helper Scripts]
**Vulnerability:** The Python visualization script accepted arbitrary file paths with traversal ('..'), allowing arbitrary file read/write operations when used with crafted arguments.
**Learning:** Security measures applied to the main application (C++ CLI) were not mirrored in auxiliary scripts, creating a weak link in the toolchain. Scripts are often overlooked in security reviews but can be vectors for local privilege escalation or data exfiltration if run by privileged users.
**Prevention:** Apply the same input validation standards (like 'is_safe_path') to all components of the ecosystem, including scripts and tools.

## 2025-01-22 - [Symlink Attack Prevention in Scripts]
**Vulnerability:** The `plot_results.py` script relied on standard file I/O (`plt.savefig`) which follows symbolic links. This allowed an attacker to create a symlink (e.g., `output.png -> /etc/passwd`) and trick the script into overwriting sensitive files when run by a privileged user.
**Learning:** Standard file I/O libraries (like `open`, `matplotlib`, `pandas`) typically follow symlinks by default. Validating the path string alone (e.g., for `..`) is insufficient against symlink attacks.
**Prevention:** Explicitly check `Path(path).is_symlink()` before writing to output files in CLI tools. Refuse to write if a symlink is detected.

## 2025-01-22 - [Symlink Overwrite in C++]
**Vulnerability:** The C++ `std::ofstream` follows symbolic links by default. Writing to a user-controlled path (even inside a safe directory) allowed overwriting arbitrary files via pre-created symlinks.
**Learning:** Checking directory safety ('..') is not enough. File writers must explicitly check for symlinks if the output location could be tampered with. Standard C++ streams lack `O_NOFOLLOW`.
**Prevention:** Use `std::filesystem::is_symlink()` before opening files for writing in security-sensitive contexts, while acknowledging the residual TOCTOU risk.

## 2025-01-22 - [Insecure Default Permissions Race Condition]
**Vulnerability:** Files and directories were created with default process permissions (often world-readable) before being restricted to owner-only using `fs::permissions`. This created a race condition where sensitive data was briefly exposed.
**Learning:** Post-creation permission hardening (e.g., `chmod` after `open`) is insufficient for confidentiality because it is not atomic. The default creation mask (`umask`) determines the initial permissions.
**Prevention:** Set the process `umask` (e.g., `0077`) at application startup to ensuring all subsequently created files and directories are secure by default (atomic security).

## 2025-01-26 - [Undefined Behavior in Constructor Initialization]
**Vulnerability:** The `DiscontinuousGalerkinSolver` constructor computed a member variable (`n_modes = p_order + 1`) in the initialization list before validating `p_order` in the body. Passing `INT_MAX` caused Signed Integer Overflow (Undefined Behavior) before validation could occur.
**Learning:** In C++, member initialization happens before the constructor body is executed. Validation logic placed inside the body is "too late" if the initialization expressions themselves are unsafe.
**Prevention:** Use static helper functions to validate input arguments *before* they are used in the initialization list (e.g., `: member(validate(input))`).

## 2026-01-26 - [Secure File Opening with High-Level Libraries]
**Vulnerability:** A TOCTOU race condition existed in a Python script where `Path.is_symlink()` was checked before `plt.savefig()`. An attacker could switch the file to a symlink between the check and the write.
**Learning:** Checking for symlinks in user-space (`is_symlink`) before opening is inherently racy. High-level libraries like `matplotlib` do not expose `O_NOFOLLOW` options directly.
**Prevention:** Use low-level `os.open` with `O_NOFOLLOW | O_CREAT | O_TRUNC` and secure permissions (`0o600`) to get a file descriptor, then wrap it with `os.fdopen` and pass the file object to the high-level library. This ensures atomic security while retaining the convenience of the library.
