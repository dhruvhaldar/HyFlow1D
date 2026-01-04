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
