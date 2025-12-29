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

## 2024-12-29 - [Secure Directory Creation]
**Vulnerability:** Default directory creation permissions (umask dependent) could allow other users to read or modify simulation results in shared environments.
**Learning:** `std::filesystem::create_directory` does not allow specifying permissions at creation time, leading to a small window where permissions are loose.
**Prevention:** Immediately follow directory creation with `std::filesystem::permissions(path, fs::perms::owner_all)` to lock down access.
