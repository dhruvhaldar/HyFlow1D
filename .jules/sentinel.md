## 2024-12-25 - [Input Validation for Mathematical Limits]
**Vulnerability:** Unchecked polynomial order in DG solver caused a runtime crash (`std::runtime_error`) when it exceeded the hardcoded limit of the quadrature library.
**Learning:** Mathematical utility functions often have implementation limits (e.g., number of precomputed quadrature points). These limits must be propagated as input validation checks in the calling layers to ensure "Fail Securely" and avoid unexpected runtime termination.
**Prevention:** Always validate configuration parameters against the capabilities of dependent libraries/modules at initialization time.
