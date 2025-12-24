#pragma once

#include <vector>
#include <cmath>

namespace numerics {

    // Evaluate Legendre polynomial P_n(x) at x
    inline double legendre(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        double p_prev = 1.0; // P_0
        double p_curr = x;   // P_1
        double p_next = 0.0;

        for (int k = 1; k < n; ++k) {
            p_next = ((2.0 * k + 1.0) * x * p_curr - k * p_prev) / (k + 1.0);
            p_prev = p_curr;
            p_curr = p_next;
        }
        return p_curr;
    }

    // Evaluate derivative of Legendre polynomial P'_n(x) at x
    inline double legendre_derivative(int n, double x) {
        if (n == 0) return 0.0;
        if (n == 1) return 1.0;

        // Relation: (x^2 - 1) P'_n(x) = n (x P_n(x) - P_{n-1}(x))
        double Pn = legendre(n, x);
        double Pn_minus_1 = legendre(n - 1, x);

        if (std::abs(x) > 1.0 - 1e-12) {
            // Handle edge cases x = +/- 1
            // P'_n(1) = n(n+1)/2
            // P'_n(-1) = (-1)^(n-1) * n(n+1)/2
            double val = n * (n + 1.0) / 2.0;
            return (x > 0) ? val : (val * std::pow(-1, n - 1));
        }

        return (n * (x * Pn - Pn_minus_1)) / (x * x - 1.0);
    }

    // Get Gauss-Legendre quadrature nodes and weights for N points
    // Returns a pair of vectors: {nodes, weights}
    // Nodes are in [-1, 1]
    std::pair<std::vector<double>, std::vector<double>> gauss_legendre(int n);

}
