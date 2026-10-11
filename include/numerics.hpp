#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "solver.hpp"

namespace numerics {

    // Physical conservation law flux f(u)
    inline double flux(double u, EquationType eq, double a = 1.0) {
        if (eq == EquationType::LinearAdvection) {
            return a * u;
        } else {
            return 0.5 * u * u;
        }
    }

    // Maximum local characteristic wave speed |f'(u)|
    inline double wave_speed(double u, EquationType eq, double a = 1.0) {
        if (eq == EquationType::LinearAdvection) {
            return std::abs(a);
        } else {
            return std::abs(u);
        }
    }

    // Local Lax-Friedrichs (Rusanov) numerical interface flux:
    // F(u_L, u_R) = 0.5 * (f(u_L) + f(u_R)) - 0.5 * lambda_max * (u_R - u_L)
    inline double rusanov_flux(double u_left, double u_right, EquationType eq, double a = 1.0) {
        double f_l = flux(u_left, eq, a);
        double f_r = flux(u_right, eq, a);
        double lambda = std::max(wave_speed(u_left, eq, a), wave_speed(u_right, eq, a));
        return 0.5 * (f_l + f_r) - 0.5 * lambda * (u_right - u_left);
    }

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
