#include "numerics.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace numerics {

    double legendre(int n, double x) {
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

    double legendre_derivative(int n, double x) {
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

    std::pair<std::vector<double>, std::vector<double>> gauss_legendre(int n) {
        if (n < 1) throw std::runtime_error("Number of quadrature points must be >= 1");

        std::vector<double> nodes(n);
        std::vector<double> weights(n);

        // Standard values for small n to avoid complex root finding logic in PoC
        // For arbitrary n, we would use Newton-Raphson on P_n(x).
        // For this PoC, we will implement N=1 to N=5 which is sufficient for up to P=4 DG.
        
        switch (n) {
            case 1:
                nodes = {0.0};
                weights = {2.0};
                break;
            case 2:
                nodes = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
                weights = {1.0, 1.0};
                break;
            case 3:
                nodes = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
                weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
                break;
            case 4:
                {
                    double x1 = std::sqrt(3.0/7.0 - 2.0/7.0 * std::sqrt(6.0/5.0));
                    double x2 = std::sqrt(3.0/7.0 + 2.0/7.0 * std::sqrt(6.0/5.0));
                    double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
                    double w2 = (18.0 - std::sqrt(30.0)) / 36.0;
                    nodes = {-x2, -x1, x1, x2};
                    weights = {w2, w1, w1, w2};
                }
                break;
            case 5:
                {
                    double x1 = 1.0/3.0 * std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0));
                    double x2 = 1.0/3.0 * std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0));
                    double w1 = (322.0 + 13.0*std::sqrt(70.0)) / 900.0;
                    double w2 = (322.0 - 13.0*std::sqrt(70.0)) / 900.0;
                    nodes = {-x2, -x1, 0.0, x1, x2};
                    weights = {w2, w1, 128.0/225.0, w1, w2};
                }
                break;
            default:
                 // Fallback: Just return 1 point (order 0 accuracy) or throw.
                 // Ideally we implement the root finder.
                 // For the purpose of this task (P=3 or P=4), N=5 is enough.
                 throw std::runtime_error("Gauss-Legendre > 5 points not implemented in this PoC.");
        }

        return {nodes, weights};
    }

}
