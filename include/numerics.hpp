#pragma once

#include <vector>
#include <cmath>

namespace numerics {

    // Evaluate Legendre polynomial P_n(x) at x
    double legendre(int n, double x);

    // Evaluate derivative of Legendre polynomial P'_n(x) at x
    double legendre_derivative(int n, double x);

    // Get Gauss-Legendre quadrature nodes and weights for N points
    // Returns a pair of vectors: {nodes, weights}
    // Nodes are in [-1, 1]
    std::pair<std::vector<double>, std::vector<double>> gauss_legendre(int n);

}
