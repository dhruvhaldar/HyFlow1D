#include <iostream>
#include <cassert>
#include <cmath>
#include "numerics.hpp"

int main() {
    std::cout << "Running Numerics Tests..." << std::endl;

    // Test Legendre P_0(x) = 1
    assert(std::abs(numerics::legendre(0, 0.5) - 1.0) < 1e-9);
    
    // Test Legendre P_1(x) = x
    assert(std::abs(numerics::legendre(1, 0.5) - 0.5) < 1e-9);

    // Test Legendre P_2(x) = 0.5 * (3x^2 - 1)
    double p2_05 = 0.5 * (3 * 0.5 * 0.5 - 1);
    assert(std::abs(numerics::legendre(2, 0.5) - p2_05) < 1e-9);

    // Test Quadrature Integration of f(x) = x^2 over [-1, 1] -> Result should be 2/3
    // Use N=2 (exact for polynomials up to degree 2*2-1 = 3)
    auto [nodes, weights] = numerics::gauss_legendre(2);
    double integral = 0.0;
    for(size_t i=0; i<nodes.size(); ++i) {
        double x = nodes[i];
        integral += weights[i] * (x * x);
    }
    assert(std::abs(integral - 2.0/3.0) < 1e-9);

    std::cout << "All Numerics Tests Passed!" << std::endl;
    return 0;
}
