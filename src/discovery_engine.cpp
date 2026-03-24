#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

int main() {
    std::cout << "--- High-Precision Math Discovery Engine ---\n";
    const double TARGET = M_PI / 4.0;
    const double EPSILON = 1e-14;

    // Searching for c1*atan(1/x1) + c2*atan(1/x2) + c3*atan(1/x3) = pi/4
    // We want larger x values for faster convergence in the BigInt phase.
    for (int x1 = 2; x1 < 20; ++x1) {
        double v1 = std::atan(1.0/x1);
        for (int x2 = x1 + 1; x2 < 200; ++x2) {
            double v2 = std::atan(1.0/x2);
            for (int x3 = x2 + 1; x3 < 2000; ++x3) {
                double v3 = std::atan(1.0/x3);
                for (int c1 = 1; c1 <= 15; ++c1) {
                    for (int c2 = -15; c2 <= 15; ++c2) {
                        if (c2 == 0) continue;
                        double partial = c1 * v1 + c2 * v2;
                        // Optimization: if partial is already far from TARGET, 
                        // x3 might not be able to pull it back.
                        for (int c3 = -15; c3 <= 15; ++c3) {
                            if (c3 == 0) continue;
                            
                            double res = partial + c3 * v3;
                            if (std::abs(res - TARGET) < EPSILON) {
                                std::cout << "FOUND VALID FORMULA:\n";
                                std::cout << "pi/4 = " << c1 << "*atan(1/" << x1 << ") + (" 
                                          << c2 << ")*atan(1/" << x2 << ") + (" 
                                          << c3 << ")*atan(1/" << x3 << ")\n";
                                // Verify with higher precision if possible, but for discovery this is good.
                                return 0;
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "No formula found in this range.\n";
    return 0;
}
