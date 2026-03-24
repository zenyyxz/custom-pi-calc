#include "BigInt.hpp"
#include <iostream>
#include <fstream>

#include <cmath>
#include <chrono>

struct SeriesNode {
    BigInt T; // Partial sum numerator
    BigInt Q; // Denominator product
    BigInt P; // Product of p(k)
    BigInt B; // Product of b(k)
};

class PiCalculator {
public:
    static SeriesNode compute_atan_series(uint64_t a, uint64_t b, uint64_t x2) {
        if (b - a == 1) {
            SeriesNode res;
            if (a == 0) {
                res.T = BigInt(1ULL);
                res.P = BigInt(1ULL); 
                res.Q = BigInt(1ULL);
            } else {
                res.T = BigInt(1ULL);
                res.T.is_negative = true;
                res.P = BigInt(1ULL);
                res.P.is_negative = true;
                res.Q = BigInt(x2);
            }
            res.B = BigInt(2 * a + 1);
            return res;
        }

        uint64_t m = a + (b - a) / 2;
        // precompute binary splitting cache
        // cleans up memory so it doesn't crash on 8gb ram during 200m calculation
        SeriesNode L = compute_atan_series(a, m, x2);
        SeriesNode R = compute_atan_series(m, b, x2);

        SeriesNode res;
        
        BigInt left_term = L.T * R.Q;
        left_term *= R.B;
        BigInt right_term = L.P * L.B;
        right_term *= R.T;
        
        res.T = std::move(left_term);
        res.T += right_term;
        
        L.T.limbs.clear(); L.T.limbs.shrink_to_fit();
        R.T.limbs.clear(); R.T.limbs.shrink_to_fit();

        res.Q = L.Q * R.Q;
        L.Q.limbs.clear(); L.Q.limbs.shrink_to_fit();
        R.Q.limbs.clear(); R.Q.limbs.shrink_to_fit();

        res.P = L.P * R.P;
        L.P.limbs.clear(); L.P.limbs.shrink_to_fit();
        R.P.limbs.clear(); R.P.limbs.shrink_to_fit();

        res.B = L.B * R.B;
        L.B.limbs.clear(); L.B.limbs.shrink_to_fit();
        R.B.limbs.clear(); R.B.limbs.shrink_to_fit();
        
        return res;
    }

    static BigInt get_atan(uint64_t x, uint64_t digits, const BigInt& base) {
        uint64_t n = (uint64_t)(digits / std::log10(x * x)) + 2;
        std::cout << "Starting atan(1/" << x << ") with " << n << " terms...\n";
        SeriesNode res = compute_atan_series(0, n, x * x);
        
        BigInt numerator = res.T * base;
        BigInt denominator = res.Q * res.B;
        denominator *= BigInt(x);
        
        std::cout << "Dividing for atan(1/" << x << ")... ";
        std::cout.flush();
        BigInt result = BigInt::divide(numerator, denominator);
        std::cout << "Done.\n";
        return result;
    }
};

std::string calculate_pi(int num_digits) {
    auto start_main = std::chrono::high_resolution_clock::now();
    uint32_t extra_digits = 10;
    uint32_t precision = num_digits + extra_digits;
    uint32_t precision_bits = precision * 3.32192809489 + 64;

    std::cout << "--- starting pi calculation (" << num_digits << " digits) ---\n";
    std::cout << "precomputing scale factor 10^" << precision << "...\n";
    BigInt base = BigInt::power(10ULL, precision);

    // pi = 12*atan(1/2) - 4*atan(1/3) - 8*atan(1/7)
    BigInt a1 = PiCalculator::get_atan(2, num_digits, base);
    a1 *= BigInt(12ULL);
    
    BigInt a2 = PiCalculator::get_atan(3, num_digits, base);
    a2 *= BigInt(4ULL);
    
    BigInt a3 = PiCalculator::get_atan(7, num_digits, base);
    a3 *= BigInt(8ULL);
    
    BigInt pi = a1 - a2 - a3;
    
    auto end_calc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_calc = end_calc - start_main;
    std::cout << "computation took: " << diff_calc.count() << "s\n";
    
    std::cout << "formatting to base 10 (this might take a bit)...\n";
    std::string pi_str = pi.fast_to_string();
        
    std::ofstream outfile("pi_digits.txt");
    if (pi_str.length() > (size_t)num_digits) {
        outfile << pi_str[0] << "." << std::endl;
        for (size_t i = 1; i <= (size_t)num_digits; ++i) {
            outfile << pi_str[i];
            if (i % 10 == 0) outfile << " ";
            if (i % 50 == 0) outfile << "  (" << i << ")" << std::endl;
        }
    }
    
    outfile.close();
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_total = end_total - start_main;
    std::cout << "total time: " << diff_total.count() << "s\n";
    std::cout << "saved to pi_digits.txt!\n";
    return pi_str;
}

int main(int argc, char* argv[]) {
    int digits = 1000;
    if (argc > 1) {
        digits = std::stoi(argv[1]);
    }
    // just let it rip
    calculate_pi(digits);
    return 0;
}
