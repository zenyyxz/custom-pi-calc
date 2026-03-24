#pragma once

#include <vector>
#include <cstdint>

// custom fast fourier transform class over finite fields
class NTT {
public:
    // NTT Primes and their primitive roots
    struct Prime {
        uint64_t p;
        uint64_t g; // primitive root
    };

    // using 4 primes to handle 2^26 array sizes without overflowing crt
    static constexpr Prime P1 = {469762049, 3};
    static constexpr Prime P2 = {1811939329, 13};
    static constexpr Prime P3 = {2013265921, 31};
    static constexpr Prime P4 = {2281701377, 3};

    static uint64_t power(uint64_t base, uint64_t exp, uint64_t mod);
    static uint64_t modInverse(uint64_t n, uint64_t m);
    
    // Core NTT
    static void ntt(std::vector<uint64_t>& a, bool invert, const Prime& p);
    
    // Multi-prime NTT multiplication with CRT
    static void multiply_crt(std::vector<__uint128_t>& res, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);

private:
    static void bit_reverse(std::vector<uint64_t>& a, uint64_t n);
};
