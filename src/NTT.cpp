#include "NTT.hpp"
#include <algorithm>

uint64_t NTT::power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__uint128_t)res * base % mod;
        base = (__uint128_t)base * base % mod;
        exp /= 2;
    }
    return res;
}

uint64_t NTT::modInverse(uint64_t n, uint64_t mod) {
    return power(n, mod - 2, mod);
}

void NTT::bit_reverse(std::vector<uint64_t>& a, uint64_t n) {
    for (uint64_t i = 1, j = 0; i < n; i++) {
        uint64_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

void NTT::ntt(std::vector<uint64_t>& a, bool invert, const Prime& p) {
    uint64_t n = a.size();
    bit_reverse(a, n);

    for (uint64_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = power(p.g, (p.p - 1) / len, p.p);
        if (invert) wlen = modInverse(wlen, p.p);
        for (uint64_t i = 0; i < n; i += len) {
            uint64_t w = 1;
            for (uint64_t j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = (__uint128_t)a[i + j + len / 2] * w % p.p;
                a[i + j] = (u + v) % p.p;
                a[i + j + len / 2] = (u + p.p - v) % p.p;
                w = (__uint128_t)w * wlen % p.p;
            }
        }
    }

    if (invert) {
        uint64_t n_inv = modInverse(n, p.p);
        for (uint64_t& x : a) x = (__uint128_t)x * n_inv % p.p;
    }
}

void NTT::multiply_crt(std::vector<__uint128_t>& res, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    uint64_t n = 1;
    while (n < a.size() + b.size()) n <<= 1;

    auto work = [&](const Prime& p) {
        std::vector<uint64_t> fa(n, 0), fb(n, 0);
        for (size_t i = 0; i < a.size(); ++i) fa[i] = a[i] % p.p;
        for (size_t i = 0; i < b.size(); ++i) fb[i] = b[i] % p.p;
        ntt(fa, false, p);
        ntt(fb, false, p);
        for (uint64_t i = 0; i < n; i++) fa[i] = (__uint128_t)fa[i] * fb[i] % p.p;
        ntt(fa, true, p);
        return fa;
    };

    std::vector<uint64_t> r1 = work(P1);
    std::vector<uint64_t> r2 = work(P2);
    std::vector<uint64_t> r3 = work(P3);
    std::vector<uint64_t> r4 = work(P4);

    uint64_t invP1_P2 = modInverse(P1.p % P2.p, P2.p);
    
    uint64_t P1P2_modP3 = ((__uint128_t)P1.p * P2.p) % P3.p;
    uint64_t invP1P2_P3 = modInverse(P1P2_modP3, P3.p);
    
    __uint128_t P1P2 = (__uint128_t)P1.p * P2.p;
    uint64_t P1P2_modP4 = P1P2 % P4.p;
    uint64_t P1P2P3_modP4 = ((__uint128_t)P1P2_modP4 * P3.p) % P4.p;
    uint64_t invP1P2P3_P4 = modInverse(P1P2P3_modP4, P4.p);

    res.resize(n);
    for (uint64_t i = 0; i < n; i++) {
        uint64_t r1_modP2 = r1[i] % P2.p;
        uint64_t k1 = (__uint128_t)(r2[i] + P2.p - r1_modP2) * invP1_P2 % P2.p;
        __uint128_t x2 = r1[i] + (__uint128_t)k1 * P1.p; // mod P1*P2
        
        uint64_t x2_modP3 = x2 % P3.p;
        uint64_t k2 = (__uint128_t)(r3[i] + P3.p - x2_modP3) * invP1P2_P3 % P3.p;
        __uint128_t x3 = x2 + (__uint128_t)k2 * P1P2; // mod P1*P2*P3
        
        uint64_t x3_modP4 = x3 % P4.p;
        uint64_t k3 = (__uint128_t)(r4[i] + P4.p - x3_modP4) * invP1P2P3_P4 % P4.p;
        __uint128_t x4 = x3 + (__uint128_t)k3 * P1P2 * P3.p; // target element
        
        res[i] = x4;
    }
}
