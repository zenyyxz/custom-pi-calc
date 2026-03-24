#include "BigInt.hpp"
#include "NTT.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#else
#error "This implementation requires x86_64 for intrinsic assembly instructions."
#endif

// Forward declarations
void ntt_multiply_adapter(std::vector<uint64_t>& res, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);
static std::string convert_rec(const BigInt& x, const std::vector<BigInt>& p10, int k, bool pad);

BigInt::BigInt() : is_negative(false) {
    limbs.push_back(0);
}

BigInt::BigInt(const std::vector<uint64_t>& l) : is_negative(false), limbs(l) {
    trim();
}

BigInt::BigInt(uint64_t value) : is_negative(false) {
    limbs.push_back(value);
}

BigInt::BigInt(const std::string& str) : is_negative(false) {
    if (str.empty()) {
        limbs.push_back(0);
        return;
    }
    
    size_t i = 0;
    if (str[i] == '-') {
        is_negative = true;
        i++;
    } else if (str[i] == '+') {
        i++;
    }
    
    limbs.push_back(0);
    // Base 10 parse
    for (; i < str.length(); ++i) {
        if (str[i] < '0' || str[i] > '9') {
            throw std::invalid_argument("Invalid character in BigInt string");
        }
        
        // Multiply current by 10 and add digit
        BigInt ten(10ULL);
        *this *= ten;
        BigInt digit(static_cast<uint64_t>(str[i] - '0'));
        if (is_negative) {
            digit.is_negative = true;
            *this += digit; // Actually subtracts abs value if negative, handled by += logic
        } else {
            *this += digit;
        }
    }
    trim();
}

void BigInt::trim() {
    while (limbs.size() > 1 && limbs.back() == 0) {
        limbs.pop_back();
    }
    if (limbs.size() == 1 && limbs[0] == 0) {
        is_negative = false;
    }
}

int BigInt::compare_abs(const BigInt& other) const {
    if (limbs.size() != other.limbs.size()) {
        return limbs.size() < other.limbs.size() ? -1 : 1;
    }
    for (size_t i = limbs.size(); i-- > 0; ) {
        if (limbs[i] != other.limbs[i]) {
            return limbs[i] < other.limbs[i] ? -1 : 1;
        }
    }
    return 0;
}

int BigInt::compare(const BigInt& other) const {
    if (is_negative != other.is_negative) {
        return is_negative ? -1 : 1;
    }
    int cmp = compare_abs(other);
    return is_negative ? -cmp : cmp;
}

void BigInt::add_limbs(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    size_t n = std::max(a.size(), b.size());
    result.resize(n);
    unsigned char carry = 0;
    
    for (size_t i = 0; i < n; ++i) {
        uint64_t a_val = (i < a.size()) ? a[i] : 0;
        uint64_t b_val = (i < b.size()) ? b[i] : 0;
        carry = _addcarry_u64(carry, a_val, b_val, (unsigned long long*)&result[i]);
    }
    if (carry) {
        result.push_back(carry);
    }
}

void BigInt::sub_limbs(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    // Requires a >= b
    size_t n = a.size();
    result.resize(n);
    unsigned char borrow = 0;
    
    for (size_t i = 0; i < n; ++i) {
        uint64_t a_val = a[i];
        uint64_t b_val = (i < b.size()) ? b[i] : 0;
        borrow = _subborrow_u64(borrow, a_val, b_val, (unsigned long long*)&result[i]);
    }
}

// naive multiplication
// uses built-in 128 bit ints to avoid carry bit drops.
void BigInt::mul_limbs_naive(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    result.assign(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.size(); ++j) {
            unsigned __int128 p = (unsigned __int128)a[i] * b[j] + result[i+j] + carry;
            result[i+j] = (uint64_t)p;
            carry = (uint64_t)(p >> 64);
        }
        result[i + b.size()] = carry;
    }
}

// straightforward o(n^2) division, useful fallback for tiny divisors
uint64_t BigInt::div_small(uint64_t divisor) {
    if (divisor == 0) throw std::domain_error("Division by zero");
    uint64_t rem = 0;
    for (size_t i = limbs.size(); i-- > 0; ) {
        unsigned __int128 n = ((unsigned __int128)rem << 64) | limbs[i];
        limbs[i] = (uint64_t)(n / divisor);
        rem = (uint64_t)(n % divisor);
    }
    trim();
    return rem;
}

void BigInt::shift_left_limbs(uint64_t n) {
    if (n == 0) return;
    if (limbs.size() == 1 && limbs[0] == 0) return;
    limbs.resize(limbs.size() + n);
    for (size_t i = limbs.size(); i-- > n; ) {
        limbs[i] = limbs[i - n];
    }
    for (size_t i = 0; i < n; ++i) limbs[i] = 0;
}

void BigInt::shift_right_limbs(uint64_t n) {
    if (n == 0) return;
    if (n >= limbs.size()) {
        limbs = {0};
    } else {
        for (size_t i = 0; i < limbs.size() - n; ++i) {
            limbs[i] = limbs[i + n];
        }
        limbs.resize(limbs.size() - n);
    }
}

void BigInt::shift_left_bits(uint32_t n) {
    if (n == 0) return;
    if (limbs.size() == 1 && limbs[0] == 0) return;
    uint32_t limbs_shift = n / 64;
    uint32_t bits_shift = n % 64;
    shift_left_limbs(limbs_shift);
    if (bits_shift == 0) return;
    
    uint64_t carry = 0;
    for (size_t i = 0; i < limbs.size(); ++i) {
        uint64_t next_carry = limbs[i] >> (64 - bits_shift);
        limbs[i] = (limbs[i] << bits_shift) | carry;
        carry = next_carry;
    }
    if (carry) {
        limbs.push_back(carry);
    }
}

BigInt BigInt::divide_long(const BigInt& numerator, const BigInt& denominator) {
    if (denominator.limbs.size() == 1 && denominator.limbs[0] == 0) throw std::domain_error("Division by zero");
    if (numerator.compare_abs(denominator) < 0) return BigInt(0ULL);

    if (denominator.limbs.size() == 1) {
        BigInt q = numerator;
        q.div_small(denominator.limbs[0]);
        return q;
    }

    uint32_t shift = __builtin_clzll(denominator.limbs.back());
    BigInt num = numerator;
    BigInt den = denominator;
    num.shift_left_bits(shift);
    den.shift_left_bits(shift);

    BigInt q, r;
    q.limbs.resize(num.limbs.size());
    
    for (size_t i = num.limbs.size(); i-- > 0; ) {
        r.shift_left_limbs(1);
        r.limbs[0] = num.limbs[i];
        r.trim();
        
        if (r.compare_abs(den) < 0) {
            q.limbs[i] = 0;
            continue;
        }

        uint64_t d_top = den.limbs.back();
        unsigned __int128 r_top = r.limbs.back();
        if (r.limbs.size() > den.limbs.size()) {
            r_top = (r_top << 64) | r.limbs[r.limbs.size() - 2];
        } else if (r.limbs.size() == den.limbs.size()) {
            // Need to make sure r_top uses previous limb if valid, or just current.
            // If sizes are equal, r_top doesn't have a higher limb than den's top.
            // But we might need the 2-limb view anyway.
        }
        
        uint64_t q_guess;
        if (r.limbs.size() > den.limbs.size() && r.limbs.back() == d_top) {
            q_guess = 0xFFFFFFFFFFFFFFFFULL;
        } else {
            q_guess = (uint64_t)(r_top / d_top);
        }
        
        BigInt test = den * BigInt(q_guess);
        
        while (test.compare_abs(r) > 0) {
            q_guess--;
            test -= den;
        }
        
        q.limbs[i] = q_guess;
        r -= test;
    }
    q.trim();
    return q;
}

BigInt BigInt::divide(const BigInt& numerator, const BigInt& denominator) {
    if (denominator.limbs.size() == 1 && denominator.limbs[0] == 0) throw std::domain_error("Division by zero");
    if (numerator.compare_abs(denominator) < 0) return BigInt(0ULL);
    
    if (denominator.limbs.size() < 128) {
        return divide_long(numerator, denominator);
    }
    
    BigInt abs_num = numerator; abs_num.is_negative = false;
    BigInt abs_den = denominator; abs_den.is_negative = false;
    
    size_t n = abs_num.limbs.size();
    size_t m = abs_den.limbs.size();
    size_t K = n + 1;
    
    size_t top_limbs = std::min((size_t)64, m);
    // newton raphson magic happens here
    // X_{k+1} = X_k * (2 - D * X_k)
    BigInt B_top;
    B_top.limbs = std::vector<uint64_t>(abs_den.limbs.begin() + (abs_den.limbs.size() - top_limbs), abs_den.limbs.end());
    
    size_t k_top = top_limbs * 2;
    BigInt two_ktop;
    two_ktop.limbs.assign(k_top + 1, 0);
    two_ktop.limbs.back() = 1;
    
    BigInt X = divide_long(two_ktop, B_top);
    
    if (K >= m + top_limbs) {
        X.shift_left_limbs(K - m - top_limbs);
    } else {
        X.shift_right_limbs(m + top_limbs - K);
    }
    
    BigInt two_K;
    two_K.limbs.assign(K + 1, 0);
    two_K.limbs.back() = 1;
    
    while (true) {
        BigInt B_X = abs_den * X;
        if (B_X == two_K) break;
        
        BigInt err;
        bool err_neg = false;
        if (two_K >= B_X) {
            err = two_K - B_X;
        } else {
            err = B_X - two_K;
            err_neg = true;
        }
        
        BigInt adj = X * err;
        adj.shift_right_limbs(K);
        
        if (adj.limbs.size() == 1 && adj.limbs[0] == 0) break;
        
        if (err_neg) X -= adj;
        else X += adj;
    }
    
    BigInt Q = abs_num * X;
    Q.shift_right_limbs(K);
    
    while ((Q + BigInt(1ULL)) * abs_den <= abs_num) {
        Q += BigInt(1ULL);
    }
    while (Q * abs_den > abs_num) {
        Q -= BigInt(1ULL);
    }
    
    if (numerator.is_negative != denominator.is_negative) Q.is_negative = true;
    Q.trim();
    return Q;
}

// helper doing the recursive chunking
static std::string convert_rec(const BigInt& x, const std::vector<BigInt>& p10, int k, bool pad) {
    if (k < 7) { 
        std::string s = x.to_string();
        if (pad) {
            size_t target_len = 1ULL << k;
            if (s.length() < target_len) {
                s = std::string(target_len - s.length(), '0') + s;
            }
        }
        return s;
    }
    
    if (x < p10[k-1]) {
        std::string s = convert_rec(x, p10, k-1, false);
        if (pad) {
            size_t target_len = 1ULL << k;
            if (s.length() < target_len) {
                s = std::string(target_len - s.length(), '0') + s;
            }
        }
        return s;
    }
    
    BigInt q = BigInt::divide(x, p10[k-1]);
    BigInt r = x - q * p10[k-1];
    
    std::string sq = convert_rec(q, p10, k-1, pad);
    std::string sr = convert_rec(r, p10, k-1, true);
    
    return sq + sr;
}

BigInt BigInt::power(uint64_t base_val, uint64_t exp) {
    BigInt res(1ULL);
    BigInt b(base_val);
    while (exp > 0) {
        if (exp % 2 == 1) res *= b;
        b *= b;
        exp /= 2;
    }
    return res;
}

std::string BigInt::fast_to_string() const {
    if (limbs.size() == 1 && limbs[0] == 0) return "0";
    
    BigInt temp = *this;
    bool is_negative = temp.is_negative;
    temp.is_negative = false;

    // generate base powers of 10 scaling up to the size of the number
    std::vector<BigInt> p10;
    p10.push_back(BigInt(10ULL)); 
    int k = 0;
    while (p10.back() <= temp) {
        p10.push_back(p10.back() * p10.back());
        k++;
    }
    
    std::string s = convert_rec(temp, p10, k, false);
    if (is_negative) s = "-" + s;
    return s;
}

BigInt& BigInt::operator+=(const BigInt& other) {
    if (is_negative == other.is_negative) {
        std::vector<uint64_t> res;
        add_limbs(res, limbs, other.limbs);
        limbs = std::move(res);
    } else {
        int cmp = compare_abs(other);
        if (cmp >= 0) {
            std::vector<uint64_t> res;
            sub_limbs(res, limbs, other.limbs);
            limbs = std::move(res);
        } else {
            std::vector<uint64_t> res;
            sub_limbs(res, other.limbs, limbs);
            limbs = std::move(res);
            is_negative = other.is_negative;
        }
    }
    trim();
    return *this;
}

BigInt& BigInt::operator-=(const BigInt& other) {
    if (is_negative != other.is_negative) {
        std::vector<uint64_t> res;
        add_limbs(res, limbs, other.limbs);
        limbs = std::move(res);
    } else {
        int cmp = compare_abs(other);
        if (cmp >= 0) {
            std::vector<uint64_t> res;
            sub_limbs(res, limbs, other.limbs);
            limbs = std::move(res);
        } else {
            std::vector<uint64_t> res;
            sub_limbs(res, other.limbs, limbs);
            limbs = std::move(res);
            is_negative = !is_negative;
        }
    }
    trim();
    return *this;
}

void BigInt::mul_limbs_karatsuba(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    size_t n = a.size();
    size_t m = b.size();
    
    // Threshold to switch to naive multiplication
    if (n < 64 || m < 64) {
        mul_limbs_naive(result, a, b);
        return;
    }
    
    size_t k = std::max(n, m) / 2;
    
    // Split: a = a_h * 2^(64k) + a_l, b = b_h * 2^(64k) + b_l
    std::vector<uint64_t> a_l, a_h, b_l, b_h;
    if (n > k) {
        a_l.assign(a.begin(), a.begin() + k);
        a_h.assign(a.begin() + k, a.end());
    } else {
        a_l = a;
        a_h = {0};
    }
    
    if (m > k) {
        b_l.assign(b.begin(), b.begin() + k);
        b_h.assign(b.begin() + k, b.end());
    } else {
        b_l = b;
        b_h = {0};
    }
    
    // 3 multiplications
    std::vector<uint64_t> z0, z1, z2;
    mul_limbs_karatsuba(z0, a_l, b_l); // z0 = a_l * b_l
    mul_limbs_karatsuba(z2, a_h, b_h); // z2 = a_h * b_h
    
    std::vector<uint64_t> sum_a, sum_b;
    add_limbs(sum_a, a_l, a_h);
    add_limbs(sum_b, b_l, b_h);
    mul_limbs_karatsuba(z1, sum_a, sum_b); // z1 = (a_l+a_h)*(b_l+b_h)
    
    // Middle term: z1 - z2 - z0
    // We can do this with vector arithmetic
    // res = z2 * 2^(128k) + (z1 - z2 - z0) * 2^(64k) + z0
    
    // ... For simplicity and safety, let's use BigInt wrappers for these vector ops
    // but the final version will be optimized vector-only.
    // For now, let's just make z1 = z1 - z2 - z0
    BigInt b_z0; b_z0.limbs = z0;
    BigInt b_z1; b_z1.limbs = z1;
    BigInt b_z2; b_z2.limbs = z2;
    
    b_z1 -= b_z2;
    b_z1 -= b_z0;
    
    // Shift and add
    // b_z2 << 2k, b_z1 << k
    BigInt final_res = b_z2;
    final_res.limbs.insert(final_res.limbs.begin(), 2 * k, 0);
    
    BigInt mid = b_z1;
    mid.limbs.insert(mid.limbs.begin(), k, 0);
    
    final_res += mid;
    final_res += b_z0;
    
    result = final_res.limbs;
}

void ntt_multiply_adapter(std::vector<uint64_t>& res, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    // Convert 64-bit limbs to 32-bit limbs to avoid overflow in NTT convolution
    std::vector<uint64_t> a32, b32;
    for (uint64_t v : a) {
        a32.push_back(v & 0xFFFFFFFF);
        a32.push_back(v >> 32);
    }
    for (uint64_t v : b) {
        b32.push_back(v & 0xFFFFFFFF);
        b32.push_back(v >> 32);
    }

    std::vector<__uint128_t> r32;
    NTT::multiply_crt(r32, a32, b32);

    // Carry and reconstruct 64-bit limbs
    unsigned __int128 carry = 0;
    std::vector<uint64_t> result32;
    for (size_t i = 0; i < r32.size(); ++i) {
        unsigned __int128 val = r32[i] + carry;
        result32.push_back((uint64_t)(val & 0xFFFFFFFF));
        carry = val >> 32;
    }
    while (carry) {
        result32.push_back((uint64_t)(carry & 0xFFFFFFFF));
        carry >>= 32;
    }

    res.clear();
    for (size_t i = 0; i < result32.size(); i += 2) {
        uint64_t low = result32[i];
        uint64_t high = (i + 1 < result32.size()) ? result32[i + 1] : 0;
        res.push_back(low | (high << 32));
    }
}

BigInt& BigInt::operator*=(const BigInt& other) {
    if ((limbs.size() == 1 && limbs[0] == 0) || (other.limbs.size() == 1 && other.limbs[0] == 0)) {
        limbs = {0};
        is_negative = false;
        return *this;
    }

    std::vector<uint64_t> res;
    if (limbs.size() > 512 || other.limbs.size() > 512) {
        ntt_multiply_adapter(res, limbs, other.limbs);
    } else if (limbs.size() >= 64 && other.limbs.size() >= 64) {
        mul_limbs_karatsuba(res, limbs, other.limbs);
    } else {
        mul_limbs_naive(res, limbs, other.limbs);
    }
    limbs = std::move(res);
    is_negative = is_negative != other.is_negative;
    trim();
    return *this;
}

BigInt operator+(BigInt lhs, const BigInt& rhs) { lhs += rhs; return lhs; }
BigInt operator-(BigInt lhs, const BigInt& rhs) { lhs -= rhs; return lhs; }
BigInt operator*(BigInt lhs, const BigInt& rhs) { lhs *= rhs; return lhs; }

bool operator==(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) == 0; }
bool operator!=(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) != 0; }
bool operator<(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) < 0; }
bool operator<=(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) <= 0; }
bool operator>(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) > 0; }
bool operator>=(const BigInt& lhs, const BigInt& rhs) { return lhs.compare(rhs) >= 0; }

std::string BigInt::to_string() const {
    if (limbs.size() == 1 && limbs[0] == 0) return "0";
    
    BigInt temp = *this;
    temp.is_negative = false;
    std::string res = "";
    
    // Naive base 10 conversion (Very slow for large numbers, will optimize later)
    BigInt ten(10ULL);
    while (temp > BigInt(0ULL)) {
        // Implement division by 10
        // For now, this is missing so to_string will fail. We need div/mod by a small uint64_t.
        // Let's implement an in-place div_mod_10
        uint64_t rem = 0;
        for (size_t i = temp.limbs.size(); i-- > 0; ) {
            unsigned __int128 n = ((unsigned __int128)rem << 64) | temp.limbs[i];
            temp.limbs[i] = (uint64_t)(n / 10);
            rem = (uint64_t)(n % 10);
        }
        temp.trim();
        res += (char)('0' + rem);
    }
    
    if (is_negative) res += "-";
    std::reverse(res.begin(), res.end());
    return res;
}
