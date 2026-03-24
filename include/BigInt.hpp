#pragma once

#include <vector>
#include <cstdint>
#include <string>


class BigInt {
public:
    std::vector<uint64_t> limbs; // little-endian base 2^64
    bool is_negative;

    // generic constructors
    // generic constructors
    BigInt();
    BigInt(uint64_t val);
    BigInt(const std::string& str);
    BigInt(const std::vector<uint64_t>& l);
    BigInt(const BigInt& other) = default;

    // assignment operators and core math
    BigInt& operator=(const BigInt& other) = default;
    
    BigInt& operator+=(const BigInt& other);
    BigInt& operator-=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);
    
    // newton's method division (O(M(N)) speed)
    // divide_long is the base case fallback using regular division
    uint64_t div_small(uint64_t divisor); 
    static BigInt divide_long(const BigInt& numerator, const BigInt& denominator);
    static BigInt divide(const BigInt& numerator, const BigInt& denominator);
    
    // Fast exponentiation
    static BigInt power(uint64_t base, uint64_t exp);

    // Core low-level routines
    static void add_limbs(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);
    static void sub_limbs(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);
    static void mul_limbs_naive(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);
    static void mul_limbs_karatsuba(std::vector<uint64_t>& result, const std::vector<uint64_t>& a, const std::vector<uint64_t>& b);

    // helper functions
    void trim(); 
    void shift_left_limbs(uint64_t n);  // moves the whole chunks
    void shift_right_limbs(uint64_t n); // drops lowest chunks
    void shift_left_bits(uint32_t n);   // shifts individual bits
    int compare_abs(const BigInt& other) const;
    int compare(const BigInt& other) const;
    std::string to_string() const; 
    
    // ridiculously fast base 10 conversion
    std::string fast_to_string() const;

    friend BigInt operator+(BigInt lhs, const BigInt& rhs);
    friend BigInt operator-(BigInt lhs, const BigInt& rhs);
    friend BigInt operator*(BigInt lhs, const BigInt& rhs);
    
    friend bool operator==(const BigInt& lhs, const BigInt& rhs);
    friend bool operator!=(const BigInt& lhs, const BigInt& rhs);
    friend bool operator<(const BigInt& lhs, const BigInt& rhs);
    friend bool operator<=(const BigInt& lhs, const BigInt& rhs);
    friend bool operator>(const BigInt& lhs, const BigInt& rhs);
    friend bool operator>=(const BigInt& lhs, const BigInt& rhs);
};
