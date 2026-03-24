#include "BigInt.hpp"
#include <iostream>
#include <cassert>

void test_addition() {
    BigInt a("12345678901234567890");
    BigInt b("98765432109876543210");
    BigInt c = a + b;
    assert(c.to_string() == "111111111011111111100");
    std::cout << "test_addition passed.\n";
}

void test_subtraction() {
    BigInt a("100000000000000000000");
    BigInt b("1");
    BigInt c = a - b;
    assert(c.to_string() == "99999999999999999999");
    
    BigInt d("10");
    BigInt e("20");
    BigInt f = d - e;
    assert(f.to_string() == "-10");
    std::cout << "test_subtraction passed.\n";
}

void test_multiplication() {
    BigInt a("999999999");
    BigInt b("999999999");
    BigInt c = a * b;
    assert(c.to_string() == "999999998000000001");
    std::cout << "test_multiplication passed.\n";
}

void test_boundary_overflow() {
    BigInt a(0xFFFFFFFFFFFFFFFFULL);
    BigInt b(1);
    BigInt c = a + b;
    // 2^64
    assert(c.to_string() == "18446744073709551616");
    std::cout << "test_boundary_overflow passed.\n";
}

void test_division() {
    BigInt a("10000000000000000000000000000000000000000");
    BigInt b("3");
    BigInt c = BigInt::divide(a, b);
    assert(c.to_string() == "3333333333333333333333333333333333333333");
    
    BigInt x("123456789123456789123456789");
    BigInt y("987654321");
    BigInt z = BigInt::divide(x, y);
    assert(z.to_string() == "124999998985937499");
    std::cout << "test_division passed.\n";
}

void test_fast_string() {
    BigInt a("99999999999999999999999999999999999999999999999999999999999999999999999999999999");
    std::string s = a.fast_to_string();
    if (s != "99999999999999999999999999999999999999999999999999999999999999999999999999999999") {
        std::cout << "FAST_STRING_IS: " << s << "\n";
        assert(false);
    }
    std::cout << "test_fast_string passed.\n";
}

int main() {
    test_addition();
    test_subtraction();
    test_multiplication();
    test_boundary_overflow();
    test_division();
    test_fast_string();
    std::cout << "All tests passed!\n";
    return 0;
}
