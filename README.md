# pi-from-scratch

wrote a custom c++ pi calculator to see if i could compute 200 million digits without using gmp or any big math libraries lol. just pure c++20 and pain.

## how it works
basically, calculating pi to 200m digits with normal math is impossible because it scales at O(N^2). naive multiplication and division would literally take years on my computer.

so i had to implement some crazy stuff:
- **NTT (Number Theoretic Transform):** implemented a custom fast fourier transform over finite fields using four different prime moduli. had to use 128-bit ints for the chinese remainder theorem reconstruction because 64-bit kept overflowing and giving me garbage arrays.
- **Newton-Raphson division:** instead of doing normal long division, it guesses the reciprocal `1/X` and iteratively doubles the precision using the NTT multiplier. division is basically instant now.
- **recursive base 10 conversion:** converting a gigabyte-sized binary blob to a decimal string natively takes forever. i generate powers of 10 (`10^64`, `10^128` etc) and recursively split the array in half.

## how to run

```bash
mkdir build && cd build
cmake ..
make -j
```

to generate digits (e.g. 10000 digits):
```bash
./pi_calc 10000
```
the digits will dump straight into `pi_digits.txt`. 

there's also a test suite if you want to run it: `./test_bigint` (checks all the edge cases and string conversions so i know the math engine actually works).

anyway hope you like it!
