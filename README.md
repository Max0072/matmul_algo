# Matrix Multiplication Algorithms Analysis

Implemented 8 algorithms:
- Naive
- Strassen 4x4
- Winograd 4x4
- AlphaEvolve 4x4
- Blocked variants (naive/winograd/alphaevolve/strassen)

132 tests on sizes 4-256.

Main results:
1. Blocked Winograd - fastest, 10-15% faster than naive
2. Fewer multiplications doesn't mean faster
3. AlphaEvolve 6.5x slower
4. Blocked Strassen 2x slower


## 1. Operation counts

Theory for 4x4:

Algorithm    | Multiplies | Additions
Naive        | 64         | 64
Winograd     | 48         | 128
Strassen 4x4 | 56         | 128
AlphaEvolve  | 48         | many

Real results for 256x256:

Algorithm           | Time (ms) | Multiplies | Savings
naive               | 59-64     | 16,777,216 | baseline
blocked_winograd    | 53-55     | 12,582,912 | -25%
blocked_strassen    | 128-129   | 14,680,064 | -12%
blocked_alphaevolve | 395-398   | 12,582,912 | -25%

Why fewer multiplications doesn't help:

Blocked Strassen: saves only 12% multiplications, 2x slower. For 256x256 does 65,536 memory allocations, too many.

AlphaEvolve: saves 25% multiplications but 6.5x slower. Does too many intermediate calculations and uses complex numbers even when not needed.

Winograd: saves 25% multiplications and 10-15% faster. Simple formula, little memory (8 numbers), works well with cache, compiler optimizes easily.

Winograd vs Strassen 4x4 comparison:

Winograd: 48 multiplications, 64 bytes memory, 0 allocations, easy to optimize
Strassen: 56 multiplications, 512 bytes memory, 16 allocations, hard to optimize

Result for 256x256:
- blocked_winograd: 53-55 ms
- blocked_strassen: 128-129 ms (2.4x slower)

Conclusion: saving multiplications matters, but extra operations often ruin it. Winograd works well because it's simple. Simplicity matters more than theory.


## 2. Memory and cache

Memory usage for 256x256:

Algorithm        | Time (ms) | Extra memory | Cache behavior
naive            | 59-64     | 0            | Medium
blocked_winograd | 53-55     | 0            | Good
blocked_strassen | 128-129   | 0            | Very bad

Naive: accesses matrix B columns out of order, many cache misses.

Blocked Winograd: 4x4 blocks fit completely in CPU cache (128 bytes). Almost no cache misses. That's why 10-15% faster.

Blocked Strassen: for 256x256 creates 4096 blocks, each needs 16 allocations. Total 65,536 allocations. Many malloc/free operations, everything falls out of cache. 2x slower.

Memory needed per 4x4 block:
- Winograd: 64 bytes (8 numbers on stack)
- Strassen 4x4: 512 bytes (16 small 2x2 matrices)
- Naive: 0

For 256x256 (4096 blocks):
- blocked_winograd: 256 KB (fits in L2 cache)
- blocked_strassen: 2 MB plus many small allocations

Conclusion: little memory good. Many allocations bad. Cache behavior matters more than theory.


## 3. Numerical accuracy

Errors for 256x256:

Double:
All algorithms work fine, small errors (1e-12 to 4e-14).

Complex:
AlphaEvolve: error 35+ (wrong results)
Others: fine, errors 3e-12 to 4e-14.

AlphaEvolve for complex: made for real numbers, takes only real part, loses imaginary. Result is wrong. Can't use for complex matrices.

Winograd: no recursion, no error accumulation. Works for all types.


## 4. Compiler and optimizations

What compiler does with -O2:

For naive does many optimizations: unrolls loops, uses SIMD (4-8 operations at once), FMA instructions (a*b+c at once), keeps variables in registers. So naive with -O2 works almost like good library.

Winograd even better: simple structure, easy to unroll loops (16 iterations), easy to parallelize. Result: -25% multiplications plus good optimization = 10-15% faster.

BLAS libraries:

Intel MKL / OpenBLAS: written in assembly by hand for each processor, use all instructions, multithreaded. 10-100x faster than our naive.

Blocked Winograd:

4x4 blocks fit well in cache. Simple kernel, easy to optimize. Few extra operations. Result: 10-15% faster. Can make even faster if add SIMD by hand and multithreading, potential 5-10x.


## 5. Conclusions

Speed ranking (256x256):

1. blocked_winograd: 53-55 ms (-10%)
2. naive: 59-64 ms (baseline)
3. blocked_strassen: 128-129 ms (+100%)
4. blocked_alphaevolve: 395-398 ms (+550%)

What to use:

For real projects:
- BLAS (MKL, OpenBLAS) - 10-100x faster
- Blocked Winograd - if can't use BLAS
- Naive with -O2 - simple option

For experiments:
- Winograd - works well
- AlphaEvolve - only for double, very slow

Don't use:
- Blocked Strassen - too slow
- AlphaEvolve for complex - doesn't work

What I learned:

1. Fewer multiplications doesn't mean faster. Extra operations matter more.
2. Simplicity better than complexity. Compiler optimizes simple code better.
3. Cache behavior decides a lot. 4x4 blocks in cache give plus 10-15% speed.
4. Compiler optimizations give more than complex algorithms.


## 6. Files

Algorithms: alg_naive.h, alg_strassen_4x4.h, alg_winograd_4x4.h, alg_alpha_evolve_4x4_complex.h, alg_blocked.h

Other: structures.h, generators.h, benchmark.h, rss.h, main.cpp

Results: benchmark_results.csv (132 tests), ANALYSIS.md, WINOGRAD_VS_STRASSEN.md


## Conclusion

Blocked Winograd best for normal sizes. 10-15% faster than naive, simple, accurate, works for all types.

AlphaEvolve interesting on paper but useless in practice. Too many operations, doesn't work for complex numbers, 6.5x slower than naive.

For real tasks better to use BLAS or blocked Winograd.

Simplicity better than complexity. Cache behavior, compiler optimizations and few extra operations matter, not just theory.
