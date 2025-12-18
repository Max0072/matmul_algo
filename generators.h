//
// Created by Max Shikunov on 18/12/2025.
//

#ifndef UNTITLED3_GENERATORS_H
#define UNTITLED3_GENERATORS_H
#include "structures.h"
#include <random>
#include <type_traits>
#include <complex>

template <class T>
T rand_scalar(std::mt19937_64& rng, double lo, double hi) {
    std::uniform_real_distribution<double> dist(lo, hi);
    if constexpr (std::is_integral_v<T>) {
        std::uniform_int_distribution<long long> idist((long long)lo, (long long)hi);
        return (T)idist(rng);
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return {dist(rng), dist(rng)};
    } else {
        return (T)dist(rng);
    }
}

template <class T>
Matrix<T> gen_random(int r, int c, uint64_t seed=42, double lo=-1.0, double hi=1.0) {
    std::mt19937_64 rng(seed);
    Matrix<T> A(r,c);
    for (int i=0;i<r;++i)
        for (int j=0;j<c;++j)
            A(i,j)=rand_scalar<T>(rng, lo, hi);
    return A;
}

template <class T>
Matrix<T> gen_symmetric(int n, uint64_t seed=42, double lo=-1.0, double hi=1.0) {
    std::mt19937_64 rng(seed);
    Matrix<T> A(n,n);
    for (int i=0;i<n;++i) {
        for (int j=i;j<n;++j) {
            T x = rand_scalar<T>(rng, lo, hi);
            A(i,j)=x;
            A(j,i)=x;
        }
    }
    return A;
}

// почти разреженная: с вероятностью p_zero ставим 0
template <class T>
Matrix<T> gen_almost_sparse(int r, int c, double p_zero=0.9, uint64_t seed=42,
                            double lo=-1.0, double hi=1.0) {
    std::mt19937_64 rng(seed);
    std::bernoulli_distribution zero(p_zero);
    Matrix<T> A(r,c);
    for (int i=0;i<r;++i)
        for (int j=0;j<c;++j)
            A(i,j) = zero(rng) ? T{} : rand_scalar<T>(rng, lo, hi);
    return A;
}
#endif //UNTITLED3_GENERATORS_H
