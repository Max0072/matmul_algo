//
// Created by Max Shikunov on 18/12/2025.
//

#ifndef UNTITLED3_ALG_NAIVE_H
#define UNTITLED3_ALG_NAIVE_H

#include "structures.h"

template <class T>
void mul_naive_view(MatrixView<const T>& A,
                    MatrixView<const T>& B,
                    MatrixView<T>& C,
                    OpCounter* cnt = nullptr) {

    assert(A.cols == B.rows);
    assert(A.rows == C.rows and B.cols == C.cols);

    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            T sum = T{};
            for (int k = 0; k < A.cols; ++k) {
                sum = add(sum,mul(A(i, k), B(k, j), cnt),cnt);
            }
            C(i, j) = sum;
        }
    }
}

template <class T>
void mul_naive(Matrix<T>& Am,
               Matrix<T>& Bm,
               Matrix<T>& Cm,
               OpCounter* cnt = nullptr) {

    Cm.resize(Am.rows, Bm.cols);;
    auto A = view(Am);
    auto B = view(Bm);
    auto C = view(Cm);

    assert(A.cols == B.rows);
    assert(A.rows == C.rows and B.cols == C.cols);

    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            T sum = T{};
            for (int k = 0; k < A.cols; ++k) {
                sum = add(sum,mul(A(i, k), B(k, j), cnt),cnt);
            }
            C(i, j) = sum;
        }
    }
}

#endif //UNTITLED3_ALG_NAIVE_H
