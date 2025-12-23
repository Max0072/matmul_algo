//
// Created on 18/12/2025.
//

#ifndef UNTITLED3_ALG_STRASSEN_H
#define UNTITLED3_ALG_STRASSEN_H

#include "alg_naive.h"
#include "alg_strassen_4x4.h"
#include "structures.h"

inline int next_pow2(int n){
    int p=1;
    while(p<n) p<<=1;
    return p;
}

template <class T>
inline void mat_zero(MatrixView<T> C) {
    for (int i = 0; i < C.rows; ++i)
        for (int j = 0; j < C.cols; ++j)
            C(i, j) = T{};
}

template <class T>
inline void mat_copy(MatrixView<const T> A, MatrixView<T> B) {
    assert(A.rows == B.rows and A.cols == B.cols);
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            B(i, j) = A(i, j);
}

template <class T>
void mat_add(MatrixView<const T> A, MatrixView<const T> B, MatrixView<T> C, OpCounter* cnt=nullptr){
    assert(A.rows == B.rows && A.cols == B.cols);
    assert(C.rows == A.rows && C.cols == A.cols);
    for(int i = 0; i < A.rows; ++i)
        for(int j = 0; j < A.cols; ++j)
            C(i, j) = add(A(i, j),B(i, j),cnt);
}

template <class T>
void mat_sub(MatrixView<const T> A, MatrixView<const T> B, MatrixView<T> C, OpCounter* cnt=nullptr){
    assert(A.rows == B.rows && A.cols == B.cols);
    assert(C.rows == A.rows && C.cols == A.cols);
    for(int i = 0; i < A.rows; ++i)
        for(int j = 0; j < A.cols; ++j)
            C(i, j) = sub(A(i, j),B(i, j),cnt);
}

template <class T>
void mat_add_inplace(MatrixView<T> C, MatrixView<const T> A, OpCounter* cnt=nullptr){
    assert(C.rows == A.rows and C.cols == A.cols);
    for(int i = 0; i < C.rows; ++i)
        for(int j = 0; j < C.cols; ++j)
            C(i, j)=add(C(i,j), A(i,j),cnt);
}


template <class T>
inline void mat_sub_inplace(MatrixView<T> C, MatrixView<const T> A, OpCounter* cnt=nullptr) {
    assert(C.rows == A.rows and C.cols == A.cols);
    for (int i = 0; i < C.rows; ++i)
        for (int j = 0; j < C.cols; ++j)
            C(i,j) = sub(C(i,j), A(i,j), cnt);
}

template <class T>
void strassen_rec(MatrixView<const T> A,
                  MatrixView<const T> B,
                  MatrixView<T> C,
                  int n,
                  int THRESH,
                  OpCounter* cnt) {


    assert(A.rows == n and A.cols == n);
    assert(B.rows == n and B.cols == n);
    assert(C.rows == n and C.cols == n);

    // Базовый случай: используем специализированный Strassen 4x4 или naive
    if (n <= 4) {
        if (n == 4) {
            mul_strassen_4x4_view(A, B, C, cnt);
        } else {
            mul_naive_view(A, B, C, cnt);
        }
        return;
    }

    int h = n / 2;
    auto A11 = subview(A, 0, 0, h, h); // Matrix, point_y, point_x, length_y, length_x
    auto A12 = subview(A, 0, h, h, h);
    auto A21 = subview(A, h, 0, h, h);
    auto A22 = subview(A, h, h, h, h);

    auto B11 = subview(B, 0, 0, h, h);
    auto B12 = subview(B, 0, h, h, h);
    auto B21 = subview(B, h, 0, h, h);
    auto B22 = subview(B, h, h, h, h);

    auto C11 = subview(C, 0, 0, h, h);
    auto C12 = subview(C, 0, h, h, h);
    auto C21 = subview(C, h, 0, h, h);
    auto C22 = subview(C, h, h, h, h);

// temps: ONLY 3 matrices
    Matrix<T> S1m(h,h), S2m(h,h), Pm(h,h);
    auto S1 = view(S1m);
    auto S2 = view(S2m);
    auto P  = view(Pm);

// zero result blocks (important, because we accumulate)
    mat_zero(C11); mat_zero(C12); mat_zero(C21); mat_zero(C22);

// 1) P = M1 = (A11 + A22)(B11 + B22)
    mat_add(A11, A22, S1, cnt);
    mat_add(B11, B22, S2, cnt);
    strassen_rec<T>(view(S1m), view(S2m), P, h, THRESH, cnt);
    mat_add_inplace(C11, MatrixView<const T>(P), cnt);
    mat_add_inplace(C22, MatrixView<const T>(P), cnt);

// 2) P = M2 = (A21 + A22)B11
    mat_add(A21, A22, S1, cnt);
    strassen_rec<T>(view(S1m), B11, P, h, THRESH, cnt);
    mat_add_inplace(C21, MatrixView<const T>(P), cnt);
    mat_sub_inplace(C22, MatrixView<const T>(P), cnt);

// 3) P = M3 = A11(B12 - B22)
    mat_sub(B12, B22, S2, cnt);
    strassen_rec<T>(A11, view(S2m), P, h, THRESH, cnt);
    mat_add_inplace(C12, MatrixView<const T>(P), cnt);
    mat_add_inplace(C22, MatrixView<const T>(P), cnt);

// 4) P = M4 = A22(B21 - B11)
    mat_sub(B21, B11, S2, cnt);
    strassen_rec<T>(A22, view(S2m), P, h, THRESH, cnt);
    mat_add_inplace(C11, MatrixView<const T>(P), cnt);
    mat_add_inplace(C21, MatrixView<const T>(P), cnt);

// 5) P = M5 = (A11 + A12)B22
    mat_add(A11, A12, S1, cnt);
    strassen_rec<T>(view(S1m), B22, P, h, THRESH, cnt);
    mat_sub_inplace(C11, MatrixView<const T>(P), cnt);
    mat_add_inplace(C12, MatrixView<const T>(P), cnt);

// 6) P = M6 = (A21 - A11)(B11 + B12)
    mat_sub(A21, A11, S1, cnt);
    mat_add(B11, B12, S2, cnt);
    strassen_rec<T>(view(S1m), view(S2m), P, h, THRESH, cnt);
    mat_add_inplace(C22, MatrixView<const T>(P), cnt);

// 7) P = M7 = (A12 - A22)(B21 + B22)
    mat_sub(A12, A22, S1, cnt);
    mat_add(B21, B22, S2, cnt);
    strassen_rec<T>(view(S1m), view(S2m), P, h, THRESH, cnt);
    mat_add_inplace(C11, MatrixView<const T>(P), cnt);


}

// ---------- public wrapper: square multiply with padding ----------
template <class T>
void mul_strassen(const Matrix<T>& A,
                  const Matrix<T>& B,
                  Matrix<T>& C,
                  int THRESH = 64,
                  OpCounter* cnt=nullptr) {

    assert(A.cols == B.rows);

    int r = A.rows;
    int k = A.cols;
    int c = B.cols;

    int n  = std::max({r, k, c});
    int n2 = next_pow2(n);

    Matrix<T> Ap(n2, n2), Bp(n2, n2), Cp(n2, n2);

    mat_zero(view(Ap));
    mat_zero(view(Bp));
    mat_zero(view(Cp));

    // copy A into Ap
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < k; ++j)
            Ap(i, j) = A(i, j);

    // copy B into Bp
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < c; ++j)
            Bp(i, j) = B(i, j);

    strassen_rec<T>(view(Ap), view(Bp), view(Cp), n2, THRESH, cnt);

    // extract С
    C.resize(r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            C(i, j) = Cp(i, j);
}

#endif //UNTITLED3_ALG_STRASSEN_H
