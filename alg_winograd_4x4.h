//
// Created by Max Shikunov on 18/12/2025.
//

#ifndef UNTITLED3_ALG_WINOGRAD_4X4_H
#define UNTITLED3_ALG_WINOGRAD_4X4_H

#include "structures.h"

template <class T>
void mul_winograd_4x4(const Matrix<T>& A,
                      const Matrix<T>& B,
                      Matrix<T>& C,
                      OpCounter* cnt=nullptr) {

    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);

    C.resize(4,4);

    T p[4], q[4];

    // p[i] = -A[i,0]*A[i,1] - A[i,2]*A[i,3]
    // q[j] = -B[0,j]*B[1,j] - B[2,j]*B[3,j]
    for (int i = 0; i < 4; ++i) {
        T t1 = mul(A(i,0), A(i,1), cnt);
        T t2 = mul(A(i,2), A(i,3), cnt);
        p[i] = sub(sub(T{}, t1, cnt), t2, cnt); // 0 - t1 - t2
    }
    for (int j = 0; j < 4; ++j) {
        T t1 = mul(B(0,j), B(1,j), cnt);
        T t2 = mul(B(2,j), B(3,j), cnt);
        q[j] = sub(sub(T{}, t1, cnt), t2, cnt);
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // (A[i,0]+B[1,j])*(A[i,1]+B[0,j])
            T s1 = add(A(i,0), B(1,j), cnt);
            T s2 = add(A(i,1), B(0,j), cnt);
            T m1 = mul(s1, s2, cnt);

            // (A[i,2]+B[3,j])*(A[i,3]+B[2,j])
            T s3 = add(A(i,2), B(3,j), cnt);
            T s4 = add(A(i,3), B(2,j), cnt);
            T m2 = mul(s3, s4, cnt);

            T res = add(add(add(p[i], q[j], cnt), m1, cnt), m2, cnt);
            C(i,j) = res;
        }
    }
}

#endif //UNTITLED3_ALG_WINOGRAD_4X4_H
