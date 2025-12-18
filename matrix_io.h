//
// Created by Max Shikunov on 18/12/2025.
//

#ifndef UNTITLED3_MATRIX_IO_H
#define UNTITLED3_MATRIX_IO_H

#include "structures.h"

#include "structures.h"

template <class T>
void read_matrix(Matrix<T>& A) {
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            std::cin >> A(i, j);
        }
    }
}

template <class T>
void print_matrix(const Matrix<T>& A) {
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            std::cout << A(i, j) << ' ';
        }
        std::cout << '\n';
    }
}


#endif //UNTITLED3_MATRIX_IO_H
