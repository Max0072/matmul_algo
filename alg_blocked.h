//
// Blocked matrix multiplication для масштабирования на большие размеры
// Использует 4x4 ядра (naive/winograd/alphaevolve) для блочного умножения
//

#ifndef ALG_BLOCKED_H
#define ALG_BLOCKED_H

#include "structures.h"
#include "alg_naive.h"
#include "alg_winograd_4x4.h"
#include "alg_alpha_evolve_4x4_complex.h"
#include "alg_strassen_4x4.h"

// Вспомогательная функция: умножение 4x4 блоков naive
template <class T>
void kernel_naive_4x4(MatrixView<const T> A, MatrixView<const T> B,
                      MatrixView<T> C, OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);
    assert(C.rows == 4 && C.cols == 4);

    mul_naive_view(A, B, C, cnt);
}

// Вспомогательная функция: умножение 4x4 блоков Winograd
template <class T>
void kernel_winograd_4x4(MatrixView<const T> A, MatrixView<const T> B,
                         MatrixView<T> C, OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);
    assert(C.rows == 4 && C.cols == 4);


    // Пока используем временные матрицы
    Matrix<T> Am(4,4), Bm(4,4), Cm(4,4);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++) {
            Am(i,j) = A(i,j);
            Bm(i,j) = B(i,j);
        }

    mul_winograd_4x4(Am, Bm, Cm, cnt);

    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            C(i,j) = Cm(i,j);
}

// Вспомогательная функция: умножение 4x4 блоков AlphaEvolve
template <class T>
void kernel_alphaevolve_4x4(MatrixView<const T> A, MatrixView<const T> B,
                            MatrixView<T> C, OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);
    assert(C.rows == 4 && C.cols == 4);


    // Пока используем временные матрицы
    Matrix<T> Am(4,4), Bm(4,4), Cm(4,4);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++) {
            Am(i,j) = A(i,j);
            Bm(i,j) = B(i,j);
        }

    alphaevolve_4x4_complex(Am, Bm, Cm, cnt);

    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            C(i,j) = Cm(i,j);
}

// Вспомогательная функция: умножение 4x4 блоков Strassen
template <class T>
void kernel_strassen_4x4(MatrixView<const T> A, MatrixView<const T> B,
                         MatrixView<T> C, OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);
    assert(C.rows == 4 && C.cols == 4);

    mul_strassen_4x4_view(A, B, C, cnt);
}

// Тип ядра для блочного умножения
enum class BlockKernel {
    NAIVE,
    WINOGRAD,
    ALPHAEVOLVE,
    STRASSEN
};

// Blocked multiply: делит матрицу на блоки 4x4 и умножает их выбранным ядром
template <class T>
void mul_blocked(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                 BlockKernel kernel = BlockKernel::NAIVE,
                 OpCounter* cnt = nullptr) {

    assert(A.cols == B.rows);

    int m = A.rows;  // строки A
    int k = A.cols;  // столбцы A = строки B
    int n = B.cols;  // столбцы B

    C.resize(m, n);

    // Обнуляем результат
    for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
            C(i,j) = T{};

    // Размер блока
    const int BS = 4;

    // Количество блоков по каждому измерению
    int num_blocks_m = (m + BS - 1) / BS;  // количество блоков по строкам A/C
    int num_blocks_k = (k + BS - 1) / BS;  // количество блоков по столбцам A / строкам B
    int num_blocks_n = (n + BS - 1) / BS;  // количество блоков по столбцам B/C

    // Создаем views для работы с подматрицами
    auto A_view = view(A);
    auto B_view = view(B);
    auto C_view = view(C);

    // Блочное умножение: C = A * B
    // C[bi, bj] = sum_bp (A[bi, bp] * B[bp, bj])
    for (int bi = 0; bi < num_blocks_m; bi++) {
        for (int bj = 0; bj < num_blocks_n; bj++) {
            // Размеры текущего блока C[bi, bj]
            int c_row_start = bi * BS;
            int c_col_start = bj * BS;
            int c_rows = std::min(BS, m - c_row_start);
            int c_cols = std::min(BS, n - c_col_start);

            // Итерация по промежуточным блокам
            for (int bp = 0; bp < num_blocks_k; bp++) {
                // Размеры блока A[bi, bp]
                int a_row_start = bi * BS;
                int a_col_start = bp * BS;
                int a_rows = std::min(BS, m - a_row_start);
                int a_cols = std::min(BS, k - a_col_start);

                // Размеры блока B[bp, bj]
                int b_row_start = bp * BS;
                int b_col_start = bj * BS;
                int b_rows = std::min(BS, k - b_row_start);
                int b_cols = std::min(BS, n - b_col_start);

                // Извлекаем subviews
                auto A_block = subview(A_view, a_row_start, a_col_start, a_rows, a_cols);
                auto B_block = subview(B_view, b_row_start, b_col_start, b_rows, b_cols);
                auto C_block = subview(C_view, c_row_start, c_col_start, c_rows, c_cols);

                // Временная матрица для результата умножения блоков
                Matrix<T> temp_result(c_rows, c_cols);
                for(int ti=0; ti<c_rows; ti++)
                    for(int tj=0; tj<c_cols; tj++)
                        temp_result(ti, tj) = T{};

                auto temp_view = view(temp_result);

                // Выбираем ядро и умножаем блоки
                if (a_rows == BS && a_cols == BS && b_rows == BS && b_cols == BS && c_rows == BS && c_cols == BS) {
                    // Полные блоки 4x4 - используем выбранное ядро
                    switch(kernel) {
                        case BlockKernel::NAIVE:
                            kernel_naive_4x4(A_block, B_block, temp_view, cnt);
                            break;
                        case BlockKernel::WINOGRAD:
                            kernel_winograd_4x4(A_block, B_block, temp_view, cnt);
                            break;
                        case BlockKernel::ALPHAEVOLVE:
                            kernel_alphaevolve_4x4(A_block, B_block, temp_view, cnt);
                            break;
                        case BlockKernel::STRASSEN:
                            kernel_strassen_4x4(A_block, B_block, temp_view, cnt);
                            break;
                    }
                } else {
                    // Граничные блоки - используем naive
                    mul_naive_view(A_block, B_block, temp_view, cnt);
                }

                // Добавляем результат к блоку C[bi, bj]
                for(int ti=0; ti<c_rows; ti++) {
                    for(int tj=0; tj<c_cols; tj++) {
                        C_block(ti, tj) = add(C_block(ti, tj), temp_result(ti, tj), cnt);
                    }
                }
            }
        }
    }
}

// Wrapper функции для удобного вызова
template <class T>
void mul_blocked_naive_kernel(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                              OpCounter* cnt = nullptr) {
    mul_blocked(A, B, C, BlockKernel::NAIVE, cnt);
}

template <class T>
void mul_blocked_winograd_kernel(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                                 OpCounter* cnt = nullptr) {
    mul_blocked(A, B, C, BlockKernel::WINOGRAD, cnt);
}

template <class T>
void mul_blocked_alphaevolve_kernel(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                                    OpCounter* cnt = nullptr) {
    mul_blocked(A, B, C, BlockKernel::ALPHAEVOLVE, cnt);
}

template <class T>
void mul_blocked_strassen_kernel(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                                 OpCounter* cnt = nullptr) {
    mul_blocked(A, B, C, BlockKernel::STRASSEN, cnt);
}

#endif // ALG_BLOCKED_H
