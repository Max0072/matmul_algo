//
// Strassen 4x4 - специализированная версия без рекурсии
// 7 умножений 2×2 блоков вместо 8
//

#ifndef ALG_STRASSEN_4X4_H
#define ALG_STRASSEN_4X4_H

#include "structures.h"

// Вспомогательные функции для операций с 2×2 блоками

// Сложение 2×2 блоков: C = A + B
template<class T, class U, class V>
inline void add_2x2(U A, V B, MatrixView<T> C, OpCounter* cnt = nullptr) {
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
            C(i,j) = add(A(i,j), B(i,j), cnt);
}

// Вычитание 2×2 блоков: C = A - B
template<class T, class U, class V>
inline void sub_2x2(U A, V B, MatrixView<T> C, OpCounter* cnt = nullptr) {
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
            C(i,j) = sub(A(i,j), B(i,j), cnt);
}

// Наивное умножение 2×2 блоков: C = A * B
template<class T, class U, class V>
inline void mul_naive_2x2(U A, V B, MatrixView<T> C, OpCounter* cnt = nullptr) {
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            T sum = T{};
            for(int k = 0; k < 2; k++) {
                sum = add(sum, mul(A(i,k), B(k,j), cnt), cnt);
            }
            C(i,j) = sum;
        }
    }
}

// Strassen для 4×4 матрицы
// Разбивает на блоки 2×2 и применяет формулу Strassen
template<class T>
void mul_strassen_4x4(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
                      OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);

    C.resize(4, 4);

    auto A_view = view(A);
    auto B_view = view(B);
    auto C_view = view(C);

    // Разбиваем матрицы на блоки 2×2
    auto A11 = subview(A_view, 0, 0, 2, 2);
    auto A12 = subview(A_view, 0, 2, 2, 2);
    auto A21 = subview(A_view, 2, 0, 2, 2);
    auto A22 = subview(A_view, 2, 2, 2, 2);

    auto B11 = subview(B_view, 0, 0, 2, 2);
    auto B12 = subview(B_view, 0, 2, 2, 2);
    auto B21 = subview(B_view, 2, 0, 2, 2);
    auto B22 = subview(B_view, 2, 2, 2, 2);

    auto C11 = subview(C_view, 0, 0, 2, 2);
    auto C12 = subview(C_view, 0, 2, 2, 2);
    auto C21 = subview(C_view, 2, 0, 2, 2);
    auto C22 = subview(C_view, 2, 2, 2, 2);

    // Временные матрицы для промежуточных результатов
    Matrix<T> S1(2,2), S2(2,2), S3(2,2), S4(2,2), S5(2,2), S6(2,2), S7(2,2);
    Matrix<T> P1(2,2), P2(2,2), P3(2,2), P4(2,2), P5(2,2), P6(2,2), P7(2,2);
    Matrix<T> T1(2,2), T2(2,2);

    // Вычисляем 7 произведений Strassen
    // P1 = A11 * (B12 - B22)
    sub_2x2(B12, B22, view(S1), cnt);
    mul_naive_2x2(A11, view(S1), view(P1), cnt);

    // P2 = (A11 + A12) * B22
    add_2x2(A11, A12, view(S2), cnt);
    mul_naive_2x2(view(S2), B22, view(P2), cnt);

    // P3 = (A21 + A22) * B11
    add_2x2(A21, A22, view(S3), cnt);
    mul_naive_2x2(view(S3), B11, view(P3), cnt);

    // P4 = A22 * (B21 - B11)
    sub_2x2(B21, B11, view(S4), cnt);
    mul_naive_2x2(A22, view(S4), view(P4), cnt);

    // P5 = (A11 + A22) * (B11 + B22)
    add_2x2(A11, A22, view(S5), cnt);
    add_2x2(B11, B22, view(S6), cnt);
    mul_naive_2x2(view(S5), view(S6), view(P5), cnt);

    // P6 = (A12 - A22) * (B21 + B22)
    sub_2x2(A12, A22, view(S7), cnt);
    add_2x2(B21, B22, view(T1), cnt);
    mul_naive_2x2(view(S7), view(T1), view(P6), cnt);

    // P7 = (A11 - A21) * (B11 + B12)
    sub_2x2(A11, A21, view(T2), cnt);
    add_2x2(B11, B12, view(S1), cnt);
    mul_naive_2x2(view(T2), view(S1), view(P7), cnt);

    // Собираем результат
    // C11 = P5 + P4 - P2 + P6
    add_2x2(view(P5), view(P4), view(T1), cnt);
    sub_2x2(view(T1), view(P2), view(T2), cnt);
    add_2x2(view(T2), view(P6), C11, cnt);

    // C12 = P1 + P2
    add_2x2(view(P1), view(P2), C12, cnt);

    // C21 = P3 + P4
    add_2x2(view(P3), view(P4), C21, cnt);

    // C22 = P5 + P1 - P3 - P7
    add_2x2(view(P5), view(P1), view(T1), cnt);
    sub_2x2(view(T1), view(P3), view(T2), cnt);
    sub_2x2(view(T2), view(P7), C22, cnt);
}

// Версия с MatrixView для использования в blocked алгоритмах
template<class T>
void mul_strassen_4x4_view(MatrixView<const T> A, MatrixView<const T> B,
                           MatrixView<T> C, OpCounter* cnt = nullptr) {
    assert(A.rows == 4 && A.cols == 4);
    assert(B.rows == 4 && B.cols == 4);
    assert(C.rows == 4 && C.cols == 4);

    // Разбиваем матрицы на блоки 2×2
    auto A11 = subview(A, 0, 0, 2, 2);
    auto A12 = subview(A, 0, 2, 2, 2);
    auto A21 = subview(A, 2, 0, 2, 2);
    auto A22 = subview(A, 2, 2, 2, 2);

    auto B11 = subview(B, 0, 0, 2, 2);
    auto B12 = subview(B, 0, 2, 2, 2);
    auto B21 = subview(B, 2, 0, 2, 2);
    auto B22 = subview(B, 2, 2, 2, 2);

    auto C11 = subview(C, 0, 0, 2, 2);
    auto C12 = subview(C, 0, 2, 2, 2);
    auto C21 = subview(C, 2, 0, 2, 2);
    auto C22 = subview(C, 2, 2, 2, 2);

    // Временные матрицы
    Matrix<T> S1(2,2), S2(2,2), S3(2,2), S4(2,2), S5(2,2), S6(2,2), S7(2,2);
    Matrix<T> P1(2,2), P2(2,2), P3(2,2), P4(2,2), P5(2,2), P6(2,2), P7(2,2);
    Matrix<T> T1(2,2), T2(2,2);

    // 7 произведений Strassen
    sub_2x2(B12, B22, view(S1), cnt);
    mul_naive_2x2(A11, view(S1), view(P1), cnt);

    add_2x2(A11, A12, view(S2), cnt);
    mul_naive_2x2(view(S2), B22, view(P2), cnt);

    add_2x2(A21, A22, view(S3), cnt);
    mul_naive_2x2(view(S3), B11, view(P3), cnt);

    sub_2x2(B21, B11, view(S4), cnt);
    mul_naive_2x2(A22, view(S4), view(P4), cnt);

    add_2x2(A11, A22, view(S5), cnt);
    add_2x2(B11, B22, view(S6), cnt);
    mul_naive_2x2(view(S5), view(S6), view(P5), cnt);

    sub_2x2(A12, A22, view(S7), cnt);
    add_2x2(B21, B22, view(T1), cnt);
    mul_naive_2x2(view(S7), view(T1), view(P6), cnt);

    sub_2x2(A11, A21, view(T2), cnt);
    add_2x2(B11, B12, view(S1), cnt);
    mul_naive_2x2(view(T2), view(S1), view(P7), cnt);

    // Собираем результат
    add_2x2(view(P5), view(P4), view(T1), cnt);
    sub_2x2(view(T1), view(P2), view(T2), cnt);
    add_2x2(view(T2), view(P6), C11, cnt);

    add_2x2(view(P1), view(P2), C12, cnt);
    add_2x2(view(P3), view(P4), C21, cnt);

    add_2x2(view(P5), view(P1), view(T1), cnt);
    sub_2x2(view(T1), view(P3), view(T2), cnt);
    sub_2x2(view(T2), view(P7), C22, cnt);
}

#endif // ALG_STRASSEN_4X4_H
