//
// AlphaEvolve 4x4 - версия с комплексными числами из Python
//

#ifndef ALG_ALPHA_EVOLVE_4X4_COMPLEX_H
#define ALG_ALPHA_EVOLVE_4X4_COMPLEX_H

#include "structures.h"
#include <complex>

template <class T>
void alphaevolve_4x4_complex(const Matrix<T>& A,
                             const Matrix<T>& B,
                             Matrix<T>& C,
                             OpCounter* cnt = nullptr)
{
    C.resize(4,4);

    using Complex = std::complex<double>;

    auto mul_ = [&](const Complex& x, const Complex& y) -> Complex {
        if (cnt) cnt->mul++;
        return x * y;
    };
    auto add_ = [&](const Complex& x, const Complex& y) -> Complex {
        if (cnt) cnt->add++;
        return x + y;
    };

    // Константы
    const Complex half(0.5, 0.0);
    const Complex half_j(0.0, 0.5);
    const Complex half_p_half_j(0.5, 0.5);
    const Complex half_m_half_j(0.5, -0.5);
    const Complex neg_half(-0.5, 0.0);
    const Complex neg_half_j(0.0, -0.5);

    // Кэшируем элементы матриц как комплексные числа
    const Complex A00(A(0,0)), A01(A(0,1)), A02(A(0,2)), A03(A(0,3));
    const Complex A10(A(1,0)), A11(A(1,1)), A12(A(1,2)), A13(A(1,3));
    const Complex A20(A(2,0)), A21(A(2,1)), A22(A(2,2)), A23(A(2,3));
    const Complex A30(A(3,0)), A31(A(3,1)), A32(A(3,2)), A33(A(3,3));

    const Complex B00(B(0,0)), B01(B(0,1)), B02(B(0,2)), B03(B(0,3));
    const Complex B10(B(1,0)), B11(B(1,1)), B12(B(1,2)), B13(B(1,3));
    const Complex B20(B(2,0)), B21(B(2,1)), B22(B(2,2)), B23(B(2,3));
    const Complex B30(B(3,0)), B31(B(3,1)), B32(B(3,2)), B33(B(3,3));

    // Линейные комбинации элементов A
    Complex a[48];

    a[0] = half_p_half_j*A00 + half_p_half_j*A01 + half_m_half_j*A10 + half_m_half_j*A11 + half_m_half_j*A20 + half_m_half_j*A21 + half_m_half_j*A30 + half_m_half_j*A31;
    a[1] = half_p_half_j*A00 + (neg_half+half_j)*A03 + half_p_half_j*A10 + (neg_half+half_j)*A13 + (neg_half+neg_half_j)*A20 + half_m_half_j*A23 + half_m_half_j*A30 + half_p_half_j*A33;
    a[2] = neg_half*A01 + half*A02 + neg_half_j*A11 + half_j*A12 + half_j*A21 + neg_half_j*A22 + neg_half_j*A31 + half_j*A32;
    a[3] = neg_half_j*A00 + neg_half*A01 + half*A02 + neg_half*A03 + half_j*A10 + neg_half*A11 + half*A12 + half*A13 + neg_half_j*A20 + neg_half*A21 + half*A22 + neg_half*A23 + neg_half*A30 + neg_half_j*A31 + half_j*A32 + half_j*A33;
    a[4] = half_p_half_j*A00 + (neg_half+neg_half_j)*A01 + (neg_half+half_j)*A10 + half_m_half_j*A11 + (neg_half+half_j)*A20 + half_m_half_j*A21 + half_m_half_j*A30 + (neg_half+half_j)*A31;
    a[5] = half_m_half_j*A02 + (neg_half+neg_half_j)*A03 + half_m_half_j*A12 + (neg_half+neg_half_j)*A13 + (neg_half+half_j)*A22 + half_p_half_j*A23 + (neg_half+neg_half_j)*A32 + (neg_half+half_j)*A33;
    a[6] = half_j*A00 + half*A03 + neg_half*A10 + half_j*A13 + half*A20 + neg_half_j*A23 + neg_half*A30 + half_j*A33;
    a[7] = half_p_half_j*A00 + (neg_half+neg_half_j)*A01 + (neg_half+neg_half_j)*A10 + half_p_half_j*A11 + (neg_half+neg_half_j)*A20 + half_p_half_j*A21 + (neg_half+half_j)*A30 + half_m_half_j*A31;
    a[8] = neg_half_j*A00 + neg_half_j*A01 + neg_half*A02 + neg_half_j*A03 + half*A10 + half*A11 + neg_half_j*A12 + half*A13 + neg_half*A20 + neg_half*A21 + neg_half_j*A22 + half*A23 + half*A30 + half*A31 + half_j*A32 + neg_half*A33;
    a[9] = (neg_half+half_j)*A00 + (neg_half+neg_half_j)*A03 + half_p_half_j*A10 + (neg_half+half_j)*A13 + (neg_half+neg_half_j)*A20 + half_m_half_j*A23 + (neg_half+neg_half_j)*A30 + half_m_half_j*A33;
    a[10] = (neg_half+half_j)*A00 + half_m_half_j*A01 + (neg_half+half_j)*A10 + half_m_half_j*A11 + half_m_half_j*A20 + (neg_half+half_j)*A21 + half_p_half_j*A30 + (neg_half+neg_half_j)*A31;
    a[11] = half*A00 + half*A01 + neg_half_j*A02 + neg_half*A03 + neg_half*A10 + neg_half*A11 + half_j*A12 + half*A13 + half*A20 + half*A21 + half_j*A22 + half*A23 + neg_half_j*A30 + neg_half_j*A31 + half*A32 + neg_half_j*A33;
    a[12] = half_p_half_j*A01 + (neg_half+neg_half_j)*A02 + (neg_half+half_j)*A11 + half_m_half_j*A12 + (neg_half+half_j)*A21 + half_m_half_j*A22 + half_m_half_j*A31 + (neg_half+half_j)*A32;
    a[13] = half_m_half_j*A01 + (neg_half+half_j)*A02 + half_m_half_j*A11 + (neg_half+half_j)*A12 + half_m_half_j*A21 + (neg_half+half_j)*A22 + half_p_half_j*A31 + (neg_half+neg_half_j)*A32;
    a[14] = half_j*A00 + neg_half*A01 + half*A02 + neg_half*A03 + half*A10 + neg_half_j*A11 + half_j*A12 + half_j*A13 + half*A20 + half_j*A21 + neg_half_j*A22 + half_j*A23 + half*A30 + neg_half_j*A31 + half_j*A32 + half_j*A33;
    a[15] = (neg_half+half_j)*A02 + half_p_half_j*A03 + half_m_half_j*A12 + (neg_half+neg_half_j)*A13 + half_m_half_j*A22 + (neg_half+neg_half_j)*A23 + (neg_half+neg_half_j)*A32 + (neg_half+half_j)*A33;
    a[16] = neg_half*A00 + half_j*A01 + half_j*A02 + neg_half_j*A03 + neg_half*A10 + neg_half_j*A11 + neg_half_j*A12 + neg_half_j*A13 + neg_half*A20 + half_j*A21 + half_j*A22 + neg_half_j*A23 + neg_half_j*A30 + half*A31 + half*A32 + half*A33;
    a[17] = half_p_half_j*A00 + half_p_half_j*A01 + half_p_half_j*A10 + half_p_half_j*A11 + half_p_half_j*A20 + half_p_half_j*A21 + (neg_half+half_j)*A30 + (neg_half+half_j)*A31;
    a[18] = half_j*A00 + half_j*A01 + neg_half*A02 + half_j*A03 + half_j*A10 + half_j*A11 + neg_half*A12 + half_j*A13 + half_j*A20 + half_j*A21 + half*A22 + neg_half_j*A23 + neg_half*A30 + neg_half*A31 + half_j*A32 + half*A33;
    a[19] = half_m_half_j*A02 + half_p_half_j*A03 + half_m_half_j*A12 + half_p_half_j*A13 + half_m_half_j*A22 + half_p_half_j*A23 + half_p_half_j*A32 + (neg_half+half_j)*A33;
    a[20] = half_p_half_j*A01 + (neg_half+neg_half_j)*A02 + half_p_half_j*A11 + (neg_half+neg_half_j)*A12 + (neg_half+neg_half_j)*A21 + half_p_half_j*A22 + half_m_half_j*A31 + (neg_half+half_j)*A32;
    a[21] = half_j*A00 + neg_half_j*A01 + neg_half*A02 + neg_half_j*A03 + neg_half_j*A10 + half_j*A11 + half*A12 + half_j*A13 + neg_half_j*A20 + half_j*A21 + neg_half*A22 + neg_half_j*A23 + neg_half*A30 + half*A31 + half_j*A32 + neg_half*A33;
    a[22] = (neg_half+neg_half_j)*A00 + (neg_half+half_j)*A03 + half_m_half_j*A10 + (neg_half+neg_half_j)*A13 + half_m_half_j*A20 + (neg_half+neg_half_j)*A23 + (neg_half+half_j)*A30 + half_p_half_j*A33;
    a[23] = (neg_half+neg_half_j)*A02 + half_m_half_j*A03 + half_m_half_j*A12 + half_p_half_j*A13 + half_m_half_j*A22 + half_p_half_j*A23 + (neg_half+half_j)*A32 + (neg_half+neg_half_j)*A33;
    a[24] = neg_half*A00 + half*A01 + neg_half_j*A02 + neg_half*A03 + neg_half_j*A10 + half_j*A11 + half*A12 + neg_half_j*A13 + neg_half_j*A20 + half_j*A21 + neg_half*A22 + half_j*A23 + half_j*A30 + neg_half_j*A31 + half*A32 + neg_half_j*A33;
    a[25] = half_m_half_j*A02 + half_p_half_j*A03 + (neg_half+neg_half_j)*A12 + half_m_half_j*A13 + half_p_half_j*A22 + (neg_half+half_j)*A23 + half_p_half_j*A32 + (neg_half+half_j)*A33;
    a[26] = half_p_half_j*A01 + half_p_half_j*A02 + (neg_half+neg_half_j)*A11 + (neg_half+neg_half_j)*A12 + half_p_half_j*A21 + half_p_half_j*A22 + half_m_half_j*A31 + half_m_half_j*A32;
    a[27] = neg_half_j*A00 + neg_half_j*A01 + half*A02 + half_j*A03 + neg_half*A10 + neg_half*A11 + neg_half_j*A12 + half*A13 + neg_half*A20 + neg_half*A21 + half_j*A22 + neg_half*A23 + neg_half*A30 + neg_half*A31 + half_j*A32 + neg_half*A33;
    a[28] = (neg_half+half_j)*A00 + (neg_half+half_j)*A01 + (neg_half+neg_half_j)*A10 + (neg_half+neg_half_j)*A11 + half_p_half_j*A20 + half_p_half_j*A21 + (neg_half+neg_half_j)*A30 + (neg_half+neg_half_j)*A31;
    a[29] = half_p_half_j*A00 + half_m_half_j*A03 + (neg_half+neg_half_j)*A10 + (neg_half+half_j)*A13 + half_p_half_j*A20 + half_m_half_j*A23 + half_m_half_j*A30 + (neg_half+neg_half_j)*A33;
    a[30] = half_p_half_j*A01 + half_p_half_j*A02 + (neg_half+neg_half_j)*A11 + (neg_half+neg_half_j)*A12 + (neg_half+neg_half_j)*A21 + (neg_half+neg_half_j)*A22 + (neg_half+half_j)*A31 + (neg_half+half_j)*A32;
    a[31] = half*A00 + neg_half*A01 + neg_half_j*A02 + half*A03 + half*A10 + neg_half*A11 + neg_half_j*A12 + half*A13 + neg_half*A20 + half*A21 + neg_half_j*A22 + half*A23 + neg_half_j*A30 + half_j*A31 + half*A32 + half_j*A33;
    a[32] = half_p_half_j*A02 + half_m_half_j*A03 + (neg_half+half_j)*A12 + half_p_half_j*A13 + half_m_half_j*A22 + (neg_half+neg_half_j)*A23 + (neg_half+half_j)*A32 + half_p_half_j*A33;
    a[33] = half*A00 + half_j*A01 + neg_half_j*A02 + neg_half_j*A03 + neg_half*A10 + half_j*A11 + neg_half_j*A12 + half_j*A13 + neg_half*A20 + neg_half_j*A21 + half_j*A22 + half_j*A23 + half_j*A30 + half*A31 + neg_half*A32 + half*A33;
    a[34] = neg_half_j*A00 + half_j*A01 + neg_half*A02 + half_j*A03 + neg_half*A10 + half*A11 + half_j*A12 + half*A13 + half*A20 + neg_half*A21 + half_j*A22 + half*A23 + half*A30 + neg_half*A31 + half_j*A32 + half*A33;
    a[35] = half_m_half_j*A02 + half_p_half_j*A03 + (neg_half+half_j)*A12 + (neg_half+neg_half_j)*A13 + half_m_half_j*A22 + half_p_half_j*A23 + (neg_half+neg_half_j)*A32 + half_m_half_j*A33;
    a[36] = (neg_half+neg_half_j)*A01 + (neg_half+neg_half_j)*A02 + (neg_half+half_j)*A11 + (neg_half+half_j)*A12 + half_m_half_j*A21 + half_m_half_j*A22 + half_m_half_j*A31 + half_m_half_j*A32;
    a[37] = half*A00 + neg_half_j*A01 + neg_half_j*A02 + neg_half_j*A03 + half_j*A10 + neg_half*A11 + neg_half*A12 + half*A13 + half_j*A20 + half*A21 + half*A22 + half*A23 + neg_half_j*A30 + half*A31 + half*A32 + neg_half*A33;
    a[38] = half_m_half_j*A01 + half_m_half_j*A02 + (neg_half+neg_half_j)*A11 + (neg_half+neg_half_j)*A12 + (neg_half+neg_half_j)*A21 + (neg_half+neg_half_j)*A22 + (neg_half+neg_half_j)*A31 + (neg_half+neg_half_j)*A32;
    a[39] = neg_half*A00 + neg_half_j*A01 + neg_half_j*A02 + neg_half_j*A03 + neg_half*A10 + half_j*A11 + half_j*A12 + neg_half_j*A13 + half*A20 + half_j*A21 + half_j*A22 + half_j*A23 + half_j*A30 + half*A31 + half*A32 + neg_half*A33;
    a[40] = (neg_half+neg_half_j)*A00 + (neg_half+neg_half_j)*A01 + half_p_half_j*A10 + half_p_half_j*A11 + (neg_half+neg_half_j)*A20 + (neg_half+neg_half_j)*A21 + (neg_half+half_j)*A30 + (neg_half+half_j)*A31;
    a[41] = half_m_half_j*A00 + (neg_half+neg_half_j)*A03 + (neg_half+half_j)*A10 + half_p_half_j*A13 + (neg_half+half_j)*A20 + half_p_half_j*A23 + half_p_half_j*A30 + half_m_half_j*A33;
    a[42] = half_p_half_j*A00 + (neg_half+half_j)*A03 + half_m_half_j*A10 + half_p_half_j*A13 + half_m_half_j*A20 + half_p_half_j*A23 + half_m_half_j*A30 + half_p_half_j*A33;
    a[43] = half_j*A00 + half*A01 + neg_half*A02 + neg_half*A03 + half*A10 + half_j*A11 + neg_half_j*A12 + half_j*A13 + neg_half*A20 + half_j*A21 + neg_half_j*A22 + neg_half_j*A23 + neg_half*A30 + neg_half_j*A31 + half_j*A32 + neg_half_j*A33;
    a[44] = half_m_half_j*A02 + (neg_half+neg_half_j)*A03 + (neg_half+neg_half_j)*A12 + (neg_half+half_j)*A13 + (neg_half+neg_half_j)*A22 + (neg_half+half_j)*A23 + (neg_half+neg_half_j)*A32 + (neg_half+half_j)*A33;
    a[45] = (neg_half+half_j)*A00 + half_m_half_j*A01 + half_p_half_j*A10 + (neg_half+neg_half_j)*A11 + (neg_half+neg_half_j)*A20 + half_p_half_j*A21 + (neg_half+neg_half_j)*A30 + half_p_half_j*A31;
    a[46] = half_m_half_j*A00 + half_p_half_j*A03 + half_m_half_j*A10 + half_p_half_j*A13 + half_m_half_j*A20 + half_p_half_j*A23 + half_p_half_j*A30 + (neg_half+half_j)*A33;
    a[47] = half*A00 + half_j*A01 + half_j*A02 + neg_half_j*A03 + half_j*A10 + half*A11 + half*A12 + half*A13 + neg_half_j*A20 + half*A21 + half*A22 + neg_half*A23 + half_j*A30 + half*A31 + half*A32 + half*A33;

    // Линейные комбинации элементов B
    Complex b[48];

    b[0] = neg_half*B00 + neg_half*B10 + half*B20 + neg_half_j*B30;
    b[1] = half_j*B01 + half_j*B03 + half_j*B11 + half_j*B13 + half_j*B21 + half_j*B23 + half*B31 + half*B33;
    b[2] = half_p_half_j*B01 + (neg_half+neg_half_j)*B11 + half_p_half_j*B21 + half_m_half_j*B31;
    b[3] = neg_half_j*B00 + half_j*B02 + neg_half_j*B11 + neg_half_j*B12 + half_j*B21 + half_j*B22 + half*B30 + neg_half*B32;
    b[4] = neg_half*B00 + half*B02 + half*B03 + half*B10 + neg_half*B12 + neg_half*B13 + half*B20 + neg_half*B22 + neg_half*B23 + half_j*B30 + neg_half_j*B32 + neg_half_j*B33;
    b[5] = half*B01 + half*B03 + half*B11 + half*B13 + half*B21 + half*B23 + half_j*B31 + half_j*B33;
    b[6] = (neg_half+neg_half_j)*B01 + half_p_half_j*B11 + half_p_half_j*B21 + half_m_half_j*B31;
    b[7] = neg_half*B00 + half*B03 + half*B10 + neg_half*B13 + neg_half*B20 + half*B23 + half_j*B30 + neg_half_j*B33;
    b[8] = half*B00 + neg_half*B02 + neg_half*B03 + half*B10 + neg_half*B12 + neg_half*B13 + half*B21 + neg_half_j*B31;
    b[9] = half_j*B01 + half_j*B02 + half_j*B03 + half_j*B11 + half_j*B12 + half_j*B13 + neg_half_j*B21 + neg_half_j*B22 + neg_half_j*B23 + half*B31 + half*B32 + half*B33;
    b[10] = half_j*B01 + half_j*B03 + neg_half_j*B11 + neg_half_j*B13 + neg_half_j*B21 + neg_half_j*B23 + neg_half*B31 + neg_half*B33;
    b[11] = neg_half_j*B00 + half_j*B03 + neg_half_j*B10 + half_j*B13 + half_j*B21 + half_j*B22 + neg_half*B31 + neg_half*B32;
    b[12] = neg_half*B00 + half*B02 + half*B03 + neg_half*B10 + half*B12 + half*B13 + half*B20 + neg_half*B22 + neg_half*B23 + half_j*B30 + neg_half_j*B32 + neg_half_j*B33;
    b[13] = half_j*B00 + neg_half_j*B02 + neg_half_j*B10 + half_j*B12 + half_j*B20 + neg_half_j*B22 + neg_half*B30 + half*B32;
    b[14] = neg_half*B01 + neg_half*B10 + half*B20 + half_j*B31;
    b[15] = half_j*B00 + neg_half_j*B03 + half_j*B10 + neg_half_j*B13 + neg_half_j*B20 + half_j*B23 + half*B30 + neg_half*B33;
    b[16] = half*B01 + half*B02 + half*B10 + neg_half*B12 + half*B20 + neg_half*B22 + neg_half_j*B31 + neg_half_j*B32;
    b[17] = neg_half_j*B00 + half_j*B02 + neg_half_j*B10 + half_j*B12 + neg_half_j*B20 + half_j*B22 + half*B30 + neg_half*B32;
    b[18] = neg_half_j*B01 + neg_half_j*B03 + neg_half_j*B11 + neg_half_j*B13 + neg_half_j*B20 + half_j*B22 + half*B30 + neg_half*B32;
    b[19] = neg_half_j*B00 + half_j*B02 + half_j*B10 + neg_half_j*B12 + half_j*B20 + neg_half_j*B22 + half*B30 + neg_half*B32;
    b[20] = neg_half_j*B01 + neg_half_j*B03 + neg_half_j*B11 + neg_half_j*B13 + half_j*B21 + half_j*B23 + half*B31 + half*B33;
    b[21] = neg_half*B01 + neg_half*B02 + half*B11 + half*B12 + neg_half*B20 + half*B23 + half_j*B30 + neg_half_j*B33;
    b[22] = neg_half_j*B00 + half_j*B02 + half_j*B03 + neg_half_j*B10 + half_j*B12 + half_j*B13 + neg_half_j*B20 + half_j*B22 + half_j*B23 + half*B30 + neg_half*B32 + neg_half*B33;
    b[23] = neg_half*B00 + half*B02 + half*B03 + neg_half*B10 + half*B12 + half*B13 + neg_half*B20 + half*B22 + half*B23 + half_j*B30 + neg_half_j*B32 + neg_half_j*B33;
    b[24] = half_j*B01 + neg_half_j*B11 + neg_half_j*B20 + half_j*B22 + half_j*B23 + half*B30 + neg_half*B32 + neg_half*B33;
    b[25] = half_j*B01 + half_j*B02 + half_j*B03 + half_j*B11 + half_j*B12 + half_j*B13 + neg_half_j*B21 + neg_half_j*B22 + neg_half_j*B23 + neg_half*B31 + neg_half*B32 + neg_half*B33;
    b[26] = half*B01 + half*B02 + neg_half*B11 + neg_half*B12 + neg_half*B21 + neg_half*B22 + neg_half_j*B31 + neg_half_j*B32;
    b[27] = half_j*B01 + half_j*B02 + half_j*B03 + half_j*B11 + half_j*B12 + half_j*B13 + neg_half_j*B20 + neg_half*B30;
    b[28] = half*B01 + half*B11 + half*B21 + neg_half_j*B31;
    b[29] = half_j*B01 + half_j*B02 + neg_half_j*B11 + neg_half_j*B12 + half_j*B21 + half_j*B22 + neg_half*B31 + neg_half*B32;
    b[30] = neg_half*B00 + half*B03 + neg_half*B10 + half*B13 + neg_half*B20 + half*B23 + half_j*B30 + neg_half_j*B33;
    b[31] = half*B00 + neg_half*B02 + neg_half*B10 + half*B12 + neg_half*B21 + neg_half*B23 + half_j*B31 + half_j*B33;
    b[32] = half_j*B01 + neg_half_j*B11 + neg_half_j*B21 + half*B31;
    b[33] = neg_half*B01 + neg_half*B03 + half*B10 + neg_half*B13 + neg_half*B20 + half*B23 + neg_half_j*B31 + neg_half_j*B33;
    b[34] = half_j*B00 + neg_half_j*B10 + half_j*B21 + half_j*B22 + half_j*B23 + neg_half*B31 + neg_half*B32 + neg_half*B33;
    b[35] = neg_half_j*B01 + neg_half_j*B02 + half_j*B11 + half_j*B12 + neg_half_j*B21 + neg_half_j*B22 + neg_half*B31 + neg_half*B32;
    b[36] = neg_half*B01 + neg_half*B02 + neg_half*B03 + neg_half*B11 + neg_half*B12 + neg_half*B13 + neg_half*B21 + neg_half*B22 + neg_half*B23 + neg_half_j*B31 + neg_half_j*B32 + neg_half_j*B33;
    b[37] = half_j*B01 + half_j*B02 + half_j*B03 + neg_half_j*B10 + half_j*B12 + half_j*B13 + neg_half_j*B20 + half_j*B22 + half_j*B23 + neg_half*B31 + neg_half*B32 + neg_half*B33;
    b[38] = half_j*B00 + neg_half_j*B10 + neg_half_j*B20 + neg_half*B30;
    b[39] = neg_half_j*B00 + half_j*B03 + half_j*B11 + half_j*B13 + half_j*B21 + half_j*B23 + neg_half*B30 + half*B33;
    b[40] = half_j*B01 + half_j*B02 + half_j*B11 + half_j*B12 + neg_half_j*B21 + neg_half_j*B22 + half*B31 + half*B32;
    b[41] = half*B00 + neg_half*B03 + half*B10 + neg_half*B13 + neg_half*B20 + half*B23 + half_j*B30 + neg_half_j*B33;
    b[42] = half_j*B00 + neg_half_j*B10 + half_j*B20 + half*B30;
    b[43] = half*B00 + neg_half*B02 + neg_half*B03 + neg_half*B11 + neg_half*B12 + neg_half*B13 + half*B21 + half*B22 + half*B23 + neg_half_j*B30 + half_j*B32 + half_j*B33;
    b[44] = neg_half_j*B00 + half_j*B10 + neg_half_j*B20 + half*B30;
    b[45] = neg_half_j*B01 + neg_half_j*B02 + neg_half_j*B03 + half_j*B11 + half_j*B12 + half_j*B13 + neg_half_j*B21 + neg_half_j*B22 + neg_half_j*B23 + half*B31 + half*B32 + half*B33;
    b[46] = neg_half*B00 + half*B02 + half*B10 + neg_half*B12 + half*B20 + neg_half*B22 + half_j*B30 + neg_half_j*B32;
    b[47] = half*B00 + half*B11 + half*B21 + half_j*B30;

    // 48 умножений
    Complex m[48];
    for (int i = 0; i < 48; i++) {
        m[i] = mul_(a[i], b[i]);
    }

    // Восстановление C
    Complex C_complex[4][4];

    C_complex[0][0] = half_j*m[0] + neg_half_j*m[1] + neg_half*m[5] + half*m[8] + half_j*m[9] + (neg_half+half_j)*m[11] + half*m[14] + neg_half_j*m[15] + (neg_half+neg_half_j)*m[16] + half_j*m[17] + (neg_half+neg_half_j)*m[18] + neg_half_j*m[24] + half_j*m[26] + half_j*m[27] + half*m[28] + half_j*m[30] + neg_half_j*m[32] + half*m[34] + half*m[36] + neg_half_j*m[37] + neg_half*m[38] + (half+neg_half_j)*m[39] + neg_half_j*m[40] + neg_half*m[42] + neg_half*m[43] + neg_half*m[44] + neg_half_j*m[46] + half*m[47];

    C_complex[0][1] = neg_half_j*m[0] + half*m[2] + (neg_half+neg_half_j)*m[3] + half*m[5] + half*m[6] + neg_half*m[8] + (half+neg_half_j)*m[11] + neg_half*m[12] + half_j*m[13] + half_j*m[14] + half_j*m[15] + neg_half_j*m[17] + (half+half_j)*m[18] + half*m[20] + neg_half*m[22] + half_j*m[24] + neg_half_j*m[27] + neg_half*m[28] + neg_half_j*m[29] + half_j*m[32] + (neg_half+neg_half_j)*m[33] + neg_half*m[34] + neg_half*m[37] + half_j*m[40] + half_j*m[41] + neg_half_j*m[43] + half*m[44] + neg_half_j*m[47];

    C_complex[0][2] = neg_half*m[2] + half*m[3] + neg_half*m[5] + neg_half_j*m[8] + half_j*m[11] + half*m[12] + neg_half_j*m[13] + neg_half_j*m[14] + neg_half_j*m[15] + neg_half*m[16] + neg_half*m[18] + half_j*m[19] + neg_half*m[20] + half_j*m[21] + neg_half*m[23] + neg_half_j*m[24] + neg_half*m[25] + half_j*m[26] + half*m[27] + half_j*m[30] + neg_half*m[31] + neg_half_j*m[32] + half*m[33] + half*m[34] + half_j*m[35] + half*m[36] + neg_half_j*m[37] + neg_half*m[38] + neg_half_j*m[39] + half_j*m[43] + neg_half*m[44] + half*m[47];

    C_complex[0][3] = half_j*m[0] + neg_half_j*m[1] + half_j*m[3] + neg_half_j*m[4] + neg_half*m[6] + half*m[7] + half*m[8] + half_j*m[9] + neg_half*m[10] + neg_half*m[11] + half*m[14] + neg_half_j*m[16] + half_j*m[17] + neg_half_j*m[18] + neg_half*m[21] + half*m[22] + half*m[24] + half_j*m[27] + half*m[28] + half_j*m[29] + neg_half_j*m[31] + half_j*m[33] + half_j*m[34] + half*m[37] + half*m[39] + neg_half_j*m[40] + neg_half_j*m[41] + neg_half*m[42] + neg_half*m[43] + neg_half_j*m[45] + neg_half_j*m[46] + half_j*m[47];

    C_complex[1][0] = neg_half*m[0] + neg_half*m[1] + neg_half*m[5] + neg_half_j*m[8] + neg_half_j*m[9] + (half+neg_half_j)*m[11] + neg_half_j*m[14] + half_j*m[15] + (neg_half+half_j)*m[16] + half_j*m[17] + (neg_half+neg_half_j)*m[18] + neg_half*m[24] + half*m[26] + neg_half*m[27] + neg_half_j*m[28] + half*m[30] + neg_half*m[32] + half_j*m[34] + half*m[36] + neg_half*m[37] + neg_half*m[38] + (neg_half+neg_half_j)*m[39] + half_j*m[40] + half*m[42] + half_j*m[43] + neg_half_j*m[44] + neg_half*m[46] + neg_half_j*m[47];

    C_complex[1][1] = half*m[0] + neg_half*m[2] + (half+neg_half_j)*m[3] + half*m[5] + half*m[6] + half_j*m[8] + (neg_half+half_j)*m[11] + half*m[12] + neg_half*m[13] + neg_half*m[14] + neg_half_j*m[15] + neg_half_j*m[17] + (half+half_j)*m[18] + half_j*m[20] + neg_half*m[22] + half*m[24] + half*m[27] + half_j*m[28] + half*m[29] + half*m[32] + (half+neg_half_j)*m[33] + neg_half_j*m[34] + neg_half_j*m[37] + neg_half_j*m[40] + neg_half*m[41] + half*m[43] + half_j*m[44] + half*m[47];

    C_complex[1][2] = half*m[2] + neg_half*m[3] + neg_half*m[5] + neg_half*m[8] + neg_half_j*m[11] + neg_half*m[12] + half*m[13] + half*m[14] + half_j*m[15] + neg_half*m[16] + neg_half*m[18] + half_j*m[19] + neg_half_j*m[20] + neg_half_j*m[21] + half_j*m[23] + neg_half*m[24] + neg_half_j*m[25] + half*m[26] + half_j*m[27] + half*m[30] + neg_half*m[31] + neg_half*m[32] + neg_half*m[33] + half_j*m[34] + neg_half_j*m[35] + half*m[36] + neg_half*m[37] + neg_half*m[38] + neg_half_j*m[39] + neg_half*m[43] + neg_half_j*m[44] + neg_half_j*m[47];

    C_complex[1][3] = neg_half*m[0] + neg_half*m[1] + half_j*m[3] + neg_half*m[4] + neg_half*m[6] + neg_half*m[7] + neg_half_j*m[8] + neg_half_j*m[9] + neg_half*m[10] + half*m[11] + neg_half_j*m[14] + half_j*m[16] + half_j*m[17] + neg_half_j*m[18] + half*m[21] + half*m[22] + neg_half_j*m[24] + neg_half*m[27] + neg_half_j*m[28] + neg_half*m[29] + neg_half_j*m[31] + half_j*m[33] + neg_half*m[34] + half_j*m[37] + neg_half*m[39] + half_j*m[40] + half*m[41] + half*m[42] + half_j*m[43] + half*m[45] + neg_half*m[46] + neg_half*m[47];

    C_complex[2][0] = neg_half_j*m[0] + half_j*m[1] + half_j*m[5] + neg_half_j*m[8] + half*m[9] + (half+half_j)*m[11] + half_j*m[14] + neg_half*m[15] + (neg_half+neg_half_j)*m[16] + half*m[17] + (neg_half+half_j)*m[18] + neg_half*m[24] + half_j*m[26] + half*m[27] + neg_half*m[28] + neg_half_j*m[30] + neg_half_j*m[32] + neg_half_j*m[34] + neg_half_j*m[36] + neg_half*m[37] + neg_half_j*m[38] + (neg_half+half_j)*m[39] + neg_half*m[40] + neg_half_j*m[42] + half_j*m[43] + neg_half*m[44] + neg_half_j*m[46] + half_j*m[47];

    C_complex[2][1] = half_j*m[0] + half_j*m[2] + (neg_half+neg_half_j)*m[3] + neg_half_j*m[5] + half_j*m[6] + half_j*m[8] + (neg_half+neg_half_j)*m[11] + half_j*m[12] + half_j*m[13] + neg_half*m[14] + half*m[15] + neg_half*m[17] + (half+neg_half_j)*m[18] + neg_half*m[20] + half_j*m[22] + half*m[24] + neg_half*m[27] + half*m[28] + neg_half_j*m[29] + half_j*m[32] + (half+half_j)*m[33] + half_j*m[34] + half_j*m[37] + half*m[40] + neg_half_j*m[41] + neg_half*m[43] + half*m[44] + half*m[47];

    C_complex[2][2] = neg_half_j*m[2] + half*m[3] + half_j*m[5] + half*m[8] + half_j*m[11] + neg_half_j*m[12] + neg_half_j*m[13] + half*m[14] + neg_half*m[15] + neg_half*m[16] + neg_half*m[18] + neg_half*m[19] + half*m[20] + neg_half_j*m[21] + half*m[23] + neg_half*m[24] + half*m[25] + half_j*m[26] + half_j*m[27] + neg_half_j*m[30] + half*m[31] + neg_half_j*m[32] + neg_half*m[33] + neg_half_j*m[34] + neg_half*m[35] + neg_half_j*m[36] + neg_half*m[37] + neg_half_j*m[38] + half_j*m[39] + half*m[43] + neg_half*m[44] + half_j*m[47];

    C_complex[2][3] = neg_half_j*m[0] + half_j*m[1] + half_j*m[3] + neg_half_j*m[4] + neg_half_j*m[6] + half_j*m[7] + neg_half_j*m[8] + half*m[9] + neg_half_j*m[10] + half*m[11] + half_j*m[14] + neg_half_j*m[16] + half*m[17] + half_j*m[18] + neg_half*m[21] + neg_half_j*m[22] + half_j*m[24] + half*m[27] + neg_half*m[28] + half_j*m[29] + neg_half_j*m[31] + neg_half_j*m[33] + neg_half*m[34] + neg_half_j*m[37] + neg_half*m[39] + neg_half*m[40] + half_j*m[41] + neg_half_j*m[42] + half_j*m[43] + neg_half_j*m[45] + neg_half_j*m[46] + neg_half*m[47];

    C_complex[3][0] = neg_half_j*m[0] + neg_half_j*m[1] + half*m[5] + half_j*m[8] + half_j*m[9] + (neg_half+half_j)*m[11] + neg_half_j*m[14] + neg_half_j*m[15] + (half+half_j)*m[16] + neg_half_j*m[17] + (half+half_j)*m[18] + half*m[24] + neg_half_j*m[26] + half*m[27] + half*m[28] + half_j*m[30] + half_j*m[32] + neg_half_j*m[34] + neg_half*m[36] + half*m[37] + neg_half*m[38] + (half+neg_half_j)*m[39] + neg_half_j*m[40] + half*m[42] + neg_half_j*m[43] + neg_half*m[44] + half_j*m[46] + neg_half_j*m[47];

    C_complex[3][1] = half_j*m[0] + neg_half*m[2] + (neg_half+neg_half_j)*m[3] + neg_half*m[5] + half*m[6] + neg_half_j*m[8] + (half+neg_half_j)*m[11] + neg_half*m[12] + half_j*m[13] + neg_half*m[14] + half_j*m[15] + half_j*m[17] + (neg_half+neg_half_j)*m[18] + neg_half*m[20] + half*m[22] + neg_half*m[24] + neg_half*m[27] + neg_half*m[28] + neg_half_j*m[29] + neg_half_j*m[32] + (half+half_j)*m[33] + half_j*m[34] + half_j*m[37] + half_j*m[40] + neg_half_j*m[41] + neg_half*m[43] + half*m[44] + half*m[47];

    C_complex[3][2] = half*m[2] + half_j*m[3] + half*m[5] + neg_half*m[8] + neg_half*m[11] + half*m[12] + neg_half_j*m[13] + half*m[14] + neg_half_j*m[15] + half_j*m[16] + half_j*m[18] + half_j*m[19] + half*m[20] + half*m[21] + neg_half*m[23] + half*m[24] + half*m[25] + neg_half_j*m[26] + half_j*m[27] + half_j*m[30] + neg_half_j*m[31] + half_j*m[32] + neg_half_j*m[33] + neg_half_j*m[34] + neg_half_j*m[35] + neg_half*m[36] + half*m[37] + neg_half*m[38] + half*m[39] + half*m[43] + neg_half*m[44] + neg_half_j*m[47];

    C_complex[3][3] = neg_half_j*m[0] + neg_half_j*m[1] + half*m[3] + half_j*m[4] + neg_half*m[6] + neg_half*m[7] + half_j*m[8] + half_j*m[9] + neg_half*m[10] + half_j*m[11] + neg_half_j*m[14] + half*m[16] + neg_half_j*m[17] + half*m[18] + neg_half_j*m[21] + neg_half*m[22] + neg_half_j*m[24] + half*m[27] + half*m[28] + half_j*m[29] + neg_half*m[31] + neg_half*m[33] + neg_half*m[34] + neg_half_j*m[37] + neg_half_j*m[39] + neg_half_j*m[40] + half_j*m[41] + half*m[42] + neg_half_j*m[43] + neg_half_j*m[45] + half_j*m[46] + neg_half*m[47];

    // Извлекаем действительную часть для реальных матриц
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            C(i,j) = C_complex[i][j].real();
        }
    }
}

#endif // ALG_ALPHA_EVOLVE_4X4_COMPLEX_H
