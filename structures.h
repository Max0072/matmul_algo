//
// Created on 18/12/2025.
//

#ifndef UNTITLED3_STRUCTURES_H
#define UNTITLED3_STRUCTURES_H

#include <vector>
#include <iostream>
#include <cassert>
#include <type_traits>

///--------------------------
///        Matrix
///--------------------------
template <class T>
struct Matrix {
    int rows = 0, cols = 0;
    std::vector<T> a;

    Matrix() = default;
    Matrix(int r, int c) : rows(r), cols(c), a((size_t)r * c, T{}) {}

    void resize(int r, int c) {
        rows = r; cols = c;
        a.assign((size_t)r * (size_t)c, T{});
    }

    T* data() { return a.data(); }
    const T* data() const { return a.data(); }
    int stride() const { return cols; }

    T& operator()(int i, int j) {
        assert(0 <= i and i < rows and 0 <= j and j < cols);
        return a[i * cols + j];
    }

    const T& operator()(int i, int j) const {
        assert(0 <= i and i < rows and 0 <= j and j < cols);
        return a[i * cols + j];
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& A) {
        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < A.cols; ++j) os << A(i,j) << ' ';
            os << '\n';
        }
        return os;
    }
};


///--------------------------
///        MatrixView
///--------------------------
template <class T>
struct MatrixView {
    T* ptr = nullptr;
    int rows = 0, cols = 0, stride = 0;

    MatrixView() = default;
    MatrixView(T* p, int r, int c, int s)
            : ptr(p), rows(r), cols(c), stride(s) {}

    T& operator()(int i,int j) {
        assert(ptr != nullptr);
        assert(0 <= i and i < rows);
        assert(0 <= j and j < cols);
        return ptr[(size_t)i * stride + j];
    }
    const T& operator()(int i,int j) const {
        assert(ptr != nullptr);
        assert(0 <= i and i < rows);
        assert(0 <= j and j < cols);
        return ptr[(size_t)i * stride + j];
    }

    template <class U,
            typename std::enable_if<std::is_convertible<U*, T*>::value, int>::type = 0>
    MatrixView(const MatrixView<U>& other)
            : ptr(other.ptr), rows(other.rows), cols(other.cols), stride(other.stride) {}
};

template <class T>
MatrixView<T> view(Matrix<T>& M) {
    return { M.data(), M.rows, M.cols, M.stride() };
}

template <class T>
MatrixView<const T> view(const Matrix<T>& M) {
    return { M.data(), M.rows, M.cols, M.stride() };
}

template <class T>
auto subview(MatrixView<T> V, int r0,int c0,int r,int c) {
    assert(V.ptr != nullptr);
    assert(r0 >= 0 and c0 >= 0);
    assert(r >= 0 and c >= 0);
    assert(r0 + r <= V.rows);
    assert(c0 + c <= V.cols);

    return MatrixView<T>{ V.ptr + (size_t)r0 * V.stride + c0,
                          r,
                          c,
                          V.stride
    };
}


///--------------------------
///        Counter
///--------------------------
struct OpCounter { uint64_t mul=0, add=0; };


///--------------------------
///        Operations
///--------------------------
template <class T>
inline T add(const T& x, const T& y, OpCounter* c){ if(c) c->add++; return x+y; }

template <class T>
inline T sub(const T& x, const T& y, OpCounter* c){ if(c) c->add++; return x-y; }

template <class T>
inline T mul(const T& x, const T& y, OpCounter* c){ if(c) c->mul++; return x*y; }


///-------------------------
///        Finish
///------------------------


#endif //UNTITLED3_STRUCTURES_H
