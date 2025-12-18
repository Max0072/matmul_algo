#include "matrix_io.h"
#include "structures.h"
#include "alg_naive.h"
#include "alg_strassen.h"
#include "generators.h"
#include "rss.h"

//auto A1 = gen_random<int>(4,4,1,-5,5);
//auto A2 = gen_random<double>(4,4,2,-1,1);
//auto A3 = gen_random<std::complex<double>>(4,4,3,-1,1);

//auto S1 = gen_symmetric<double>(8,4,-1,1);
//auto Z1 = gen_almost_sparse<double>(8,8,0.95,5,-1,1);

int main() {
//    int n, m, p;
//    std::cin >> n >> m >> p;
//
//    Matrix<int> A(n, m), B(m, p), C;
//
//    read_matrix(A);
//    read_matrix(B);
//    auto A = gen_random<std::complex<double>>(4,4,3,-1,1);
//    auto B = gen_random<std::complex<double>>(4,4,4,-1,1);
//    Matrix<std::complex<double>> C;

    auto A = gen_random<int>(2,2,1,-5,5);
    auto B = gen_random<int>(2,2,2,-5,5);
    Matrix<int> C;



    OpCounter cnt;                                            // counter

    uint64_t rss0 = current_rss_bytes();
    auto t0 = std::chrono::steady_clock::now(); // time 1

//    mul_naive(A, B, C, &cnt);                              // matmul
    mul_strassen(A, B, C, /*THRESH=*/64, &cnt);

    auto t1 = std::chrono::steady_clock::now(); // time 2
    uint64_t rss1 = current_rss_bytes();

    print_matrix(C);

    std::cout << "rss_delta=" << (rss1 - rss0) << "\n";
    std::cout << "mul=" << cnt.mul << " add=" << cnt.add << "\n";

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "Time: " << ms << "\n";

    return 0;
}

//2 3 2
//1 2 3
//4 5 6
//7 8
//9 10
//11 12