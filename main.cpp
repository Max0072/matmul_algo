#include "benchmark.h"
#include "alg_naive.h"
#include "alg_strassen.h"
#include "alg_strassen_4x4.h"
#include "alg_winograd_4x4.h"
#include "alg_alpha_evolve_4x4_complex.h"
#include "alg_blocked.h"
#include <complex>

// Wrapper функции для бенчмарков
template<class T>
void wrapper_naive(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_naive(A, B, C, cnt);
}

template<class T>
void wrapper_strassen(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_strassen(A, B, C, /*THRESH=*/64, cnt);
}

template<class T>
void wrapper_winograd_4x4(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_winograd_4x4(A, B, C, cnt);
}

template<class T>
void wrapper_alphaevolve_4x4(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    alphaevolve_4x4_complex(A, B, C, cnt);
}

template<class T>
void wrapper_blocked_naive(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_blocked_naive_kernel(A, B, C, cnt);
}

template<class T>
void wrapper_blocked_winograd(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_blocked_winograd_kernel(A, B, C, cnt);
}

template<class T>
void wrapper_blocked_alphaevolve(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_blocked_alphaevolve_kernel(A, B, C, cnt);
}

template<class T>
void wrapper_strassen_4x4(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_strassen_4x4(A, B, C, cnt);
}

template<class T>
void wrapper_blocked_strassen(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, OpCounter* cnt) {
    mul_blocked_strassen_kernel(A, B, C, cnt);
}

// Функция для запуска бенчмарков для одного типа элементов
template<class T>
void run_benchmarks_for_type(
    BenchmarkSuite& suite,
    const std::string& element_type,
    const std::vector<int>& sizes,
    const std::vector<std::string>& matrix_types
) {
    std::cout << "\n=== Running benchmarks for " << element_type << " ===\n";

    for (const auto& matrix_type : matrix_types) {
        std::cout << "\nMatrix type: " << matrix_type << "\n";

        for (int size : sizes) {
            std::cout << "  Size " << size << "x" << size << "...\n";

            // Генерируем матрицы
            auto [A, B] = generate_matrices<T>(matrix_type, size, 42);

            // Вычисляем эталонный результат (naive)
            Matrix<T> C_reference;
            mul_naive(A, B, C_reference, nullptr);

            // Список алгоритмов для тестирования
            struct AlgoTest {
                std::string name;
                std::function<void(const Matrix<T>&, const Matrix<T>&, Matrix<T>&, OpCounter*)> func;
                bool only_4x4;  // Алгоритм работает только для размера 4
                bool only_pow2; // Алгоритм работает только для степеней 2
            };

            std::vector<AlgoTest> algorithms = {
                {"naive", wrapper_naive<T>, false, false},
                {"strassen", wrapper_strassen<T>, false, true},
                {"strassen_4x4", wrapper_strassen_4x4<T>, true, false},
                {"winograd_4x4", wrapper_winograd_4x4<T>, true, false},
                {"alphaevolve_4x4", wrapper_alphaevolve_4x4<T>, true, false},
                {"blocked_naive", wrapper_blocked_naive<T>, false, false},
                {"blocked_winograd", wrapper_blocked_winograd<T>, false, false},
                {"blocked_alphaevolve", wrapper_blocked_alphaevolve<T>, false, false},
                {"blocked_strassen", wrapper_blocked_strassen<T>, false, false}
            };

            for (const auto& algo : algorithms) {
                // Пропускаем 4x4 алгоритмы для других размеров
                if (algo.only_4x4 && size != 4) continue;

                // Пропускаем алгоритмы только для степеней 2
                if (algo.only_pow2) {
                    int pow2 = 1;
                    while (pow2 < size) pow2 *= 2;
                    if (pow2 != size) continue;
                }

                try {
                    auto result = run_single_benchmark<T>(
                        algo.name,
                        matrix_type,
                        element_type,
                        size,
                        A, B,
                        algo.func,
                        &C_reference
                    );

                    suite.add_result(result);
                    std::cout << "    " << algo.name << ": "
                              << std::fixed << std::setprecision(2) << result.time_ms << " ms"
                              << ", error: " << std::scientific << result.correctness_error << "\n";

                } catch (const std::exception& e) {
                    std::cerr << "    " << algo.name << ": FAILED (" << e.what() << ")\n";
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    std::cout << "Matrix Multiplication Benchmark Suite\n";
    std::cout << "======================================\n";

    BenchmarkSuite suite;

    // Конфигурация бенчмарков
    std::vector<int> sizes = {4, 8, 16, 64, 256};  // Размеры матриц
    std::vector<std::string> matrix_types = {"random", "symmetric"};  // Типы матриц

    // Запускаем бенчмарки для double
    run_benchmarks_for_type<double>(suite, "double", sizes, matrix_types);

    // Запускаем бенчмарки для complex<double>
    run_benchmarks_for_type<std::complex<double>>(suite, "complex", sizes, matrix_types);

    // Сохраняем результаты
    std::cout << "\n=== Saving results ===\n";
    suite.save_csv("benchmark_results.csv");

    // Выводим сводку
    std::cout << "\nTotal benchmarks run: " << suite.size() << "\n";
    std::cout << "Results saved to benchmark_results.csv\n";

    return 0;
}
