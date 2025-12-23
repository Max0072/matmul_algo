//
// Система бенчмарков для анализа алгоритмов умножения матриц
//

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "structures.h"
#include "generators.h"
#include "rss.h"
#include <chrono>
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>

// Результат одного бенчмарка
struct BenchmarkResult {
    std::string algorithm;      // Название алгоритма
    std::string matrix_type;    // Тип матрицы (random/symmetric/sparse)
    std::string element_type;   // Тип элементов (double/complex)
    int size;                   // Размер матрицы (n для nxn)
    double time_ms;             // Время выполнения в миллисекундах
    uint64_t memory_bytes;      // Использованная память в байтах
    uint64_t mul_count;         // Количество умножений
    uint64_t add_count;         // Количество сложений
    double correctness_error;   // Максимальная ошибка относительно naive

    // CSV заголовок
    static std::string csv_header() {
        return "algorithm,matrix_type,element_type,size,time_ms,memory_bytes,mul_count,add_count,correctness_error";
    }

    // Вывод в CSV формате
    std::string to_csv() const {
        std::ostringstream oss;
        oss << algorithm << ","
            << matrix_type << ","
            << element_type << ","
            << size << ","
            << std::fixed << std::setprecision(6) << time_ms << ","
            << memory_bytes << ","
            << mul_count << ","
            << add_count << ","
            << std::scientific << std::setprecision(10) << correctness_error;
        return oss.str();
    }

    // Человекочитаемый вывод
    void print() const {
        std::cout << "Algorithm: " << algorithm << "\n";
        std::cout << "Matrix Type: " << matrix_type << "\n";
        std::cout << "Element Type: " << element_type << "\n";
        std::cout << "Size: " << size << "x" << size << "\n";
        std::cout << "Time: " << std::fixed << std::setprecision(3) << time_ms << " ms\n";
        std::cout << "Memory: " << memory_bytes << " bytes\n";
        std::cout << "Operations: " << mul_count << " mul, " << add_count << " add\n";
        std::cout << "Error vs Naive: " << std::scientific << correctness_error << "\n";
        std::cout << "---\n";
    }
};

// Функция для вычисления max_abs_diff между двумя матрицами
template<class T>
double compute_max_diff(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows != B.rows || A.cols != B.cols) return 1e100;

    double max_d = 0.0;
    for(int i = 0; i < A.rows; i++) {
        for(int j = 0; j < A.cols; j++) {
            double d = std::abs(A(i,j) - B(i,j));
            max_d = std::max(max_d, d);
        }
    }
    return max_d;
}

// Специализация для complex<double>
template<>
inline double compute_max_diff<std::complex<double>>(
    const Matrix<std::complex<double>>& A,
    const Matrix<std::complex<double>>& B) {

    if (A.rows != B.rows || A.cols != B.cols) return 1e100;

    double max_d = 0.0;
    for(int i = 0; i < A.rows; i++) {
        for(int j = 0; j < A.cols; j++) {
            double d = std::abs(A(i,j) - B(i,j));
            max_d = std::max(max_d, d);
        }
    }
    return max_d;
}

// Запуск одного бенчмарка
template<class T>
BenchmarkResult run_single_benchmark(
    const std::string& algo_name,
    const std::string& matrix_type,
    const std::string& element_type,
    int size,
    const Matrix<T>& A,
    const Matrix<T>& B,
    std::function<void(const Matrix<T>&, const Matrix<T>&, Matrix<T>&, OpCounter*)> multiply_func,
    const Matrix<T>* C_reference = nullptr  // Для проверки корректности
) {
    BenchmarkResult result;
    result.algorithm = algo_name;
    result.matrix_type = matrix_type;
    result.element_type = element_type;
    result.size = size;

    Matrix<T> C;
    OpCounter cnt;

    // Измеряем память до и после
    uint64_t rss_before = current_rss_bytes();

    // Измеряем время
    auto t_start = std::chrono::steady_clock::now();
    multiply_func(A, B, C, &cnt);
    auto t_end = std::chrono::steady_clock::now();

    uint64_t rss_after = current_rss_bytes();

    // Сохраняем результаты
    result.time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    result.memory_bytes = (rss_after > rss_before) ? (rss_after - rss_before) : 0;
    result.mul_count = cnt.mul;
    result.add_count = cnt.add;

    // Проверяем корректность
    if (C_reference != nullptr) {
        result.correctness_error = compute_max_diff(C, *C_reference);
    } else {
        result.correctness_error = 0.0;
    }

    return result;
}

// Генерация матриц заданного типа
template<class T>
std::pair<Matrix<T>, Matrix<T>> generate_matrices(
    const std::string& matrix_type,
    int size,
    uint64_t seed = 42
) {
    if (matrix_type == "random") {
        return {gen_random<T>(size, size, seed, -1.0, 1.0),
                gen_random<T>(size, size, seed + 1, -1.0, 1.0)};
    }
    else if (matrix_type == "symmetric") {
        return {gen_symmetric<T>(size, seed, -1.0, 1.0),
                gen_symmetric<T>(size, seed + 1, -1.0, 1.0)};
    }
    else if (matrix_type == "sparse") {
        return {gen_almost_sparse<T>(size, size, 0.9, seed, -1.0, 1.0),
                gen_almost_sparse<T>(size, size, 0.9, seed + 1, -1.0, 1.0)};
    }
    else {
        // По умолчанию random
        return {gen_random<T>(size, size, seed, -1.0, 1.0),
                gen_random<T>(size, size, seed + 1, -1.0, 1.0)};
    }
}

// Класс для управления бенчмарками
class BenchmarkSuite {
private:
    std::vector<BenchmarkResult> results;

public:
    void add_result(const BenchmarkResult& result) {
        results.push_back(result);
    }

    void print_all() const {
        std::cout << "\n=== Benchmark Results ===\n\n";
        for(const auto& r : results) {
            r.print();
        }
    }

    void save_csv(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open " << filename << " for writing\n";
            return;
        }

        // Заголовок
        file << BenchmarkResult::csv_header() << "\n";

        // Данные
        for(const auto& r : results) {
            file << r.to_csv() << "\n";
        }

        file.close();
        std::cout << "Results saved to " << filename << "\n";
    }

    void clear() {
        results.clear();
    }

    size_t size() const {
        return results.size();
    }
};

#endif // BENCHMARK_H
