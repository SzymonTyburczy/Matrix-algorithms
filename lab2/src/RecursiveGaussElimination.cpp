#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iomanip>
#include <string>
#include <chrono>

#include "helperFunctions.h"
#include "HelperFunctionsLab2.h"

#include "RecursiveLUFactorization.h"
#include "matrix_Strassen.h"
#include "matrix_Binet.h"

std::vector<double> solve_pointwise_internal(Matrix &A, std::vector<double> &b,
                                             unsigned long long &flop_count, int offset, int n_total)
{
    int current_size = n_total - offset;
    if (current_size == 1)
    {
        if (fabs(A[offset][offset]) < EPS)
            throw std::runtime_error("Zero pivot");
        flop_count++;
        return std::vector<double>{b[offset] / A[offset][offset]};
    }
    int piv = offset;
    double maxv = fabs(A[offset][offset]);
    for (int i = offset + 1; i < n_total; i++)
    {
        double av = fabs(A[i][offset]);
        if (av > maxv)
        {
            maxv = av;
            piv = i;
        }
    }
    if (maxv < EPS)
        throw std::runtime_error("Matrix singular (pointwise)");
    if (piv != offset)
    {
        std::swap(A[offset], A[piv]);
        std::swap(b[offset], b[piv]);
    }
    for (int i = offset + 1; i < n_total; i++)
    {
        double factor = A[i][offset] / A[offset][offset];
        flop_count++;
        A[i][offset] = 0.0;
        for (int j = offset + 1; j < n_total; j++)
        {
            A[i][j] -= factor * A[offset][j];
            flop_count += 2;
        }
        b[i] -= factor * b[offset];
        flop_count += 2;
    }
    std::vector<double> x_sub = solve_pointwise_internal(A, b, flop_count, offset + 1, n_total);
    std::vector<double> x(current_size);
    double s = 0.0;
    for (int j = offset + 1; j < n_total; j++)
    {
        s += A[offset][j] * x_sub[j - (offset + 1)];
        flop_count += 2;
    }
    x[0] = (b[offset] - s) / A[offset][offset];
    flop_count += 2;
    for (int i = 1; i < current_size; i++)
        x[i] = x_sub[i - 1];
    return x;
}

std::vector<double> solve_using_lu(const Matrix &L, const Matrix &U, const std::vector<double> &b,
                                   unsigned long long &flop_count)
{
    std::vector<double> y = solve_lower_triangular(L, b, flop_count);

    std::vector<double> x = solve_upper_triangular(U, y, flop_count);

    return x;
}

std::vector<double> solve_block_recursive(Matrix &A, std::vector<double> &b,
                                          unsigned long long &flop_count, int offset, int block_size,
                                          MultiplyAlgorithm algo)
{
    int n_total = A.size();
    int current_size = n_total - offset;

    if (current_size <= block_size)
    {
        return solve_pointwise_internal(A, b, flop_count, offset, n_total);
    }

    int piv = offset;
    double maxv = fabs(A[offset][offset]);
    for (int i = offset + 1; i < n_total; i++)
    {
        double av = fabs(A[i][offset]);
        if (av > maxv)
        {
            maxv = av;
            piv = i;
        }
    }
    if (maxv < EPS)
        throw std::runtime_error("Matrix singular (block pivot)");

    if (piv != offset)
    {
        std::swap(A[offset], A[piv]);
        std::swap(b[offset], b[piv]);
    }

    int b_size = block_size;
    int r_size = current_size - b_size;

    int r11 = offset, c11 = offset;
    int r12 = offset, c12 = offset + b_size;
    int r21 = offset + b_size, c21 = offset;
    int r22 = offset + b_size, c22 = offset + b_size;

    Matrix A11 = get_submatrix(A, r11, c11, b_size, b_size);

    LU_Result lu11 = recursive_lu_factorization(A11, flop_count, algo);

    Matrix X = createMatrix(b_size, r_size);

    std::vector<double> A12_col(b_size);
    for (int j = 0; j < r_size; ++j)
    {
        for (int i = 0; i < b_size; ++i)
            A12_col[i] = A[r12 + i][c12 + j];

        std::vector<double> x_col = solve_using_lu(lu11.L, lu11.U, A12_col, flop_count);

        for (int i = 0; i < b_size; ++i)
            X[i][j] = x_col[i];
    }

    std::vector<double> b1 = get_subvector(b, r11, b_size);
    std::vector<double> y = solve_using_lu(lu11.L, lu11.U, b1, flop_count);

    Matrix A21_X = createMatrix(r_size, r_size);

    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:

        multiply_strassen_inplace(A21_X, 0, 0,
                                  A, r21, c21,
                                  X, 0, 0,
                                  r_size, b_size, r_size,
                                  flop_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(A21_X, 0, 0,
                               A, r21, c21,
                               X, 0, 0,
                               r_size, b_size, r_size,
                               flop_count);
        break;
    case MultiplyAlgorithm::ITERATIVE:
    default:

        iterativeMultiply_inplace(A21_X, 0, 0,
                                  A, r21, c21,
                                  X, 0, 0,
                                  r_size, b_size, r_size,
                                  flop_count);
        break;
    }

    subtractMatrices_inplace(A, r22, c22,
                             A, r22, c22,
                             A21_X, 0, 0,
                             r_size, r_size,
                             flop_count);

    for (int i = 0; i < r_size; ++i)
    {
        double sum = 0.0;
        for (int k = 0; k < b_size; ++k)
        {
            sum += A[r21 + i][c21 + k] * y[k];
            flop_count += 2;
        }
        b[r21 + i] -= sum;
        flop_count++;
    }

    std::vector<double> x2 = solve_block_recursive(A, b, flop_count, r22, block_size, algo);
    std::vector<double> x1(b_size);

    for (int i = 0; i < b_size; ++i)
    {
        double sum = 0.0;
        for (int k = 0; k < r_size; ++k)
        {
            sum += X[i][k] * x2[k];
            flop_count += 2;
        }
        x1[i] = y[i] - sum;
        flop_count++;
    }

    std::vector<double> x_final(current_size);
    for (int i = 0; i < b_size; ++i)
        x_final[i] = x1[i];
    for (int i = 0; i < r_size; ++i)
        x_final[b_size + i] = x2[i];

    return x_final;
}

std::vector<double> solve_block_gauss(Matrix A, std::vector<double> b,
                                      unsigned long long &flop_count, int block_size,
                                      MultiplyAlgorithm algo)
{
    if (block_size <= 0)
        throw std::invalid_argument("Block size must be > 0");
    if (A.size() != b.size())
        throw std::invalid_argument("Matrix/vector size mismatch");

    return solve_block_recursive(A, b, flop_count, 0, block_size, algo);
}
/*
int main()
{
    int N = 8;
    int BLOCK_SIZE = 2;

    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::ITERATIVE;
    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::BINET;
    MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::STRASSEN;

    std::string algo_name;
    switch (algo_to_use)
    {
    case MultiplyAlgorithm::STRASSEN:
        algo_name = "Strassen";
        break;
    case MultiplyAlgorithm::BINET:
        algo_name = "Binet";
        break;
    default:
        algo_name = "Iterative";
        break;
    }

    Matrix A = createMatrix(N, N, true);
    std::vector<double> b(N);
    for (int i = 0; i < N; ++i)
    {
        double sum = 0;
        for (int j = 0; j < N; ++j)
            sum += A[i][j];
        b[i] = sum;
    }

    std::cout << "Test: Blokowa eliminacja Gaussa (rozmiar "
              << N << "x" << N << ", blok " << BLOCK_SIZE << ")" << std::endl;
    std::cout << "Uzywany algorytm mnozenia: " << algo_name << std::endl;

    unsigned long long total_flops = 0;
    double mem_usage_kb = 0;
    getPeakPrivateUsageKB();

    auto start_time = std::chrono::high_resolution_clock::now();
    try
    {
        std::vector<double> x = solve_block_gauss(A, b, total_flops, BLOCK_SIZE, algo_to_use);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_ms = end_time - start_time;
        mem_usage_kb = getPeakPrivateUsageKB();

        std::cout << "\nSolution x (should be close to 1.0):" << std::endl;
        printVector(x);
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Total number of operations (FLOPS): " << total_flops << std::endl;
        std::cout << "Execution time: " << duration_ms.count() << " ms" << std::endl;
        std::cout << "Peak memory usage: " << mem_usage_kb << " KB" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
*/