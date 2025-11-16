#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <windows.h>
#include <psapi.h>
#include <fstream>
#include <cstdio>
#include <iomanip>

#include "helperFunctions.h"

using Matrix = std::vector<std::vector<double>>;

void recursiveMultiply(Matrix &C, int rC, int cC,
                       const Matrix &A, int rA, int cA,
                       const Matrix &B, int rB, int cB,
                       int m, int k, int p, unsigned long long &op_count)
{
    if (m <= 2 || k <= 2 || p <= 2)
    {
        iterativeMultiply_inplace(C, rC, cC, A, rA, cA, B, rB, cB, m, k, p, op_count);
        return;
    }

    int m_split = m / 2;
    int k_split = k / 2;
    int p_split = p / 2;

    int m_rem = m - m_split;
    int k_rem = k - k_split;
    int p_rem = p - p_split;

    Matrix c11_p1 = createMatrix(m_split, p_split, false);
    Matrix c11_p2 = createMatrix(m_split, p_split, false);
    Matrix c12_p1 = createMatrix(m_split, p_rem, false);
    Matrix c12_p2 = createMatrix(m_split, p_rem, false);
    Matrix c21_p1 = createMatrix(m_rem, p_split, false);
    Matrix c21_p2 = createMatrix(m_rem, p_split, false);
    Matrix c22_p1 = createMatrix(m_rem, p_rem, false);
    Matrix c22_p2 = createMatrix(m_rem, p_rem, false);

    // 1. C11_p1 = A11 * B11
    recursiveMultiply(c11_p1, 0, 0, A, rA, cA, B, rB, cB, m_split, k_split, p_split, op_count);
    // 2. C11_p2 = A12 * B21
    recursiveMultiply(c11_p2, 0, 0, A, rA, cA + k_split, B, rB + k_split, cB, m_split, k_rem, p_split, op_count);

    // 3. C12_p1 = A11 * B12
    recursiveMultiply(c12_p1, 0, 0, A, rA, cA, B, rB, cB + p_split, m_split, k_split, p_rem, op_count);
    // 4. C12_p2 = A12 * B22
    recursiveMultiply(c12_p2, 0, 0, A, rA, cA + k_split, B, rB + k_split, cB + p_split, m_split, k_rem, p_rem, op_count);

    // 5. C21_p1 = A21 * B11
    recursiveMultiply(c21_p1, 0, 0, A, rA + m_split, cA, B, rB, cB, m_rem, k_split, p_split, op_count);
    // 6. C21_p2 = A22 * B21
    recursiveMultiply(c21_p2, 0, 0, A, rA + m_split, cA + k_split, B, rB + k_split, cB, m_rem, k_rem, p_split, op_count);

    // 7. C22_p1 = A21 * B12
    recursiveMultiply(c22_p1, 0, 0, A, rA + m_split, cA, B, rB, cB + p_split, m_rem, k_split, p_rem, op_count);
    // 8. C22_p2 = A22 * B22
    recursiveMultiply(c22_p2, 0, 0, A, rA + m_split, cA + k_split, B, rB + k_split, cB + p_split, m_rem, k_rem, p_rem, op_count);

    // C11 = c11_p1 + c11_p2
    addMatrices_inplace(C, rC, cC, c11_p1, 0, 0, c11_p2, 0, 0, m_split, p_split, op_count);

    // C12 = c12_p1 + c12_p2
    addMatrices_inplace(C, rC, cC + p_split, c12_p1, 0, 0, c12_p2, 0, 0, m_split, p_rem, op_count);

    // C21 = c21_p1 + c21_p2
    addMatrices_inplace(C, rC + m_split, cC, c21_p1, 0, 0, c21_p2, 0, 0, m_rem, p_split, op_count);

    // C22 = c22_p1 + c22_p2
    addMatrices_inplace(C, rC + m_split, cC + p_split, c22_p1, 0, 0, c22_p2, 0, 0, m_rem, p_rem, op_count);
}

Matrix multiply_recursive_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    // op_count = 0;
    size_t m = A.size();
    size_t k = A[0].size();
    size_t p = B[0].size();

    Matrix C = createMatrix(m, p, false);

    recursiveMultiply(C, 0, 0, A, 0, 0, B, 0, 0, m, k, p, op_count);

    return C;
}

void multiply_binet_inplace(Matrix &C, int rC, int cC,
                            const Matrix &A, int rA, int cA,
                            const Matrix &B, int rB, int cB,
                            int m, int k, int p, unsigned long long &op_count)
{
    recursiveMultiply(C, rC, cC, A, rA, cA, B, rB, cB, m, k, p, op_count);
}

/*
int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "Recursive Matrix Multiplication - Binet" << std::endl;
    std::vector<int> vector_matrices = {2, 3, 5, 7, 20, 50, 100, 200};

    std::cout << " We consider following matrices sizes: " << std::endl;
    for (int n : vector_matrices)
    {
        printf("%d ", n);
    }
    std::cout << std::endl;

    unsigned long long general_op_count = 0;
    std::vector<resultsinCSV> results = {};

    double baselineMemoryKB = getPeakPrivateUsageKB();

    for (int n : vector_matrices)
    {
        Matrix A = createMatrix(static_cast<size_t>(n), static_cast<size_t>(n), true);
        Matrix B = createMatrix(static_cast<size_t>(n), static_cast<size_t>(n), true);

        auto start = std::chrono::high_resolution_clock::now();

        Matrix C = multiply_recursive_wrapper(A, B, general_op_count);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;

        std::cout << "Matrix size: " << n << "x" << n
                  << ", Operations count: " << general_op_count
                  << ", Duration: " << duration.count() << " ms" << std::endl;

        double peakMemoryKB = getPeakPrivateUsageKB();
        double algorithmMemoryKB = peakMemoryKB - baselineMemoryKB;

        printf("Memory Used (Algorithm): %.2f KB\n", algorithmMemoryKB);

        results.push_back({n, "Recursive Binet", general_op_count, duration.count(), algorithmMemoryKB});

        std::cout << std::endl
                  << std::endl;
        general_op_count = 0;
    }

    std::ofstream csvfile("matrix_multiplication_results_BINET.csv");
    csvfile << "Size,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.size << "," << res.algorithm << "," << res.operations << ","
                << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    return 0;
}*/