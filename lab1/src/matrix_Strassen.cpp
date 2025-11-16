#include "helperFunctions.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <fstream>
#include <iomanip>

// Normally, Strassen's algorithm is not efficient for very small matrices due to overhead.
// const int STRASSEN_LEAF_SIZE = 16;
const int STRASSEN_LEAF_SIZE = 0; // for demonstration, set to 0 to always use Strassen's method

using Matrix = std::vector<std::vector<double>>;

void strassenRecursive(Matrix &C, int rC, int cC,
                       const Matrix &A, int rA, int cA,
                       const Matrix &B, int rB, int cB,
                       int size, unsigned long long &op_count)
{
    if (size <= STRASSEN_LEAF_SIZE)
    {
        iterativeMultiply_inplace(C, rC, cC, A, rA, cA, B, rB, cB, size, size, size, op_count);
        return;
    }

    if (size == 1)
    {
        C[rC][cC] = A[rA][cA] * B[rB][cB];
        op_count++;
        return;
    }

    if (size % 2 == 0)
    {
        int n_split = size / 2;

        // s1 = B12 - B22  (B: rows rB.., columns cB+n_split..)
        Matrix s1 = subtractMatrices(B, rB, cB + n_split, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        // s2 = A11 + A12  (A: rows rA.., columns cA..)
        Matrix s2 = addMatrices(A, rA, cA, A, rA, cA + n_split, n_split, n_split, op_count);
        // s3 = A21 + A22
        Matrix s3 = addMatrices(A, rA + n_split, cA, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        // s4 = B21 - B11
        Matrix s4 = subtractMatrices(B, rB + n_split, cB, B, rB, cB, n_split, n_split, op_count);
        // s5 = A11 + A22
        Matrix s5 = addMatrices(A, rA, cA, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        // s6 = B11 + B22
        Matrix s6 = addMatrices(B, rB, cB, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        // s7 = A12 - A22
        Matrix s7 = subtractMatrices(A, rA, cA + n_split, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        // s8 = B21 + B22
        Matrix s8 = addMatrices(B, rB + n_split, cB, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        // s9 = A11 - A21
        Matrix s9 = subtractMatrices(A, rA, cA, A, rA + n_split, cA, n_split, n_split, op_count);
        // s10 = B11 + B12
        Matrix s10 = addMatrices(B, rB, cB, B, rB, cB + n_split, n_split, n_split, op_count);

        Matrix p1 = createMatrix(n_split, n_split);
        Matrix p2 = createMatrix(n_split, n_split);
        Matrix p3 = createMatrix(n_split, n_split);
        Matrix p4 = createMatrix(n_split, n_split);
        Matrix p5 = createMatrix(n_split, n_split);
        Matrix p6 = createMatrix(n_split, n_split);
        Matrix p7 = createMatrix(n_split, n_split);

        // p1 = A11 * s1   (A11 * (B12 - B22))
        strassenRecursive(p1, 0, 0, A, rA, cA, s1, 0, 0, n_split, op_count);
        // p2 = s2 * B22   ((A11 + A12) * B22)
        strassenRecursive(p2, 0, 0, s2, 0, 0, B, rB + n_split, cB + n_split, n_split, op_count);
        // p3 = s3 * B11   ((A21 + A22) * B11)
        strassenRecursive(p3, 0, 0, s3, 0, 0, B, rB, cB, n_split, op_count);
        // p4 = A22 * s4   (A22 * (B21 - B11))
        strassenRecursive(p4, 0, 0, A, rA + n_split, cA + n_split, s4, 0, 0, n_split, op_count);
        // p5 = s5 * s6    ((A11 + A22) * (B11 + B22))
        strassenRecursive(p5, 0, 0, s5, 0, 0, s6, 0, 0, n_split, op_count);
        // p6 = s7 * s8    ((A12 - A22) * (B21 + B22))
        strassenRecursive(p6, 0, 0, s7, 0, 0, s8, 0, 0, n_split, op_count);
        // p7 = s9 * s10   ((A11 - A21) * (B11 + B12))
        strassenRecursive(p7, 0, 0, s9, 0, 0, s10, 0, 0, n_split, op_count);

        // C11 = p5 + p4 - p2 + p6
        addMatrices_inplace(C, rC, cC, p5, 0, 0, p4, 0, 0, n_split, n_split, op_count);       // C11 += p5 + p4
        subtractMatrices_inplace(C, rC, cC, C, rC, cC, p2, 0, 0, n_split, n_split, op_count); // C11 -= p2
        addMatrices_inplace(C, rC, cC, C, rC, cC, p6, 0, 0, n_split, n_split, op_count);      // C11 += p6

        // C12 = p1 + p2
        addMatrices_inplace(C, rC, cC + n_split, p1, 0, 0, p2, 0, 0, n_split, n_split, op_count); // C12 = p1 + p2
        // C21 = p3 + p4
        addMatrices_inplace(C, rC + n_split, cC, p3, 0, 0, p4, 0, 0, n_split, n_split, op_count); // C21 = p3 + p4

        // C22 = p5 + p1 - p3 - p7
        addMatrices_inplace(C, rC + n_split, cC + n_split, p5, 0, 0, p1, 0, 0, n_split, n_split, op_count);                           // C22 += p5 + p1
        subtractMatrices_inplace(C, rC + n_split, cC + n_split, C, rC + n_split, cC + n_split, p3, 0, 0, n_split, n_split, op_count); // C22 -= p3
        subtractMatrices_inplace(C, rC + n_split, cC + n_split, C, rC + n_split, cC + n_split, p7, 0, 0, n_split, n_split, op_count); // C22 -= p7
    }
    else
    {
        int n1 = size - 1;

        // recursion (n-1)x(n-1) top left A11*B11
        strassenRecursive(C, rC, cC, A, rA, cA, B, rB, cB, n1, op_count);

        Matrix C11_p2 = createMatrix(n1, n1);
        // extra term: A11_right * B_bottomLeft
        iterativeMultiply_inplace(C11_p2, 0, 0, A, rA, cA + n1, B, rB + n1, cB, n1, 1, n1, op_count);

        // add to C11
        addMatrices_inplace(C, rC, cC, C, rC, cC, C11_p2, 0, 0, n1, n1, op_count);

        // C12 (top-right) += A11 * B_rightCol
        iterativeMultiply_inplace(C, rC, cC + n1, A, rA, cA, B, rB, cB + n1, n1, n1, 1, op_count);

        Matrix C12_p2 = createMatrix(n1, 1);
        // additional for C12: A11_right * B_bottomRight
        iterativeMultiply_inplace(C12_p2, 0, 0, A, rA, cA + n1, B, rB + n1, cB + n1, n1, 1, 1, op_count);

        addMatrices_inplace(C, rC, cC + n1, C, rC, cC + n1, C12_p2, 0, 0, n1, 1, op_count);

        // C21 (bottom-left) += A_bottomRow * B_left
        iterativeMultiply_inplace(C, rC + n1, cC, A, rA + n1, cA, B, rB, cB, 1, n1, n1, op_count);

        Matrix C21_p2 = createMatrix(1, n1);
        // additional for C21: A_bottomRight * B_topLeftColumn
        iterativeMultiply_inplace(C21_p2, 0, 0, A, rA + n1, cA + n1, B, rB + n1, cB, 1, 1, n1, op_count);

        addMatrices_inplace(C, rC + n1, cC, C, rC + n1, cC, C21_p2, 0, 0, 1, n1, op_count);

        // C22 (bottom-right) += A_bottomRow * B_rightCol
        iterativeMultiply_inplace(C, rC + n1, cC + n1, A, rA + n1, cA, B, rB, cB + n1, 1, n1, 1, op_count);

        Matrix C22_p2 = createMatrix(1, 1);
        iterativeMultiply_inplace(C22_p2, 0, 0, A, rA + n1, cA + n1, B, rB + n1, cB + n1, 1, 1, 1, op_count);

        C[rC + n1][cC + n1] += C22_p2[0][0];
        op_count++;
    }
}

Matrix multiply_strassen_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }
    size_t n = A.size();
    if (n != A[0].size() || n != B.size() || n != B[0].size())
    {
        throw std::invalid_argument("Input matrices are not square and of the same size N.");
    }

    // op_count = 0;
    Matrix C = createMatrix(n, n);

    strassenRecursive(C, 0, 0, A, 0, 0, B, 0, 0, n, op_count);

    return C;
}
/*
int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "Strassen's Matrix Multiplication" << std::endl;
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
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        auto start = std::chrono::high_resolution_clock::now();

        Matrix C = multiply_strassen_wrapper(A, B, general_op_count);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;

        double peakMemoryKB = getPeakPrivateUsageKB();
        double algorithmMemoryKB = peakMemoryKB - baselineMemoryKB;

        std::cout << "Matrix size: " << n << "x" << n
                  << ", Operations count (Strassen): " << general_op_count
                  << ", Duration: " << duration.count() << " ms" << std::endl;

        printf("Memory Used (Algorithm): %.2f KB\n", algorithmMemoryKB);

        results.emplace_back(n, "Strassen", general_op_count, duration.count(), algorithmMemoryKB);

        std::cout << std::endl
                  << std::endl;
        general_op_count = 0;
    }

    std::ofstream csvfile("matrix_multiplication_results_STRASSEN.csv");
    csvfile << "Size,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.size << "," << res.algorithm << "," << res.operations << ","
                << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    return 0;
}*/