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

// small alias for matrix type
using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, bool random = false);
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);

// created a submatrix
Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end)
{
    int rows = r_end - r_start;
    int cols = c_end - c_start;
    Matrix sub = createMatrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            sub[i][j] = M[r_start + i][c_start + j];
        }
    }
    return sub;
}
// joining 4 submatrices into one
void joinMatrices(Matrix &C, const Matrix &c11, const Matrix &c12,
                  const Matrix &c21, const Matrix &c22, int m_split, int p_split)
{

    for (int i = 0; i < c11.size(); ++i)
    {
        for (int j = 0; j < c11[0].size(); ++j)
        {
            C[i][j] = c11[i][j];
        }
    }

    for (int i = 0; i < c12.size(); ++i)
    {
        if (!c12[i].empty())
        {
            for (int j = 0; j < c12[0].size(); ++j)
            {
                C[i][j + p_split] = c12[i][j];
            }
        }
    }

    for (int i = 0; i < c21.size(); ++i)
    {
        if (!c21[i].empty())
        {
            for (int j = 0; j < c21[0].size(); ++j)
            {
                C[i + m_split][j] = c21[i][j];
            }
        }
    }

    for (int i = 0; i < c22.size(); ++i)
    {
        if (!c22[i].empty())
        {
            for (int j = 0; j < c22[0].size(); ++j)
            {
                C[i + m_split][j + p_split] = c22[i][j];
            }
        }
    }
}

Matrix recursiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int m = A.size();
    int k = (m > 0) ? A[0].size() : 0;
    int p = (B.size() > 0) ? B[0].size() : 0;

    if (m == 0 || k == 0 || p == 0)
    {
        return createMatrix(m, p);
    }

    if (m <= 2 || k <= 2 || p <= 2)
    {
        Matrix C_iter = createMatrix(m, p);
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < p; ++j)
            {
                double sum = 0.0;
                for (int l = 0; l < k; ++l)
                {
                    sum += A[i][l] * B[l][j];
                    op_count++;
                }
                C_iter[i][j] = sum;
                if (k > 1)
                {
                    op_count += (k - 1);
                }
            }
        }
        return C_iter;
    }

    int m_split = m / 2;
    int k_split = k / 2;
    int p_split = p / 2;

    Matrix a11 = subMatrix(A, 0, m_split, 0, k_split);
    Matrix a12 = subMatrix(A, 0, m_split, k_split, k);
    Matrix a21 = subMatrix(A, m_split, m, 0, k_split);
    Matrix a22 = subMatrix(A, m_split, m, k_split, k);

    Matrix b11 = subMatrix(B, 0, k_split, 0, p_split);
    Matrix b12 = subMatrix(B, 0, k_split, p_split, p);
    Matrix b21 = subMatrix(B, k_split, k, 0, p_split);
    Matrix b22 = subMatrix(B, k_split, k, p_split, p);

    Matrix c11_p1 = recursiveMultiply(a11, b11, op_count);
    Matrix c11_p2 = recursiveMultiply(a12, b21, op_count);

    Matrix c12_p1 = recursiveMultiply(a11, b12, op_count);
    Matrix c12_p2 = recursiveMultiply(a12, b22, op_count);

    Matrix c21_p1 = recursiveMultiply(a21, b11, op_count);
    Matrix c21_p2 = recursiveMultiply(a22, b21, op_count);

    Matrix c22_p1 = recursiveMultiply(a21, b12, op_count);
    Matrix c22_p2 = recursiveMultiply(a22, b22, op_count);

    Matrix c11 = addMatrices(c11_p1, c11_p2, op_count);
    Matrix c12 = addMatrices(c12_p1, c12_p2, op_count);
    Matrix c21 = addMatrices(c21_p1, c21_p2, op_count);
    Matrix c22 = addMatrices(c22_p1, c22_p2, op_count);

    Matrix C = createMatrix(m, p);
    joinMatrices(C, c11, c12, c21, c22, m_split, p_split);

    return C;
}

Matrix multiply_recursive_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    op_count = 0;
    return recursiveMultiply(A, B, op_count);
}

// Function to add two matrices (pretty obvious i guess)
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || A[0].empty())
        return B;
    if (B.empty() || B[0].empty())
        return A;

    int n = A.size();
    int m = A[0].size();

    if (n != B.size() || m != B[0].size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for addition.");
    }

    Matrix C = createMatrix(n, m);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
            op_count++;
        }
    }
    return C;
}

Matrix createMatrix(int rows, int cols, bool random)
{
    Matrix mat(rows, std::vector<double>(cols, 0.0));
    if (random && rows > 0 && cols > 0)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(0.00000001, 1.0);

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                mat[i][j] = dis(gen);
            }
        }
    }
    return mat;
}

void printMatrix(const std::vector<std::vector<double>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (const auto &val : row)
        {
            printf("%0.4f ", val);
        }
        std::cout << std::endl;
    }
}

double printMemoryUsage()
{
    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        double peakMemoryKB = pmc.PeakWorkingSetSize / 1024.0;
        std::cout << "Memory Used: " << peakMemoryKB << " KB" << std::endl;
        return peakMemoryKB;
    }
    return 0;
}

struct resultsinCSV
{
    int size;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;

    resultsinCSV(int s, std::string algo, unsigned long long ops, double d, double mem)
        : size(s), algorithm(algo), operations(ops), duration_ms(d), memory_kb(mem) {}
};

int main(int argc, char *argv[])
{
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

    for (int n : vector_matrices)
    {
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        auto start = std::chrono::high_resolution_clock::now();

        Matrix C = multiply_recursive_wrapper(A, B, general_op_count);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;
        std::cout << "Matrix size: " << n << "x" << n << ", Operations count: " << general_op_count << ", Duration: " << duration.count() / 1000.0 << " s" << std::endl;

        double memory = printMemoryUsage();

        results.push_back({n, "Recursive Binet", general_op_count, duration.count(), memory});

        std::cout << std::endl
                  << std::endl;
        general_op_count = 0;
    }
    // Matrix A = createMatrix(3, 3, true);
    // Matrix B = createMatrix(3, 3, true);
    // std::cout << "Matrix A:" << std::endl;
    // printMatrix(A);
    // std::cout << "Matrix B:" << std::endl;
    // printMatrix(B);
    // Matrix C = multiply_recursive_wrapper(A, B, general_op_count);
    // std::cout << "Resultant Matrix C (A x B):" << std::endl;
    // printMatrix(C);
    // std::cout << "Operations count: " << general_op_count << std::endl;

    // Writing results to CSV
    std::ofstream csvfile("matrix_multiplication_results_BINET.csv");
    csvfile << "Size,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.size << "," << res.algorithm << "," << res.operations << ","
                << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    return 0;
}