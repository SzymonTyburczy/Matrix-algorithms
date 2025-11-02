#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <windows.h>
#include <psapi.h>
#include <iomanip>
#include <fstream>
#include <cstdio>

using Matrix = std::vector<std::vector<double>>;

// moze lepiej nie uzywac tego, bo zaklamuje szybkość metody strassena
// dlatego może ustawie to na 0, mimo że sensownie byłoby ustawić na 16
// const int STRASSEN_LEAF_SIZE = 16;
const int STRASSEN_LEAF_SIZE = 0;

Matrix createMatrix(int rows, int cols, bool random = false);
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix subtractMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix iterativeMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix strassenRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count);
void joinMatrices(Matrix &C, const Matrix &c11, const Matrix &c12,
                  const Matrix &c21, const Matrix &c22, int m_split, int p_split);
Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end);

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

Matrix iterativeMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int m = A.size();
    int k = (m > 0) ? A[0].size() : 0;
    int p = (B.size() > 0) ? B[0].size() : 0;

    if (k != (int)B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    Matrix C = createMatrix(m, p);
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
            C[i][j] = sum;
            if (k > 1)
            {
                op_count += (k - 1);
            }
        }
    }
    return C;
}

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
        throw std::invalid_argument("Niezgodne wymiary macierzy do dodawania.");
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

Matrix subtractMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || A[0].empty())
    {
        Matrix C = B;
        for (auto &row : C)
            for (auto &val : row)
                val = -val;
        op_count += C.size() * C[0].size();
        return C;
    }
    if (B.empty() || B[0].empty())
        return A;

    int n = A.size();
    int m = A[0].size();

    if (n != B.size() || m != B[0].size())
    {
        throw std::invalid_argument("Niezgodne wymiary macierzy do odejmowania.");
    }

    Matrix C = createMatrix(n, m);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
            op_count++;
        }
    }
    return C;
}

Matrix strassenRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int n = A.size();
    /*
    if (n <= STRASSEN_LEAF_SIZE)
    {
        return iterativeMultiply(A, B, op_count);
    }
    */
    if (n == 1)
    {
        Matrix C = createMatrix(1, 1);
        C[0][0] = A[0][0] * B[0][0];
        op_count++;
        return C;
    }

    if (n % 2 == 0)
    {
        int n_split = n / 2;

        Matrix a11 = subMatrix(A, 0, n_split, 0, n_split);
        Matrix a12 = subMatrix(A, 0, n_split, n_split, n);
        Matrix a21 = subMatrix(A, n_split, n, 0, n_split);
        Matrix a22 = subMatrix(A, n_split, n, n_split, n);

        Matrix b11 = subMatrix(B, 0, n_split, 0, n_split);
        Matrix b12 = subMatrix(B, 0, n_split, n_split, n);
        Matrix b21 = subMatrix(B, n_split, n, 0, n_split);
        Matrix b22 = subMatrix(B, n_split, n, n_split, n);

        Matrix s1 = subtractMatrices(b12, b22, op_count);
        Matrix s2 = addMatrices(a11, a12, op_count);
        Matrix s3 = addMatrices(a21, a22, op_count);
        Matrix s4 = subtractMatrices(b21, b11, op_count);
        Matrix s5 = addMatrices(a11, a22, op_count);
        Matrix s6 = addMatrices(b11, b22, op_count);
        Matrix s7 = subtractMatrices(a12, a22, op_count);
        Matrix s8 = addMatrices(b21, b22, op_count);
        Matrix s9 = subtractMatrices(a11, a21, op_count);
        Matrix s10 = addMatrices(b11, b12, op_count);

        Matrix p1 = strassenRecursive(a11, s1, op_count);
        Matrix p2 = strassenRecursive(s2, b22, op_count);
        Matrix p3 = strassenRecursive(s3, b11, op_count);
        Matrix p4 = strassenRecursive(a22, s4, op_count);
        Matrix p5 = strassenRecursive(s5, s6, op_count);
        Matrix p6 = strassenRecursive(s7, s8, op_count);
        Matrix p7 = strassenRecursive(s9, s10, op_count);

        Matrix c11_p1 = addMatrices(p5, p4, op_count);
        Matrix c11_p2 = subtractMatrices(c11_p1, p2, op_count);
        Matrix c11 = addMatrices(c11_p2, p6, op_count);

        Matrix c12 = addMatrices(p1, p2, op_count);

        Matrix c21 = addMatrices(p3, p4, op_count);

        Matrix c22_p1 = addMatrices(p5, p1, op_count);
        Matrix c22_p2 = subtractMatrices(c22_p1, p3, op_count);
        Matrix c22 = subtractMatrices(c22_p2, p7, op_count);

        Matrix C = createMatrix(n, n);
        joinMatrices(C, c11, c12, c21, c22, n_split, n_split);

        return C;
    }

    else
    {
        int n1 = n - 1;

        Matrix A11 = subMatrix(A, 0, n1, 0, n1);
        Matrix a12 = subMatrix(A, 0, n1, n1, n); // kolumna
        Matrix a21 = subMatrix(A, n1, n, 0, n1); // wiersz
        Matrix a22 = subMatrix(A, n1, n, n1, n); // skalar 1x1

        Matrix B11 = subMatrix(B, 0, n1, 0, n1);
        Matrix b12 = subMatrix(B, 0, n1, n1, n); // kolumna
        Matrix b21 = subMatrix(B, n1, n, 0, n1); // wiersz
        Matrix b22 = subMatrix(B, n1, n, n1, n); // skalar 1x1

        Matrix C11_p1 = strassenRecursive(A11, B11, op_count); // REKURENCJA STRASSENA
        Matrix C11_p2 = iterativeMultiply(a12, b21, op_count); // O(N^2)
        Matrix C11 = addMatrices(C11_p1, C11_p2, op_count);

        Matrix C12_p1 = iterativeMultiply(A11, b12, op_count); // O(N^2)
        Matrix C12_p2 = iterativeMultiply(a12, b22, op_count); // O(N)
        Matrix C12 = addMatrices(C12_p1, C12_p2, op_count);

        Matrix C21_p1 = iterativeMultiply(a21, B11, op_count); // O(N^2)
        Matrix C21_p2 = iterativeMultiply(a22, b21, op_count); // O(N)
        Matrix C21 = addMatrices(C21_p1, C21_p2, op_count);

        Matrix C22_p1 = iterativeMultiply(a21, b12, op_count); // O(N)
        Matrix C22_p2 = iterativeMultiply(a22, b22, op_count); // O(1)
        Matrix C22 = addMatrices(C22_p1, C22_p2, op_count);

        Matrix C = createMatrix(n, n);
        joinMatrices(C, C11, C12, C21, C22, n1, n1);

        return C;
    }
}

Matrix multiply_strassen_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Niezgodne wymiary macierzy do mnozenia.");
    }
    if (A.size() != A[0].size() || B.size() != B[0].size() || A.size() != B.size())
    {
        throw std::invalid_argument("Macierze wejsciowe nie sa kwadratowe o tym samym wymiarze N.");
    }

    op_count = 0;
    return strassenRecursive(A, B, op_count);
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

    for (int n : vector_matrices)
    {
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        auto start = std::chrono::high_resolution_clock::now();

        Matrix C = multiply_strassen_wrapper(A, B, general_op_count);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;

        std::cout << "Matrix size: " << n << "x" << n
                  << ", Operations count (Strassen): " << general_op_count
                  << ", Duration: " << duration.count() << " ms" << std::endl;

        double memory = printMemoryUsage();
        results.emplace_back(n, "Strassen", general_op_count, duration.count(), memory);

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
}