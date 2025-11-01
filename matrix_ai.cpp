#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <windows.h>
#include <psapi.h>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;

struct resultsinCSV
{
    int n_level; // 0 for 4x5, 1 for 8x10, etc.
    std::string dimensions;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
};

Matrix createMatrix(int rows, int cols, bool random = false);
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end);

bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon = 1e-9);
void joinMatrices(Matrix &C, const Matrix &c11, const Matrix &c12,
                  const Matrix &c21, const Matrix &c22, int m_split, int p_split);

double printMemoryUsage();

Matrix naiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix matrix_ai(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix multiply_ai_recursive(const Matrix &A, const Matrix &B, unsigned long long &op_count);

Matrix createMatrix(int rows, int cols, bool random)
{
    Matrix mat(rows, std::vector<double>(cols, 0.0));
    if (random)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(0.0, 1.0);
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

bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon)
{
    if (A.size() != B.size())
        return false;
    if (A.empty())
        return B.empty();
    if (A[0].size() != B[0].size())
        return false;
    int rows = A.size();
    int cols = A[0].size();
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (std::abs(A[i][j] - B[i][j]) > epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || A[0].empty())
        return B;
    if (B.empty() || B[0].empty())
        return A;
    size_t n = A.size();
    size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for addition.");
    }

    Matrix C = createMatrix(n, m);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
            op_count++;
        }
    }
    return C;
}

Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end)
{
    if (r_start >= r_end || c_start >= c_end)
    {
        return createMatrix(0, 0);
    }
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

    for (size_t i = 0; i < c11.size(); ++i)
    {
        if (c11[i].empty())
            continue;
        for (size_t j = 0; j < c11[0].size(); ++j)
        {
            C[i][j] = c11[i][j];
        }
    }
    for (size_t i = 0; i < c12.size(); ++i)
    {
        if (c12[i].empty())
            continue;
        for (size_t j = 0; j < c12[0].size(); ++j)
        {
            if (i < C.size() && (j + p_split) < C[i].size())
                C[i][j + p_split] = c12[i][j];
        }
    }

    for (size_t i = 0; i < c21.size(); ++i)
    {
        if (c21[i].empty())
            continue;
        for (size_t j = 0; j < c21[0].size(); ++j)
        {
            if ((i + m_split) < C.size() && j < C[i + m_split].size())
                C[i + m_split][j] = c21[i][j];
        }
    }
    for (size_t i = 0; i < c22.size(); ++i)
    {
        if (c22[i].empty())
            continue;
        for (size_t j = 0; j < c22[0].size(); ++j)
        {
            if ((i + m_split) < C.size() && (j + p_split) < C[i + m_split].size())
                C[i + m_split][j + p_split] = c22[i][j];
        }
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

//  Naive matrix multiplication (iterative)
Matrix naiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    size_t m = A.size();
    size_t k = (m > 0) ? A[0].size() : 0;
    size_t p = (B.size() > 0) ? B[0].size() : 0;

    if (m == 0 || p == 0)
    {
        return createMatrix(m, p);
    }
    if (k == 0 || k != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    Matrix C = createMatrix(m, p);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < p; ++j)
        {
            double sum = 0.0;
            for (size_t l = 0; l < k; ++l)
            {
                sum += A[i][l] * B[l][j];
                op_count++; // Multiplication
            }
            C[i][j] = sum;
            if (k > 1)
            {
                op_count += (k - 1); // Additions
            }
        }
    }
    return C;
}

//  matrix_ai (4x5 * 5x5 = 4x5)
Matrix matrix_ai(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    assert(A.size() == 4 && A[0].size() == 5 && "Matrix A must be 4x5");
    assert(B.size() == 5 && B[0].size() == 5 && "Matrix B must be 5x5");

    op_count = 0;
    std::vector<double> H(76);

    // Compute H values (as in article)
    H[0] = A[2][1] * (-B[1][0] - B[1][4] - B[2][0]);
    op_count += 3;
    H[1] = (A[1][1] + A[1][4] - A[2][4]) * (-B[1][4] - B[4][0]);
    op_count += 4;
    H[2] = (-A[2][0] - A[3][0] + A[3][1]) * (-B[0][0] + B[1][4]);
    op_count += 4;
    H[3] = (A[0][1] + A[0][3] + A[2][3]) * (-B[1][4] - B[3][0]);
    op_count += 4;
    H[4] = (A[0][4] + A[1][1] + A[1][4]) * (-B[1][3] + B[4][0]);
    op_count += 4;
    H[5] = (-A[1][1] - A[1][4] - A[3][4]) * (B[1][2] + B[4][0]);
    op_count += 4;
    H[6] = (-A[0][0] + A[3][0] - A[3][1]) * (B[0][0] + B[1][3]);
    op_count += 4;
    H[7] = (A[2][1] - A[2][2] - A[3][2]) * (-B[1][2] + B[2][0]);
    op_count += 4;
    H[8] = (-A[0][1] - A[0][3] + A[3][3]) * (B[1][2] + B[3][0]);
    op_count += 4;
    H[9] = (A[1][1] + A[1][4]) * B[4][0];
    op_count += 2;
    H[10] = (-A[1][0] - A[3][0] + A[3][1]) * (-B[0][0] + B[1][1]);
    op_count += 4;
    H[11] = (A[3][0] - A[3][1]) * B[0][0];
    op_count += 2;
    H[12] = (A[0][1] + A[0][3] + A[1][3]) * (B[1][1] + B[3][0]);
    op_count += 4;
    H[13] = (A[0][2] - A[2][1] + A[2][2]) * (B[1][3] + B[2][0]);
    op_count += 4;
    H[14] = (-A[0][1] - A[0][3]) * B[3][0];
    op_count += 2;
    H[15] = (-A[2][1] + A[2][2]) * B[2][0];
    op_count += 2;
    H[16] = (A[0][1] + A[0][3] - A[1][0] + A[1][1] - A[1][2] + A[1][3] - A[2][1] + A[2][2] - A[3][0] + A[3][1]) * B[1][1];
    op_count += 10;
    H[17] = A[1][0] * (B[0][0] + B[0][1] + B[4][1]);
    op_count += 3;
    H[18] = -A[1][2] * (B[2][0] + B[2][1] + B[4][1]);
    op_count += 3;
    H[19] = (-A[0][4] + A[1][0] + A[1][2] - A[1][4]) * (-B[0][0] - B[0][1] + B[0][3] - B[4][1]);
    op_count += 8;
    H[20] = (A[1][0] + A[1][2] - A[1][4]) * B[4][1];
    op_count += 3;
    H[21] = (A[0][2] - A[0][3] - A[1][3]) * (B[0][0] + B[0][1] - B[0][3] - B[2][0] - B[2][1] + B[2][3] + B[3][3]);
    op_count += 9;
    H[22] = A[0][2] * (-B[2][0] + B[2][3] + B[3][3]);
    op_count += 3;
    H[23] = A[0][4] * (-B[3][3] - B[4][0] + B[4][3]);
    op_count += 3;
    H[24] = -A[0][0] * (B[0][0] - B[0][3]);
    op_count += 2;
    H[25] = (-A[0][2] + A[0][3] + A[0][4]) * B[3][3];
    op_count += 3;
    H[26] = (A[0][2] - A[2][0] + A[2][2]) * (B[0][0] - B[0][3] + B[0][4] + B[2][4]);
    op_count += 6;
    H[27] = -A[2][3] * (-B[2][4] - B[3][0] - B[3][4]);
    op_count += 3;
    H[28] = A[2][0] * (B[0][0] + B[0][4] + B[2][4]);
    op_count += 3;
    H[29] = (A[2][0] - A[2][2] + A[2][3]) * B[2][4];
    op_count += 3;
    H[30] = (-A[0][3] - A[0][4] - A[2][3]) * (-B[3][3] - B[4][0] + B[4][3] - B[4][4]);
    op_count += 7;
    H[31] = (A[1][0] + A[3][0] + A[3][3]) * (B[0][2] - B[3][0] - B[3][1] - B[3][2]);
    op_count += 6;
    H[32] = A[3][2] * (-B[2][0] - B[2][2]);
    op_count += 2;
    H[33] = A[3][3] * (-B[0][2] + B[3][0] + B[3][2]);
    op_count += 3;
    H[34] = -A[3][4] * (B[0][2] + B[4][0] + B[4][2]);
    op_count += 3;
    H[35] = (A[1][2] - A[1][4] - A[3][4]) * (B[2][0] + B[2][1] + B[2][2] + B[4][1]);
    op_count += 6;
    H[36] = (-A[3][0] - A[3][3] + A[3][4]) * B[0][2];
    op_count += 3;
    H[37] = (-A[1][2] - A[2][0] + A[2][2] - A[2][3]) * (B[2][4] + B[3][0] + B[3][1] + B[3][4]);
    op_count += 7;
    H[38] = (-A[2][0] - A[3][0] - A[3][3] + A[3][4]) * (B[0][2] + B[4][0] + B[4][2] + B[4][4]);
    op_count += 7;
    H[39] = (-A[0][2] + A[0][3] + A[0][4] - A[3][3]) * (-B[2][0] - B[2][2] + B[2][3] + B[3][3]);
    op_count += 7;
    H[40] = (-A[0][0] + A[3][0] - A[3][4]) * (B[0][2] + B[2][0] + B[2][2] - B[2][3] + B[4][0] + B[4][2] - B[4][3]);
    op_count += 9;
    H[41] = (-A[1][0] + A[1][4] - A[2][4]) * (-B[0][0] - B[0][1] - B[0][4] + B[3][0] + B[3][1] + B[3][4] - B[4][1]);
    op_count += 10;
    H[42] = A[1][3] * (B[3][0] + B[3][1]);
    op_count += 2;
    H[43] = (A[1][2] + A[2][1] - A[2][2]) * (B[1][1] - B[2][0]);
    op_count += 4;
    H[44] = (-A[2][2] + A[2][3] - A[3][2]) * (B[2][4] + B[3][0] + B[3][2] + B[3][4] + B[4][0] + B[4][2] + B[4][4]);
    op_count += 9;
    H[45] = -A[2][4] * (-B[4][0] - B[4][4]);
    op_count += 2;
    H[46] = (A[1][0] - A[1][4] - A[2][0] + A[2][4]) * (B[0][0] + B[0][1] + B[0][4] - B[3][0] - B[3][1] - B[3][4]);
    op_count += 9;
    H[47] = (-A[1][2] + A[2][2]) * (B[1][1] + B[2][1] + B[2][4] + B[3][0] + B[3][1] + B[3][4]);
    op_count += 7;
    H[48] = (-A[0][0] - A[0][2] + A[0][3] + A[0][4] - A[1][0] - A[1][2] + A[1][3] + A[1][4]) * (-B[0][0] - B[0][1] + B[0][3]);
    op_count += 11;
    H[49] = (-A[0][3] - A[1][3]) * (B[1][1] - B[2][0] - B[2][1] + B[2][3] - B[3][1] + B[3][3]);
    op_count += 7;
    H[50] = A[1][1] * (B[1][0] + B[1][1] - B[4][0]);
    op_count += 3;
    H[51] = A[3][1] * (B[0][0] + B[1][0] + B[1][2]);
    op_count += 3;
    H[52] = -A[0][1] * (-B[1][0] + B[1][3] + B[3][0]);
    op_count += 3;
    H[53] = (A[0][1] + A[0][3] - A[1][1] - A[1][4] - A[2][1] + A[2][2] - A[3][1] + A[3][2] - A[3][3] - A[3][4]) * B[1][2];
    op_count += 10;
    H[54] = (A[0][3] - A[3][3]) * (-B[1][2] + B[2][0] + B[2][2] - B[2][3] + B[3][2] - B[3][3]);
    op_count += 7;
    H[55] = (A[0][0] - A[0][4] - A[3][0] + A[3][4]) * (B[2][0] + B[2][2] - B[2][3] + B[4][0] + B[4][2] - B[4][3]);
    op_count += 9;
    H[56] = (-A[2][0] - A[3][0]) * (-B[0][2] - B[0][4] - B[1][4] - B[4][0] - B[4][2] - B[4][4]);
    op_count += 7;
    H[57] = (-A[0][3] - A[0][4] - A[2][3] - A[2][4]) * (-B[4][0] + B[4][3] - B[4][4]);
    op_count += 6;
    H[58] = (-A[2][2] + A[2][3] - A[3][2] + A[3][3]) * (B[3][0] + B[3][2] + B[3][4] + B[4][0] + B[4][2] + B[4][4]);
    op_count += 9;
    H[59] = (A[1][4] + A[3][4]) * (B[1][2] - B[2][0] - B[2][1] - B[2][2] - B[4][1] - B[4][2]);
    op_count += 7;
    H[60] = (A[0][3] + A[2][3]) * (B[0][0] - B[0][3] + B[0][4] - B[1][4] - B[3][3] + B[3][4] - B[4][0] + B[4][3] - B[4][4]);
    op_count += 10;
    H[61] = (A[1][0] + A[3][0]) * (B[0][1] + B[0][2] + B[1][1] - B[3][0] - B[3][1] - B[3][2]);
    op_count += 7;
    H[62] = (-A[2][2] - A[3][2]) * (-B[1][2] - B[2][2] - B[2][4] - B[3][0] - B[3][2] - B[3][4]);
    op_count += 7;
    H[63] = (A[0][0] - A[0][2] - A[0][3] + A[2][0] - A[2][2] - A[2][3]) * (B[0][0] - B[0][3] + B[0][4]);
    op_count += 8;
    H[64] = (-A[0][0] + A[3][0]) * (-B[0][2] + B[0][3] + B[1][3] - B[4][0] - B[4][2] + B[4][3]);
    op_count += 7;
    H[65] = (A[0][0] - A[0][1] + A[0][2] - A[0][4] - A[1][1] - A[1][4] - A[2][1] + A[2][2] - A[3][0] + A[3][1]) * B[1][3];
    op_count += 10;
    H[66] = (A[1][4] - A[2][4]) * (B[0][0] + B[0][1] + B[0][4] - B[1][4] - B[3][0] - B[3][1] - B[3][4] + B[4][1] + B[4][4]);
    op_count += 10;
    H[67] = (A[0][0] + A[0][2] - A[0][3] - A[0][4] - A[3][0] - A[3][2] + A[3][3] + A[3][4]) * (-B[2][0] - B[2][2] + B[2][3]);
    op_count += 10;
    H[68] = (-A[0][2] + A[0][3] - A[1][2] + A[1][3]) * (-B[1][3] - B[2][0] - B[2][1] + B[2][3] - B[4][1] + B[4][3]);
    op_count += 9;
    H[69] = (A[1][2] - A[1][4] + A[3][2] - A[3][4]) * (-B[2][0] - B[2][1] - B[2][2]);
    op_count += 6;
    H[70] = (-A[2][0] + A[2][2] - A[2][3] + A[2][4] - A[3][0] + A[3][2] - A[3][3] + A[3][4]) * (-B[4][0] - B[4][2] - B[4][4]);
    op_count += 10;
    H[71] = (-A[1][0] - A[1][3] - A[3][0] - A[3][3]) * (B[3][0] + B[3][1] + B[3][2]);
    op_count += 6;
    H[72] = (A[0][2] - A[0][3] - A[0][4] + A[1][2] - A[1][3] - A[1][4]) * (B[0][0] + B[0][1] - B[0][3] + B[1][3] + B[4][1] - B[4][3]);
    op_count += 11;
    H[73] = (A[1][0] - A[1][2] + A[1][3] - A[2][0] + A[2][2] - A[2][3]) * (B[3][0] + B[3][1] + B[3][4]);
    op_count += 8;
    H[74] = -(A[0][1] + A[0][3] - A[1][1] - A[1][4] - A[2][0] + A[2][1] + A[2][3] + A[2][4] - A[3][0] + A[3][1]) * B[1][4];
    op_count += 10;
    H[75] = (A[0][2] + A[2][2]) * (-B[0][0] + B[0][3] - B[0][4] + B[1][3] + B[2][3] - B[2][4]);
    op_count += 7;

    // Compute the resulting matrix C (4x5)
    Matrix C = createMatrix(4, 5);

    C[0][0] = -H[9] + H[11] + H[13] - H[14] - H[15] + H[52] + H[4] - H[65] - H[6];
    op_count += 8;
    C[1][0] = H[9] + H[10] - H[11] + H[12] + H[14] + H[15] - H[16] - H[43] + H[50];
    op_count += 8;
    C[2][0] = H[9] - H[11] + H[14] + H[15] - H[0] + H[1] + H[2] - H[3] + H[74];
    op_count += 8;
    C[3][0] = -H[9] + H[11] - H[14] - H[15] + H[51] + H[53] - H[5] - H[7] + H[8];
    op_count += 8;
    C[0][1] = H[12] + H[14] + H[19] + H[20] - H[21] + H[22] + H[24] - H[42] + H[48] + H[49];
    op_count += 9;
    C[1][1] = -H[10] + H[11] - H[12] - H[14] - H[15] + H[16] + H[17] - H[18] - H[20] + H[42] + H[43];
    op_count += 10;
    C[2][1] = -H[15] - H[18] - H[20] - H[27] - H[28] - H[37] + H[41] + H[43] - H[46] + H[47];
    op_count += 9;
    C[3][1] = H[10] - H[11] - H[17] + H[20] - H[31] + H[32] - H[33] - H[35] + H[61] - H[69];
    op_count += 9;
    C[0][2] = H[14] + H[22] + H[23] + H[33] - H[36] + H[39] - H[40] + H[54] - H[55] - H[8];
    op_count += 9;
    C[1][2] = -H[9] + H[18] + H[31] + H[34] + H[35] + H[36] - H[42] - H[59] - H[5] - H[71];
    op_count += 9;
    C[2][2] = -H[15] - H[27] + H[32] + H[36] - H[38] + H[44] - H[45] + H[62] - H[70] - H[7];
    op_count += 9;
    C[3][2] = H[9] + H[14] + H[15] - H[32] + H[33] - H[34] - H[36] - H[53] + H[5] + H[7] - H[8];
    op_count += 10;
    C[0][3] = -H[9] + H[11] + H[13] - H[15] + H[22] + H[23] + H[24] + H[25] + H[4] - H[65] - H[6];
    op_count += 10;
    C[1][3] = H[9] + H[17] - H[18] + H[19] - H[21] - H[23] - H[25] - H[4] - H[68] + H[72];
    op_count += 9;
    C[2][3] = -H[13] + H[15] - H[22] - H[25] + H[26] + H[28] + H[30] + H[45] - H[57] + H[75];
    op_count += 9;
    C[3][3] = H[11] + H[24] + H[25] - H[32] - H[34] - H[39] + H[40] + H[64] - H[67] - H[6];
    op_count += 9;
    C[0][4] = H[14] + H[23] + H[24] + H[26] - H[27] + H[29] + H[30] - H[3] + H[60] + H[63];
    op_count += 9;
    C[1][4] = -H[9] - H[17] - H[1] - H[29] - H[37] + H[41] - H[42] + H[45] + H[66] + H[73];
    op_count += 9;
    C[2][4] = -H[9] + H[11] - H[14] + H[27] + H[28] - H[1] - H[29] - H[2] + H[45] + H[3] - H[74];
    op_count += 10;
    C[3][4] = -H[11] - H[28] + H[29] - H[33] + H[34] + H[38] + H[2] - H[44] + H[56] + H[58];
    op_count += 9;

    return C;
}

Matrix multiply_ai_recursive(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    const int M_BASE = 4;
    const int K_BASE = 5;
    const int P_BASE = 5;

    int M = A.size();
    int K = (M > 0) ? A[0].size() : 0;
    int P = (B.size() > 0) ? B[0].size() : 0;

    if (M == M_BASE && K == K_BASE && B.size() == K_BASE && P == P_BASE)
    {
        unsigned long long temp_ops = 0;
        Matrix C = matrix_ai(A, B, temp_ops);
        op_count += temp_ops;
        return C;
    }

    int m_split = M / 2;
    int k_split = K / 2;
    int p_split = P / 2;

    // Check for odd dimensions
    if (M % 2 != 0 || K % 2 != 0 || P % 2 != 0)
    {
        std::cerr << "Error: Matrix dimensions (" << M << "x" << K << ") * ("
                  << K << "x" << P << "). "
                  << "Matrix does not have appropriate dimensions for this algorithm." << std::endl;
        throw std::invalid_argument("Matrix does not have appropriate dimensions for this algorithm.");
    }

    // Subdivide A into 4 blocks
    Matrix a11 = subMatrix(A, 0, m_split, 0, k_split);
    Matrix a12 = subMatrix(A, 0, m_split, k_split, K);
    Matrix a21 = subMatrix(A, m_split, M, 0, k_split);
    Matrix a22 = subMatrix(A, m_split, M, k_split, K);

    Matrix b11 = subMatrix(B, 0, k_split, 0, p_split);
    Matrix b12 = subMatrix(B, 0, k_split, p_split, P);
    Matrix b21 = subMatrix(B, k_split, K, 0, p_split);
    Matrix b22 = subMatrix(B, k_split, K, p_split, P);

    // 8 recursive calls
    Matrix c11_p1 = multiply_ai_recursive(a11, b11, op_count);
    Matrix c11_p2 = multiply_ai_recursive(a12, b21, op_count);

    Matrix c12_p1 = multiply_ai_recursive(a11, b12, op_count);
    Matrix c12_p2 = multiply_ai_recursive(a12, b22, op_count);

    Matrix c21_p1 = multiply_ai_recursive(a21, b11, op_count);
    Matrix c21_p2 = multiply_ai_recursive(a22, b21, op_count);

    Matrix c22_p1 = multiply_ai_recursive(a21, b12, op_count);
    Matrix c22_p2 = multiply_ai_recursive(a22, b22, op_count);

    Matrix c11 = addMatrices(c11_p1, c11_p2, op_count);
    Matrix c12 = addMatrices(c12_p1, c12_p2, op_count);
    Matrix c21 = addMatrices(c21_p1, c21_p2, op_count);
    Matrix c22 = addMatrices(c22_p1, c22_p2, op_count);

    Matrix C = createMatrix(M, P);
    joinMatrices(C, c11, c12, c21, c22, m_split, p_split);

    return C;
}

int main()
{
    std::cout << "Benchmark Matrix Multiplication (Recursive AI vs Naive)" << std::endl;

    std::vector<resultsinCSV> results;
    unsigned long long op_count_naive = 0;
    unsigned long long op_count_ai = 0;

    // Loop through recursion levels (n=0 to base, n=1 to 8x10, n=2 to 16x20, etc.)
    int max_n_level = 5; // Test up to size (4*2^5) x (5*2^5) = 128x160

    for (int n = 0; n <= max_n_level; ++n)
    {
        // Compute dimensions for the given level n
        int M = 4 * static_cast<int>(std::pow(2, n));
        int K = 5 * static_cast<int>(std::pow(2, n));
        int P = 5 * static_cast<int>(std::pow(2, n));

        std::string dims_str = "(" + std::to_string(M) + "x" + std::to_string(K) + ") * (" +
                               std::to_string(K) + "x" + std::to_string(P) + ")";

        std::cout << "\n--- Testing (Level n=" << n << "): " << dims_str << " ---" << std::endl;

        Matrix A = createMatrix(M, K, true);
        Matrix B = createMatrix(K, P, true);

        op_count_naive = 0;
        Matrix C_benchmark = naiveMultiply(A, B, op_count_naive);

        op_count_ai = 0;
        Matrix C_ai;
        if (n == 0)
        {
            C_ai = matrix_ai(A, B, op_count_ai);
        }
        else
        {
            C_ai = multiply_ai_recursive(A, B, op_count_ai);
        }

        bool ai_ok = areMatricesEqual(C_benchmark, C_ai, 1e-9);
        std::cout << "Correctness: " << (ai_ok ? "PASSED" : "FAILED") << std::endl;

        if (!ai_ok)
        {
            std::cout << "Critical error, AI algorithm is not working correctly. Stopping benchmark." << std::endl;
            break;
        }

        op_count_naive = 0;
        auto start_naive = std::chrono::high_resolution_clock::now();
        naiveMultiply(A, B, op_count_naive);
        auto end_naive = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_naive = end_naive - start_naive;
        double memory_naive = printMemoryUsage(); // Memory (note: this is the memory of the entire process)

        std::cout << "Naive:     Operations: " << std::setw(12) << op_count_naive
                  << ", Time: " << std::fixed << std::setprecision(4) << duration_naive.count() << " ms" << std::endl;
        results.push_back({n, dims_str, "Naive", op_count_naive, duration_naive.count(), memory_naive});

        op_count_ai = 0;
        auto start_ai = std::chrono::high_resolution_clock::now();
        if (n == 0)
        {
            matrix_ai(A, B, op_count_ai);
        }
        else
        {
            multiply_ai_recursive(A, B, op_count_ai);
        }
        auto end_ai = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_ai = end_ai - start_ai;
        double memory_ai = printMemoryUsage();

        std::cout << "AI (Rec):  Operations: " << std::setw(12) << op_count_ai
                  << ", Time: " << std::fixed << std::setprecision(4) << duration_ai.count() << " ms" << std::endl;
        results.push_back({n, dims_str, "AI_Recursive", op_count_ai, duration_ai.count(), memory_ai});
    }

    std::ofstream csvfile("ai_recursive_benchmark.csv");
    csvfile << "n_level,Dimensions,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.n_level << ",\"" << res.dimensions << "\"," << res.algorithm << ","
                << res.operations << "," << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    std::cout << "\nBenchmark finished. Results saved to 'ai_recursive_benchmark.csv'" << std::endl;

    return 0;
}