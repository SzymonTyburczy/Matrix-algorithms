#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <cassert>   // Do asercji (sprawdzania)
#include <windows.h> // Do pomiaru pamięci
#include <psapi.h>   // Do pomiaru pamięci

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, bool random = false);
bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon = 1e-9);

Matrix naiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix matrix_ai(const Matrix &A, const Matrix &B, unsigned long long &op_count);

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

//  Mnożenie naiwne (Iteracyjne)
Matrix naiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int m = A.size();
    int k = (m > 0) ? A[0].size() : 0;
    int p = (B.size() > 0) ? B[0].size() : 0;

    if (m == 0 || p == 0)
    {
        return createMatrix(m, p);
    }
    if (k == 0 || k != (int)B.size())
    {
        throw std::invalid_argument("Niezgodne wymiary macierzy do mnozenia iteracyjnego.");
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
                op_count++; // Mnożenie
            }
            C[i][j] = sum;
            if (k > 1)
            {
                op_count += (k - 1); // Dodawania
            }
        }
    }
    return C;
}

//  matrix_ai (4x5 * 5x5 = 4x5)
Matrix matrix_ai(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    assert(A.size() == 4 && A[0].size() == 5 && "Macierz A musi byc 4x5");
    assert(B.size() == 5 && B[0].size() == 5 && "Macierz B musi byc 5x5");

    op_count = 0;
    std::vector<double> H(76);

    // Obliczanie 76 wartości pośrednich H
    // Ręczne zliczanie operacji (+, -, *)
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

    // Obliczanie macierzy wynikowej C (4x5)
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

int main()
{
    unsigned long long op_count_naive = 0;
    unsigned long long op_count_ai = 0;

    std::cout << "\n--- Testowanie algorytmu matrix_ai (4x5 * 5x5) ---" << std::endl;
    Matrix A_ai = createMatrix(4, 5, true);
    Matrix B_ai = createMatrix(5, 5, true);

    Matrix C_benchmark_ai = naiveMultiply(A_ai, B_ai, op_count_naive);
    Matrix C_ai = matrix_ai(A_ai, B_ai, op_count_ai);

    bool ai_ok = areMatricesEqual(C_benchmark_ai, C_ai, 1e-9);
    std::cout << "Poprawnosc matrix_ai: " << (ai_ok ? "ZALICZONY" : "BLAD") << std::endl;
    std::cout << "Operacje (Naiwne): " << op_count_naive << std::endl;
    std::cout << "Operacje (matrix_ai): " << op_count_ai << std::endl;

    return 0;
}
