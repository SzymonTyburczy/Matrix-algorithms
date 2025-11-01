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

// --- Aliasy i Struktury ---

using Matrix = std::vector<std::vector<double>>;

// Struktura do przechowywania wyników
struct resultsinCSV
{
    int size;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
};

// --- Deklaracje Funkcji ---

// Funkcje pomocnicze
Matrix createMatrix(int rows, int cols, bool random = false);
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix subtractMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end);
void joinMatrices(Matrix &C, const Matrix &c11, const Matrix &c12,
                  const Matrix &c21, const Matrix &c22, int m_split, int p_split);
bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon = 1e-9);
double printMemoryUsage();

// Algorytmy mnożenia
Matrix naiveMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix binetRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix strassenRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix new_AI(const Matrix &A, const Matrix &B, unsigned long long &op_count);

// Funkcje "opakowujące" (wrappery)
Matrix multiply_binet_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix multiply_strassen_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count);

// --- Implementacje Funkcji Pomocniczych ---

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

Matrix subMatrix(const Matrix &M, int r_start, int r_end, int c_start, int c_end)
{
    // Zabezpieczenie przed wycięciem poza zakres
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

// Solidna wersja joinMatrices, która obsługuje różne rozmiary podmacierzy
void joinMatrices(Matrix &C, const Matrix &c11, const Matrix &c12,
                  const Matrix &c21, const Matrix &c22, int m_split, int p_split)
{

    for (int i = 0; i < c11.size(); ++i)
    {
        if (c11[i].empty())
            continue;
        for (int j = 0; j < c11[0].size(); ++j)
        {
            C[i][j] = c11[i][j];
        }
    }
    for (int i = 0; i < c12.size(); ++i)
    {
        if (c12[i].empty())
            continue;
        for (int j = 0; j < c12[0].size(); ++j)
        {
            if (i < C.size() && (j + p_split) < C[i].size())
                C[i][j + p_split] = c12[i][j];
        }
    }
    for (int i = 0; i < c21.size(); ++i)
    {
        if (c21[i].empty())
            continue;
        for (int j = 0; j < c21[0].size(); ++j)
        {
            if ((i + m_split) < C.size() && j < C[i + m_split].size())
                C[i + m_split][j] = c21[i][j];
        }
    }
    for (int i = 0; i < c22.size(); ++i)
    {
        if (c22[i].empty())
            continue;
        for (int j = 0; j < c22[0].size(); ++j)
        {
            if ((i + m_split) < C.size() && (j + p_split) < C[i + m_split].size())
                C[i + m_split][j + p_split] = c22[i][j];
        }
    }
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
        op_count += C.size() * (C.empty() ? 0 : C[0].size());
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

// --- ALGORYTM 1: Mnożenie naiwne (Iteracyjne) ---
// (Odpowiednik `naive_multiplication` z Pythona)
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

// --- ALGORYTM 2: Mnożenie Bineta (Rekurencyjne, 8 wywołań) ---
// (Odpowiednik `binet_multiplication` z Pythona)
Matrix binetRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int m = A.size();
    int k = (m > 0) ? A[0].size() : 0;
    int p = (B.size() > 0) ? B[0].size() : 0;

    if (m == 0 || k == 0 || p == 0)
    {
        return createMatrix(m, p);
    }

    // Warunek bazowy (zamiast `A_rows == 1 ...` używamy małego rozmiaru)
    if (m <= 2 || k <= 2 || p <= 2)
    {
        return naiveMultiply(A, B, op_count);
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

    Matrix c11_p1 = binetRecursive(a11, b11, op_count);
    Matrix c11_p2 = binetRecursive(a12, b21, op_count);

    Matrix c12_p1 = binetRecursive(a11, b12, op_count);
    Matrix c12_p2 = binetRecursive(a12, b22, op_count);

    Matrix c21_p1 = binetRecursive(a21, b11, op_count);
    Matrix c21_p2 = binetRecursive(a22, b21, op_count);

    Matrix c22_p1 = binetRecursive(a21, b12, op_count);
    Matrix c22_p2 = binetRecursive(a22, b22, op_count);

    Matrix c11 = addMatrices(c11_p1, c11_p2, op_count);
    Matrix c12 = addMatrices(c12_p1, c12_p2, op_count);
    Matrix c21 = addMatrices(c21_p1, c21_p2, op_count);
    Matrix c22 = addMatrices(c22_p1, c22_p2, op_count);

    Matrix C = createMatrix(m, p);
    joinMatrices(C, c11, c12, c21, c22, m_split, p_split);

    return C;
}

Matrix multiply_binet_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Niezgodne wymiary macierzy do mnozenia.");
    }
    op_count = 0;
    return binetRecursive(A, B, op_count);
}

// --- ALGORYTM 3: Mnożenie Strassena (Rekurencyjne, 7 wywołań) ---
// (Odpowiednik `strassen_multiplication` z Pythona)
Matrix strassenRecursive(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    int n = A.size();

    // Warunek bazowy (leaf size)
    const int STRASSEN_LEAF_SIZE = 16;
    if (n <= STRASSEN_LEAF_SIZE)
    {
        return naiveMultiply(A, B, op_count);
    }

    // Przypadek parzysty
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

        Matrix c11 = addMatrices(subtractMatrices(addMatrices(p5, p4, op_count), p2, op_count), p6, op_count);
        Matrix c12 = addMatrices(p1, p2, op_count);
        Matrix c21 = addMatrices(p3, p4, op_count);
        Matrix c22 = subtractMatrices(subtractMatrices(addMatrices(p5, p1, op_count), p3, op_count), p7, op_count);

        Matrix C = createMatrix(n, n);
        joinMatrices(C, c11, c12, c21, c22, n_split, n_split);
        return C;
    }
    // Przypadek nieparzysty (logika `remove_row` z Pythona)
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
        Matrix C11_p2 = naiveMultiply(a12, b21, op_count);
        Matrix C11 = addMatrices(C11_p1, C11_p2, op_count);

        Matrix C12_p1 = naiveMultiply(A11, b12, op_count);
        Matrix C12_p2 = naiveMultiply(a12, b22, op_count);
        Matrix C12 = addMatrices(C12_p1, C12_p2, op_count);

        Matrix C21_p1 = naiveMultiply(a21, B11, op_count);
        Matrix C21_p2 = naiveMultiply(a22, b21, op_count);
        Matrix C21 = addMatrices(C21_p1, C21_p2, op_count);

        Matrix C22_p1 = naiveMultiply(a21, b12, op_count);
        Matrix C22_p2 = naiveMultiply(a22, b22, op_count);
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

// --- ALGORYTM 4: new_AI (4x5 * 5x5 = 4x5) ---
// (Bezpośrednie tłumaczenie `new_AI` z Pythona)
Matrix new_AI(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    // Sprawdzenie poprawności wymiarów
    assert(A.size() == 4 && A[0].size() == 5 && "Macierz A musi byc 4x5");
    assert(B.size() == 5 && B[0].size() == 5 && "Macierz B musi byc 5x5");

    op_count = 0;
    std::vector<double> H(77); // H = np.zeros(77)

    // Obliczanie 77 wartości pośrednich H
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
    H[76] = 0; // H[76] nie jest używane, ale jest w tablicy

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

// --- Funkcja Główna (main) ---
int main()
{
    // --- Część 1: Benchmark dla algorytmów n x n (Binet, Strassen) ---
    std::cout << "Benchmark Mnozenia Macierzy Kwadratowych (Binet vs Strassen)" << std::endl;
    std::vector<int> vector_matrices = {2, 3, 5, 7, 20, 50, 100, 200}; // Z twojego kodu Strassena

    std::cout << " Rozwazane rozmiary macierzy: ";
    for (int n : vector_matrices)
    {
        printf("%d ", n);
    }
    std::cout << std::endl
              << std::endl;

    unsigned long long op_count_naive = 0;
    unsigned long long op_count_binet = 0;
    unsigned long long op_count_strassen = 0;
    std::vector<resultsinCSV> results;

    for (int n : vector_matrices)
    {
        if (n == 0)
            continue;
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        std::cout << "--- Testowanie rozmiaru: " << n << "x" << n << " ---" << std::endl;

        // --- Test Poprawności ---
        Matrix C_benchmark = naiveMultiply(A, B, op_count_naive);
        Matrix C_binet = multiply_binet_wrapper(A, B, op_count_binet);
        Matrix C_strassen = multiply_strassen_wrapper(A, B, op_count_strassen);

        bool binet_ok = areMatricesEqual(C_benchmark, C_binet, 1e-9);
        bool strassen_ok = areMatricesEqual(C_benchmark, C_strassen, 1e-9);
        std::cout << "Poprawnosc Binet: " << (binet_ok ? "ZALICZONY" : "BLAD") << std::endl;
        std::cout << "Poprawnosc Strassen: " << (strassen_ok ? "ZALICZONY" : "BLAD") << std::endl;

        // --- Pomiar Wydajności Binet ---
        op_count_binet = 0;
        auto start_binet = std::chrono::high_resolution_clock::now();
        multiply_binet_wrapper(A, B, op_count_binet);
        auto end_binet = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_binet = end_binet - start_binet;
        double memory_binet = printMemoryUsage();
        results.push_back({n, "Binet", op_count_binet, duration_binet.count(), memory_binet});
        std::cout << "Binet:    Operacje: " << op_count_binet << ", Czas: " << duration_binet.count() << " ms" << std::endl;

        // --- Pomiar Wydajności Strassen ---
        op_count_strassen = 0;
        auto start_strassen = std::chrono::high_resolution_clock::now();
        multiply_strassen_wrapper(A, B, op_count_strassen);
        auto end_strassen = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_strassen = end_strassen - start_strassen;
        double memory_strassen = printMemoryUsage();
        results.push_back({n, "Strassen", op_count_strassen, duration_strassen.count(), memory_strassen});
        std::cout << "Strassen: Operacje: " << op_count_strassen << ", Czas: " << duration_strassen.count() << " ms" << std::endl;

        std::cout << std::endl;
    }

    // Zapis do CSV
    std::ofstream csvfile("matrix_multiplication_results.csv");
    csvfile << "Size,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.size << "," << res.algorithm << "," << res.operations << ","
                << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();
    std::cout << "Wyniki benchmarku zapisano do matrix_multiplication_results.csv" << std::endl;

    // --- Część 2: Test jednostkowy dla algorytmu new_AI (4x5 * 5x5) ---
    std::cout << "\n--- Testowanie algorytmu new_AI (4x5 * 5x5) ---" << std::endl;
    Matrix A_ai = createMatrix(4, 5, true);
    Matrix B_ai = createMatrix(5, 5, true);

    unsigned long long op_count_ai = 0;
    op_count_naive = 0;

    Matrix C_benchmark_ai = naiveMultiply(A_ai, B_ai, op_count_naive);
    Matrix C_ai = new_AI(A_ai, B_ai, op_count_ai);

    bool ai_ok = areMatricesEqual(C_benchmark_ai, C_ai, 1e-9);
    std::cout << "Poprawnosc new_AI: " << (ai_ok ? "ZALICZONY" : "BLAD") << std::endl;
    std::cout << "Operacje (Naiwne): " << op_count_naive << std::endl;
    std::cout << "Operacje (new_AI): " << op_count_ai << std::endl;

    return 0;
}
