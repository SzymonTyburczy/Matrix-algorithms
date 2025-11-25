#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "HelpFunctions.h"
using Matrix = std::vector<std::vector<double>>;

struct LUResult {
    Matrix L, U;
    double det;
    long long operations;
};

// Pomocnicze funkcje



Matrix multiply(const Matrix& A, const Matrix& B, long long& ops) {
    int n = A.size(), m = B[0].size(), k = B.size();
    Matrix C = zeros(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int t = 0; t < k; ++t) {
                C[i][j] += A[i][t] * B[t][j];
                ops += 2; // jedno mnożenie + jedno dodanie
            }
    return C;
}

Matrix subtract(const Matrix& A, const Matrix& B, long long& ops) {
    int n = A.size(), m = A[0].size();
    Matrix C = zeros(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) {
            C[i][j] = A[i][j] - B[i][j];
            ops++;
        }
    return C;
}

// --- Rekurencyjna faktoryzacja LU ---
LUResult luRecursive(const Matrix& A) {
    int n = A.size();
    if (n != (int)A[0].size())
        throw std::invalid_argument("Macierz musi być kwadratowa.");

    LUResult result;
    result.operations = 0;
    result.det = 1.0;

    // baza rekurencji
    if (n == 1) {
        result.L = Matrix{{1}};
        result.U = Matrix{{A[0][0]}};
        result.det = A[0][0];
        return result;
    }

    // Podział na bloki
    double a = A[0][0];
    Matrix B = zeros(1, n - 1);
    for (int j = 1; j < n; ++j) B[0][j - 1] = A[0][j];

    Matrix C = zeros(n - 1, 1);
    for (int i = 1; i < n; ++i) C[i - 1][0] = A[i][0];

    Matrix D = zeros(n - 1, n - 1);
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            D[i - 1][j - 1] = A[i][j];

    // L i U w postaci blokowej
    Matrix L = identity(n);
    for (int i = 1; i < n; ++i) {
        L[i][0] = C[i - 1][0] / a;
        result.operations++;
    }

    Matrix U = zeros(n, n);
    U[0][0] = a;
    for (int j = 1; j < n; ++j)
        U[0][j] = B[0][j - 1];

    // Oblicz Schur complement: D' = D - (C*B)/a
    Matrix CB = multiply(C, B, result.operations);
    for (int i = 0; i < n - 1; ++i)
        for (int j = 0; j < n - 1; ++j)
            CB[i][j] /= a, result.operations++;

    Matrix Dp = subtract(D, CB, result.operations);

    // Rekurencja dla D'
    LUResult sub = luRecursive(Dp);
    result.operations += sub.operations;

    // Wstawienie wyników
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            L[i][j] = sub.L[i - 1][j - 1];

    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            U[i][j] = sub.U[i - 1][j - 1];

    result.L = L;
    result.U = U;

    // Wyznacznik = iloczyn elementów diagonalnych U
    result.det = a * sub.det;
    return result;
}


int main() {
    Matrix A = {
        {4, 2, 1},
        {6, 3, 2},
        {2, 1, 3}
    };

    try {
        LUResult r = luRecursive(A);
        std::cout << "Macierz L:\n";
        printMatrix(r.L);
        std::cout << "\nMacierz U:\n";
        printMatrix(r.U);
        std::cout << "\nWyznacznik: " << r.det << "\n";
        std::cout << "Liczba operacji: " << r.operations << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Błąd: " << e.what() << "\n";
    }
}
