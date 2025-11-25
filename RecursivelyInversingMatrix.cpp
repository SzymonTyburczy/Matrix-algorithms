#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "HelpFunctions.h"
using Matrix = std::vector<std::vector<double>>;

// --- Pomocnicze funkcje macierzowe ---


// --- Główna funkcja rekurencyjna ---

Matrix inverseRecursive(const Matrix& A, double tol = 1e-12) {
    int n = A.size();
    if (n != (int)A[0].size())
        throw std::invalid_argument("Macierz musi być kwadratowa.");

    if (n == 1) {
        if (std::fabs(A[0][0]) < tol)
            throw std::runtime_error("Macierz nieodwracalna (pivot ~ 0).");
        return Matrix{{1.0 / A[0][0]}};
    }

    // --- Podział na bloki ---
    double a = A[0][0];
    if (std::fabs(a) < tol)
        throw std::runtime_error("Pivot A[0][0] zbyt mały (dodaj pivotowanie).");

    Matrix B = zeros(1, n - 1);
    for (int j = 1; j < n; ++j) B[0][j - 1] = A[0][j];

    Matrix C = zeros(n - 1, 1);
    for (int i = 1; i < n; ++i) C[i - 1][0] = A[i][0];

    Matrix D = zeros(n - 1, n - 1);
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            D[i - 1][j - 1] = A[i][j];

    // --- Komplement Schura ---
    Matrix CB = multiply(C, B);
    Matrix S = subtract(D, scalarMultiply(CB, 1.0 / a));

    // --- Rekurencja: S^{-1} ---
    Matrix S_inv = inverseRecursive(S, tol);

    // --- Bloki odwrotnej ---
    double a_inv = 1.0 / a;

    // upper_left = a_inv + a_inv * B * S_inv * C * a_inv
    Matrix temp = multiply(B, S_inv);
    Matrix temp2 = multiply(temp, C);
    double upper_left = a_inv + a_inv * temp2[0][0] * a_inv;

    // upper_right = -a_inv * B * S_inv
    Matrix upper_right = scalarMultiply(temp, -a_inv);

    // lower_left = -S_inv * C * a_inv
    Matrix lower_left = scalarMultiply(multiply(S_inv, C), -a_inv);

    Matrix lower_right = S_inv;

    // --- Złożenie całości ---
    Matrix A_inv = zeros(n, n);
    A_inv[0][0] = upper_left;
    for (int j = 1; j < n; ++j)
        A_inv[0][j] = upper_right[0][j - 1];
    for (int i = 1; i < n; ++i)
        A_inv[i][0] = lower_left[i - 1][0];
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            A_inv[i][j] = lower_right[i - 1][j - 1];

    return A_inv;
}



int main() {
    Matrix A = {
        {4, 2, 1},
        {1, 3, 0},
        {2, 0, 5}
    };

    try {
        Matrix Ainv = inverseRecursive(A);
        std::cout << "Odwrotna macierz A:\n";
        printMatrix(Ainv);
    } catch (const std::exception& e) {
        std::cerr << "Błąd: " << e.what() << "\n";
    }

    return 0;
}
