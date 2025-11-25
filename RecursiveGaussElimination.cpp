#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "HelpFunctions.h"
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

void printVector(const Vector& v) {
    for (double x : v) std::cout << x << " ";
    std::cout << "\n";
}

// --- Rekurencyjna eliminacja Gaussa ---
Vector gaussRecursive(Matrix A, Vector b, double tol = 1e-12) {
    int n = A.size();
    if (n == 0) return {};

    // baza rekurencji
    if (n == 1) {
        if (std::fabs(A[0][0]) < tol)
            throw std::runtime_error("Macierz osobliwa (pivot ~ 0).");
        return { b[0] / A[0][0] };
    }

    // pivotowanie: znajdź największy element w kolumnie 0
    int pivot = 0;
    for (int i = 1; i < n; ++i)
        if (std::fabs(A[i][0]) > std::fabs(A[pivot][0]))
            pivot = i;

    if (pivot != 0) {
        std::swap(A[0], A[pivot]);
        std::swap(b[0], b[pivot]);
    }

    double a00 = A[0][0];
    if (std::fabs(a00) < tol)
        throw std::runtime_error("Pivot równy 0 — układ nieodwracalny.");

    // Redukcja dolnych wierszy
    Matrix A_sub(n - 1, std::vector<double>(n - 1));
    Vector b_sub(n - 1);

    for (int i = 1; i < n; ++i) {
        double m = A[i][0] / a00;
        b_sub[i - 1] = b[i] - m * b[0];
        for (int j = 1; j < n; ++j)
            A_sub[i - 1][j - 1] = A[i][j] - m * A[0][j];
    }

    // Rekurencja
    Vector x_sub = gaussRecursive(A_sub, b_sub, tol);

    // Cofnięcie
    double sum = 0.0;
    for (int j = 1; j < n; ++j)
        sum += A[0][j] * x_sub[j - 1];

    double x0 = (b[0] - sum) / a00;

    // Składamy wynik
    Vector x(n);
    x[0] = x0;
    for (int i = 1; i < n; ++i)
        x[i] = x_sub[i - 1];
    return x;
}

int main() {
    Matrix A = {
        {2, 1, -1},
        {-3, -1, 2},
        {-2, 1, 2}
    };
    Vector b = {8, -11, -3};

    try {
        Vector x = gaussRecursive(A, b);
        std::cout << "Rozwiązanie x = ";
        printVector(x);
    } catch (const std::exception& e) {
        std::cerr << "Błąd: " << e.what() << "\n";
    }
}
