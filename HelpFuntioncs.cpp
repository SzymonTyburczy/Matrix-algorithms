#include "HelpFunctions.h"
#include <iostream>
Matrix zeros(int n, int m) {
    return Matrix(n, std::vector<double>(m, 0.0));
}

Matrix identity(int n) {
    Matrix I = zeros(n, n);
    for (int i = 0; i < n; ++i) I[i][i] = 1.0;
    return I;
}

Matrix multiply(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = B[0].size(), k = B.size();
    Matrix C = zeros(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int t = 0; t < k; ++t)
                C[i][j] += A[i][t] * B[t][j];
    return C;
}

Matrix subtract(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = A[0].size();
    Matrix C = zeros(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

Matrix scalarMultiply(const Matrix& A, double s) {
    int n = A.size(), m = A[0].size();
    Matrix C = zeros(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] * s;
    return C;
}


void printMatrix(const Matrix& A) {
    for (const auto& row : A) {
        for (double val : row)
            std::cout << val << " ";
        std::cout << "\n";
    }
}