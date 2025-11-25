#pragma once
#include <vector>
using Matrix = std::vector<std::vector<double>>;
Matrix zeros(int n, int m);
Matrix identity(int n);
Matrix multiply(const Matrix& A, const Matrix& B);
Matrix subtract(const Matrix& A, const Matrix& B);
Matrix scalarMultiply(const Matrix& A, double s);
void printMatrix(const Matrix& A);
