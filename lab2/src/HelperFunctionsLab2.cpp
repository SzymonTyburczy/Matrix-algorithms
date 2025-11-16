#include "HelperFunctionsLab2.h"

Matrix get_submatrix(const Matrix &A, int r_start, int c_start, int num_rows, int num_cols)
{
    Matrix sub = createMatrix(num_rows, num_cols); // Używamy Twojej funkcji
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            sub[i][j] = A[r_start + i][c_start + j];
        }
    }
    return sub;
}

// Wyciąga pod-wektor (jako wektor) z b
std::vector<double> get_subvector(const std::vector<double> &b, int r_start, int num_rows)
{
    std::vector<double> sub(num_rows);
    for (int i = 0; i < num_rows; ++i)
    {
        sub[i] = b[r_start + i];
    }
    return sub;
}

void printVector(const std::vector<double> &v)
{
    for (const auto &val : v)
    {
        printf("%0.4f\n", val);
    }
}