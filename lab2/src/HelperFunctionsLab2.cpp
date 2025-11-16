#include "HelperFunctionsLab2.h"
#include "helperFunctions.h"
#include <stdexcept>
#include <cmath>
#include <cstdio>

const double EPS = 1e-12;

Matrix get_submatrix(const Matrix &A, int r_start, int c_start, int num_rows, int num_cols)
{
    Matrix sub = createMatrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            sub[i][j] = A[r_start + i][c_start + j];
        }
    }
    return sub;
}

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

std::vector<double> solve_lower_triangular(const Matrix &L, const std::vector<double> &b, unsigned long long &op_count)
{
    size_t n = L.size();
    if (n == 0)
        return {};
    if (n != b.size() || n != L[0].size())
    {
        throw std::runtime_error("Incompatible matrix/vector dimensions in solve_lower_triangular");
    }

    std::vector<double> x(n);

    for (size_t i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (size_t j = 0; j < i; ++j)
        {
            sum += L[i][j] * x[j];
            op_count += 2;
        }

        if (fabs(L[i][i]) < EPS)
        {
            throw std::runtime_error("Division by zero: singular triangular matrix.");
        }

        x[i] = (b[i] - sum) / L[i][i];
        op_count += 2;
    }

    return x;
}

std::vector<double> solve_upper_triangular(const Matrix &U, const std::vector<double> &b, unsigned long long &op_count)
{
    size_t n = U.size();
    if (n == 0)
        return {};
    if (n != b.size() || n != U[0].size())
    {
        throw std::runtime_error("Incompatible matrix/vector dimensions in solve_upper_triangular");
    }

    std::vector<double> x(n);

    for (size_t i = n; i-- > 0;)
    {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; ++j)
        {
            sum += U[i][j] * x[j];
            op_count += 2;
        }

        if (fabs(U[i][i]) < EPS)
        {
            throw std::runtime_error("Division by zero: singular triangular matrix.");
        }
        x[i] = (b[i] - sum) / U[i][i];
        op_count += 2;
    }

    return x;
}

void copy_block_to_matrix(Matrix &Target, const Matrix &Source,
                          int r_target, int c_target)
{
    for (size_t i = 0; i < Source.size(); ++i)
    {
        for (size_t j = 0; j < Source[0].size(); ++j)
        {
            Target[r_target + i][c_target + j] = Source[i][j];
        }
    }
}