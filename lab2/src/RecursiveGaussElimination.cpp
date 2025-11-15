#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
using Matrix = std::vector<std::vector<int>>;
const double EPS = 1e-12;

std::vector<double> solve_recursive_gauss(Matrix A, std::vector<double> b)
{
    int n = A.size();
    if (n == 0)
        return {};
    if (n == 1)
    {
        if (fabs(A[0][0]) < EPS)
            throw runtime_error("Zero pivot in base case");
        return std::vector<double>{b[0] / A[0][0]};
    }

    // partial pivoting: find max abs in column 0 among rows 0..n-1
    int piv = 0;
    double maxv = fabs(A[0][0]);
    for (int i = 1; i < n; i++)
    {
        double av = fabs(A[i][0]);
        if (av > maxv)
        {
            maxv = av;
            piv = i;
        }
    }
    if (maxv < EPS)
        throw std::runtime_error("Matrix is singular (pivot ~ 0)");

    if (piv != 0)
    {
        swap(A[0], A[piv]);
        swap(b[0], b[piv]);
    }

    // eliminate below row 0
    for (int i = 1; i < n; i++)
    {
        double factor = A[i][0] / A[0][0];
        A[i][0] = 0.0;
        for (int j = 1; j < n; j++)
        {
            A[i][j] -= factor * A[0][j];
        }
        b[i] -= factor * b[0];
    }

    // build submatrix and subvector for recursion
    Matrix A_sub(n - 1, std::vector<double>(n - 1));
    std::vector<double> b_sub(n - 1);
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
            A_sub[i - 1][j - 1] = A[i][j];
        b_sub[i - 1] = b[i];
    }

    std::vector<double> x_sub = solve_recursive_gauss(A_sub, b_sub);

    std::vector<double> x(n);
    // x0 = (b0 - sum_{j=1..n-1} A0j * xj) / A00
    double s = 0.0;
    for (int j = 1; j < n; j++)
        s += A[0][j] * x_sub[j - 1];
    x[0] = (b[0] - s) / A[0][0];
    for (int i = 1; i < n; i++)
        x[i] = x_sub[i - 1];
    return x;
}

// Using the above solver to compute inverse: solve A x = e_j for each j
Matrix inverse_via_recursive_gauss(const Matrix &A)
{
    int n = A.size();
    Matrix Inv(n, std::vector<double>(n, 0.0));
    for (int j = 0; j < n; j++)
    {
        std::vector<double> e(n, 0.0);
        e[j] = 1.0;
        // we pass a copy of A to solver (it modifies)
        std::vector<double> col = solve_recursive_gauss(A, e);
        for (int i = 0; i < n; i++)
            Inv[i][j] = col[i];
    }
    return Inv;
}

int main()
{
    std::vector<int> matrix_sizes = {};
    for (int size = 2; size <= 1000; size++)
    {
        matrix_sizes.emplace_back(size);
    }

    return 0;
}