#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;

const double EPS = 1e-12;

std::vector<double> solve_recursive_internal(Matrix &A, std::vector<double> &b, long long &flop_count, int offset)
{
    int n_total = A.size();
    int current_size = n_total - offset;

    if (current_size == 0)
        return {};

    if (current_size == 1)
    {
        if (fabs(A[offset][offset]) < EPS)
            throw std::runtime_error("Zero pivot in base case");

        flop_count++;
        return std::vector<double>{b[offset] / A[offset][offset]};
    }

    int piv = offset;
    double maxv = fabs(A[offset][offset]);

    for (int i = offset + 1; i < n_total; i++)
    {
        double av = fabs(A[i][offset]);
        if (av > maxv)
        {
            maxv = av;
            piv = i;
        }
    }
    if (maxv < EPS)
        throw std::runtime_error("Matrix is singular (pivot ~ 0)");

    if (piv != offset)
    {
        std::swap(A[offset], A[piv]);
        std::swap(b[offset], b[piv]);
    }

    for (int i = offset + 1; i < n_total; i++)
    {
        double factor = A[i][offset] / A[offset][offset];
        flop_count++;

        A[i][offset] = 0.0;
        for (int j = offset + 1; j < n_total; j++)
        {
            A[i][j] -= factor * A[offset][j];
            flop_count += 2;
        }

        b[i] -= factor * b[offset];
        flop_count += 2;
    }

    std::vector<double> x_sub = solve_recursive_internal(A, b, flop_count, offset + 1);

    std::vector<double> x(current_size);
    double s = 0.0;

    for (int j = offset + 1; j < n_total; j++)
    {
        s += A[offset][j] * x_sub[j - (offset + 1)];
        flop_count += 2;
    }

    x[0] = (b[offset] - s) / A[offset][offset];
    flop_count += 2;

    for (int i = 1; i < current_size; i++)
        x[i] = x_sub[i - 1];

    return x;
}

std::vector<double> solve_recursive_gauss(Matrix A, std::vector<double> b, long long &flop_count)
{
    return solve_recursive_internal(A, b, flop_count, 0);
}

void print_vector(const std::vector<double> &v)
{
    std::cout << "[ ";
    for (const auto &val : v)
    {
        std::cout << std::fixed << std::setprecision(4) << val << " ";
    }
    std::cout << "]" << std::endl;
}

void print_matrix(const Matrix &M)
{
    for (const auto &row : M)
    {
        std::cout << "| ";
        for (const auto &val : row)
        {
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << val << " ";
        }
        std::cout << " |" << std::endl;
    }
}

int main()
{
    Matrix A = {
        {2, 1, -1},
        {-3, -1, 2},
        {-2, 1, 2}};
    std::vector<double> b = {8, -11, -3};

    std::cout << "### Test 1: Rozwiazywanie ukladu Ax = b ###" << std::endl;
    std::cout << "Macierz A (przed wywolaniem):" << std::endl;
    print_matrix(A);
    std::cout << "Wektor b (przed wywolaniem):" << std::endl;
    print_vector(b);

    long long total_flops = 0;

    try
    {
        std::vector<double> x = solve_recursive_gauss(A, b, total_flops);
        std::cout << "\nRozwiazanie x:" << std::endl;
        print_vector(x);

        std::cout << "\nLaczna liczba operacji zmiennoprzecinkowych (FLOPS): "
                  << total_flops << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Blad: " << e.what() << std::endl;
    }

    std::cout << "\n"
              << std::string(40, '-') << "\n"
              << std::endl;

    // Możemy sprawdzić, że oryginalna macierz A w main nie została zmieniona
    // std::cout << "Macierz A (po wywolaniu w main):" << std::endl;
    // print_matrix(A);

    return 0;
}