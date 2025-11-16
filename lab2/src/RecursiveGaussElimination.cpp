#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;

const double EPS = 1e-12;

std::vector<double> solve_recursive_gauss(Matrix A, std::vector<double> b, long long &flop_count)
{
    int n = A.size();
    if (n == 0)
        return {};
    if (n == 1)
    {
        if (fabs(A[0][0]) < EPS)
            throw std::runtime_error("Zero pivot in base case");

        flop_count++;
        return std::vector<double>{b[0] / A[0][0]};
    }

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
        std::swap(A[0], A[piv]);
        std::swap(b[0], b[piv]);
    }

    // --- Faza eliminacji ---
    for (int i = 1; i < n; i++)
    {
        // 1 operacja: A[i][0] / A[0][0]
        double factor = A[i][0] / A[0][0];
        flop_count++; // <-- ZLICZANIE (1 dzielenie)

        A[i][0] = 0.0;
        for (int j = 1; j < n; j++)
        {
            // 2 operacje: (factor * A[0][j]) i (... -= ...)
            A[i][j] -= factor * A[0][j];
            flop_count += 2; // <-- ZLICZANIE (1 mnożenie, 1 odejmowanie)
        }

        // 2 operacje: (factor * b[0]) i (... -= ...)
        b[i] -= factor * b[0];
        flop_count += 2; // <-- ZLICZANIE (1 mnożenie, 1 odejmowanie)
    }

    // --- Przygotowanie do rekurencji (kopiowanie, brak FLOPS) ---
    Matrix A_sub(n - 1, std::vector<double>(n - 1));
    std::vector<double> b_sub(n - 1);
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
            A_sub[i - 1][j - 1] = A[i][j];
        b_sub[i - 1] = b[i];
    }

    // Wywołanie rekurencyjne: 'flop_count' jest aktualizowany wewnątrz
    std::vector<double> x_sub = solve_recursive_gauss(A_sub, b_sub, flop_count);

    // --- Faza podstawiania wstecznego ---
    std::vector<double> x(n);
    double s = 0.0;
    for (int j = 1; j < n; j++)
    {
        // 2 operacje: (A[0][j] * x_sub[j - 1]) i (s += ...)
        s += A[0][j] * x_sub[j - 1];
        flop_count += 2; // <-- ZLICZANIE (1 mnożenie, 1 dodawanie)
    }

    // 2 operacje: (b[0] - s) i (... / A[0][0])
    x[0] = (b[0] - s) / A[0][0];
    flop_count += 2; // <-- ZLICZANIE (1 odejmowanie, 1 dzielenie)

    for (int i = 1; i < n; i++)
        x[i] = x_sub[i - 1];
    return x;
}

// ----- Funkcje pomocnicze (bez zmian) -----

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

// ----- Zmodyfikowana funkcja main -----

int main()
{
    Matrix A = {
        {2, 1, -1},
        {-3, -1, 2},
        {-2, 1, 2}};
    std::vector<double> b = {8, -11, -3};

    std::cout << "### Test 1: Rozwiazywanie ukladu Ax = b ###" << std::endl;
    std::cout << "Macierz A:" << std::endl;
    print_matrix(A);
    std::cout << "Wektor b:" << std::endl;
    print_vector(b);

    // Inicjujemy licznik FLOPS
    long long total_flops = 0;

    try
    {
        // Przekazujemy licznik do funkcji
        std::vector<double> x = solve_recursive_gauss(A, b, total_flops);
        std::cout << "\nRozwiazanie x:" << std::endl;
        print_vector(x);

        // Wyświetlamy wynik zliczania
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

    return 0;
}