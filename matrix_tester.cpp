/*
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <windows.h>
#include <psapi.h>

#include "matrix_ai.cpp"
#include "matrix_Strassen.cpp"

using Matrix = std::vector<std::vector<double>>;

bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon = 1e-9)
{
    if (A.size() != B.size())
        return false;
    if (A.empty() || B.empty())
        return true; // Dwie puste macierze są równe
    if (A[0].size() != B[0].size())
        return false;

    int rows = A.size();
    int cols = A[0].size();

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            // Sprawdź, czy absolutna różnica jest większa niż tolerancja
            if (std::abs(A[i][j] - B[i][j]) > epsilon)
            {
                // Jeśli chcesz zobaczyć, gdzie wystąpił błąd:
                // std::cerr << "Blad na pozycji [" << i << "][" << j << "]: "
                //           << A[i][j] << " != " << B[i][j] << std::endl;
                return false;
            }
        }
    }
    return true; // Wszystkie elementy są w granicach tolerancji
}

int main(int argc, char *argv[])
{
    std::cout << "Strassen's Matrix Multiplication - TEST POPRAWNOSCI" << std::endl;
    // Możesz użyć mniejszych rozmiarów do szybkiego testowania
    std::vector<int> vector_matrices = {2, 3, 5, 7, 16, 32, 50, 64};

    std::cout << " Rozwazane rozmiary macierzy: ";
    for (int n : vector_matrices)
    {
        printf("%d ", n);
    }
    std::cout << std::endl
              << std::endl;

    unsigned long long op_count_strassen = 0;
    unsigned long long op_count_iterative = 0;
    std::vector<resultsinCSV> results = {}; // Nadal zbieramy wyniki wydajności

    for (int n : vector_matrices)
    {
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        // --- TESTOWANIE POPRAWNOŚCI ---
        std::cout << "--- Rozmiar " << n << "x" << n << " ---" << std::endl;

        // 1. Oblicz wynik algorytmem Strassena
        op_count_strassen = 0;
        Matrix C_strassen = multiply_strassen_wrapper(A, B, op_count_strassen);

        // 2. Oblicz wynik algorytmem wzorcowym (iteracyjnym)
        op_count_iterative = 0;
        Matrix C_benchmark = iterativeMultiply(A, B, op_count_iterative);

        // 3. Porównaj wyniki
        if (areMatricesEqual(C_strassen, C_benchmark))
        {
            std::cout << "Wynik: TEST ZALICZONY (PASS)" << std::endl;
        }
        else
        {
            std::cout << "Blad: TEST NIEZALICZONY (FAIL)" << std::endl;
            // Możesz tu zatrzymać program lub wydrukować macierze, by zobaczyć błąd
            // printMatrix(C_strassen);
            // std::cout << "---" << std::endl;
            // printMatrix(C_benchmark);
            // break; // Zatrzymaj pętlę po pierwszym błędzie
        }

        // --- POMIAR WYDAJNOŚCI (jak poprzednio) ---
        auto start = std::chrono::high_resolution_clock::now();
        multiply_strassen_wrapper(A, B, op_count_strassen); // Mierzymy czas ponownego wywołania
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;

        std::cout << "Operacje (Strassen): " << op_count_strassen
                  << ", Czas: " << duration.count() << " ms" << std::endl;

        double memory = printMemoryUsage();
        results.emplace_back(n, "Strassen", op_count_strassen, duration.count(), memory);

        std::cout << std::endl;
    }

    // Zapis do CSV (jak poprzednio)
    std::ofstream csvfile("matrix_multiplication_results_STRASSEN.csv");
    // ... (reszta kodu zapisu do CSV) ...
    csvfile.close();

    return 0;
}*/