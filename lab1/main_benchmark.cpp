#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip> // Potrzebne dla std::setw, std::setprecision, std::fixed
#include <stdexcept>
#include <map> // Do przechowywania wyników

// Dołączenie wszystkich nagłówków z implementacjami
#include "helperFunctions.h"
#include "matrix_Strassen.h"
#include "matrix_Binet.h"
#include "matrix_ai.h"

using Matrix = std::vector<std::vector<double>>;

// Struktura do przechowywania wyników
struct BenchmarkResult
{
    std::string dimensions;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
    bool passed;
};

// Funkcja do sprawdzania poprawności
bool areMatricesEqual(const Matrix &A, const Matrix &B, double epsilon = 1e-9)
{
    if (A.size() != B.size())
        return false;
    if (A.empty())
        return B.empty();
    if (A[0].size() != B[0].size())
        return false;
    size_t rows = A.size();
    size_t cols = A[0].size();
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            if (std::abs(A[i][j] - B[i][j]) > epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

// --- Funkcja do benchmarku N x N (Strassen, Binet, Iterative) ---
// *** ZMIANA: Funkcja przyjmuje 'algorithmToRun', aby uruchomić tylko jeden test ***
void run_NxN_Benchmark(std::vector<BenchmarkResult> &results, const std::string &algorithmToRun)
{
    std::cout << "==================================================" << std::endl;
    std::cout << "   Testowanie algorytmu: " << algorithmToRun << " (Macierze N x N)" << std::endl;
    std::cout << "==================================================" << std::endl;

    std::vector<int> vector_matrices = {2, 3, 5, 7, 20, 50, 100, 200, 500, 1000};

    // Baza dla całego procesu (każdy test n=... będzie mierzony względem poprzedniego szczytu)
    // To jest teraz POPRAWNE, ponieważ mierzymy tylko jeden algorytm naraz.
    // double processBaselineKB = getPeakPrivateUsageKB();

    for (int n : vector_matrices)
    {
        std::string dim_str = std::to_string(n) + "x" + std::to_string(n);
        std::cout << "\n--- Testowanie wymiaru: " << dim_str << " ---" << std::endl;

        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        // Baza *przed* uruchomieniem algorytmu (ale po alokacji A i B)
        double matsBaselineKB = getPeakPrivateUsageKB();

        unsigned long long ops = 0;
        double peakMem = 0, deltaMem = 0;
        std::chrono::duration<double, std::milli> duration(0);

        Matrix C_benchmark; // Potrzebna do weryfikacji
        if (algorithmToRun != "iterative")
        {
            unsigned long long dummy_ops = 0;
            C_benchmark = iterativeMultiply(A, B, dummy_ops);
            // Aktualizujemy bazę, aby uwzględnić C_benchmark
            matsBaselineKB = getPeakPrivateUsageKB();
        }

        Matrix C_result;
        bool passed = true;

        auto start = std::chrono::high_resolution_clock::now();

        if (algorithmToRun == "iterative")
        {
            C_result = iterativeMultiply(A, B, ops);
            C_benchmark = C_result; // Sam jest dla siebie benchmarkiem
        }
        else if (algorithmToRun == "strassen")
        {
            C_result = multiply_strassen_wrapper(A, B, ops);
        }
        else if (algorithmToRun == "binet")
        {
            C_result = multiply_recursive_wrapper(A, B, ops);
        }

        auto end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        peakMem = getPeakPrivateUsageKB();

        // Delta = Szczyt po teście - Baza (po A i B [i C_benchmark])
        deltaMem = peakMem - matsBaselineKB;
        passed = areMatricesEqual(C_benchmark, C_result);

        results.push_back({dim_str, algorithmToRun, ops, duration.count(), deltaMem, passed});
        std::cout << std::setw(10) << algorithmToRun << ": "
                  << std::fixed << std::setprecision(4) << duration.count() << " ms, "
                  << deltaMem << " KB, Poprawny: " << (passed ? "TAK" : "NIE") << std::endl;

        // Ustawiamy bazę dla *następnej iteracji pętli* na obecny szczyt
        // processBaselineKB = peakMem;
    }
}

// --- Funkcja do benchmarku AI ---
// *** ZMIANA: Funkcja przyjmuje 'algorithmToRun', aby uruchomić tylko jeden test ***
void run_AI_Benchmark(std::vector<BenchmarkResult> &results, const std::string &algorithmToRun)
{
    std::cout << "\n\n==================================================" << std::endl;
    std::cout << "   Testowanie algorytmu: " << algorithmToRun << " (Macierze AI)" << std::endl;
    std::cout << "==================================================" << std::endl;

    int max_n_level = 7; // Test do rozmiaru (4*2^7) x (5*2^7) = 512x640
    // double processBaselineKB = getPeakPrivateUsageKB();

    for (int n = 0; n <= max_n_level; ++n)
    {
        int M = 4 * static_cast<int>(std::pow(2, n));
        int K = 5 * static_cast<int>(std::pow(2, n));
        int P = 5 * static_cast<int>(std::pow(2, n));

        std::string dim_str = "(" + std::to_string(M) + "x" + std::to_string(K) + ") * (" +
                              std::to_string(K) + "x" + std::to_string(P) + ")";
        std::cout << "\n--- Testowanie (Poziom n=" << n << "): " << dim_str << " ---" << std::endl;

        Matrix A = createMatrix(M, K, true);
        Matrix B = createMatrix(K, P, true);

        double matsBaselineKB = getPeakPrivateUsageKB();

        unsigned long long ops = 0;
        double peakMem = 0, deltaMem = 0;
        std::chrono::duration<double, std::milli> duration(0);

        Matrix C_benchmark;                   // Potrzebna do weryfikacji
        if (algorithmToRun != "iterative_ai") // Używamy innej nazwy, by nie kolidować z NxN
        {
            unsigned long long dummy_ops = 0;
            C_benchmark = iterativeMultiply(A, B, dummy_ops);
            matsBaselineKB = getPeakPrivateUsageKB();
        }

        Matrix C_result;
        bool passed = true;
        std::string algoName;

        auto start = std::chrono::high_resolution_clock::now();

        if (algorithmToRun == "iterative_ai")
        {
            algoName = "Iterative";
            C_result = iterativeMultiply(A, B, ops);
            C_benchmark = C_result;
        }
        else if (algorithmToRun == "ai")
        {
            algoName = "AI (Recursive)";
            C_result = multiply_ai_recursive_wrapper(A, B, ops);
        }

        auto end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        peakMem = getPeakPrivateUsageKB();
        deltaMem = peakMem - matsBaselineKB;
        passed = areMatricesEqual(C_benchmark, C_result);

        results.push_back({dim_str, algoName, ops, duration.count(), deltaMem, passed});
        std::cout << std::setw(10) << algoName << ": "
                  << std::fixed << std::setprecision(4) << duration.count() << " ms, "
                  << deltaMem << " KB, Poprawny: " << (passed ? "TAK" : "NIE") << std::endl;

        if (!passed)
        {
            std::cout << "KRYTYCZNY BLAD: Algorytm nie dziala poprawnie. Zatrzymywanie." << std::endl;
            break;
        }

        // processBaselineKB = peakMem;
    }
}

// --- Funkcja do zapisu CSV ---
void writeResultsToCSV(const std::string &filename, const std::vector<BenchmarkResult> &results)
{
    std::ofstream csvfile(filename);
    if (!csvfile.is_open())
    {
        std::cerr << "Nie mozna otworzyc pliku do zapisu: " << filename << std::endl;
        return;
    }

    csvfile << "Dimensions,Algorithm,Operations,Duration_ms,Memory_kb,Passed\n";
    for (const auto &res : results)
    {
        csvfile << "\"" << res.dimensions << "\","
                << res.algorithm << ","
                << res.operations << ","
                << res.duration_ms << ","
                << res.memory_kb << ","
                << (res.passed ? "Yes" : "No") << "\n";
    }
    csvfile.close();
    std::cout << "\nZapisano wyniki do pliku: " << filename << std::endl;
}

void printUsage()
{
    std::cout << "Blad: Nieprawidlowy argument." << std::endl;
    std::cout << "Uzycie: ./benchmark.exe [TYP_ALGORYTMU]" << std::endl;
    std::cout << "Dostepne typy algorytmow:" << std::endl;
    std::cout << "  iterative (dla testu NxN)" << std::endl;
    std::cout << "  strassen" << std::endl;
    std::cout << "  binet" << std::endl;
    std::cout << "  ai" << std::endl;
    std::cout << "  iterative_ai (dla testu AI)" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printUsage();
        return 1;
    }

    std::string algorithmToRun = argv[1];
    std::vector<BenchmarkResult> results;
    std::string outputFilename;

    if (algorithmToRun == "iterative")
    {
        run_NxN_Benchmark(results, "iterative");
        outputFilename = "benchmark_iterative_NxN.csv";
    }
    else if (algorithmToRun == "strassen")
    {
        run_NxN_Benchmark(results, "strassen");
        outputFilename = "benchmark_strassen.csv";
    }
    else if (algorithmToRun == "binet")
    {
        run_NxN_Benchmark(results, "binet");
        outputFilename = "benchmark_binet.csv";
    }
    else if (algorithmToRun == "ai")
    {
        run_AI_Benchmark(results, "ai");
        outputFilename = "benchmark_ai.csv";
    }
    else if (algorithmToRun == "iterative_ai")
    {
        run_AI_Benchmark(results, "iterative_ai");
        outputFilename = "benchmark_iterative_AI.csv";
    }
    else
    {
        printUsage();
        return 1;
    }

    // Zapisz wyniki do odpowiedniego pliku CSV
    if (!results.empty())
    {
        writeResultsToCSV(outputFilename, results);
    }

    std::cout << "\nBenchmark dla '" << algorithmToRun << "' zakonczony." << std::endl;
    return 0;
}