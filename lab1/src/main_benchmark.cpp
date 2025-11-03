#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <map>

#include "helperFunctions.h"
#include "matrix_Strassen.h"
#include "matrix_Binet.h"
#include "matrix_ai.h"

using Matrix = std::vector<std::vector<double>>;

struct BenchmarkResult
{
    std::string dimensions;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
    bool passed;
};

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

void run_NxN_Benchmark(std::vector<BenchmarkResult> &results, const std::string &algorithmToRun)
{
    std::cout << "   Testing algorithm: " << algorithmToRun << " (Matrix N x N)" << std::endl;

    std::vector<int> vector_matrices = {2, 3, 5, 7, 20, 50, 100, 200, 500, 1000};

    for (int n : vector_matrices)
    {
        std::string dim_str = std::to_string(n) + "x" + std::to_string(n);
        std::cout << "\n--- Testing dimension: " << dim_str << " ---" << std::endl;

        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        double matsBaselineKB = getPeakPrivateUsageKB();

        unsigned long long ops = 0;
        double peakMem = 0, deltaMem = 0;
        std::chrono::duration<double, std::milli> duration(0);

        Matrix C_benchmark;
        if (algorithmToRun != "iterative")
        {
            unsigned long long dummy_ops = 0;
            C_benchmark = iterativeMultiply(A, B, dummy_ops);
            matsBaselineKB = getPeakPrivateUsageKB();
        }

        Matrix C_result;
        bool passed = true;

        auto start = std::chrono::high_resolution_clock::now();

        if (algorithmToRun == "iterative")
        {
            C_result = iterativeMultiply(A, B, ops);
            C_benchmark = C_result;
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

        deltaMem = peakMem - matsBaselineKB;
        passed = areMatricesEqual(C_benchmark, C_result);

        results.push_back({dim_str, algorithmToRun, ops, duration.count(), deltaMem, passed});
        std::cout << std::setw(10) << algorithmToRun << ": "
                  << std::fixed << std::setprecision(4) << duration.count() << " ms, "
                  << deltaMem << " KB, Correct: " << (passed ? "YES" : "NO") << std::endl;
    }
}

void run_AI_Benchmark(std::vector<BenchmarkResult> &results, const std::string &algorithmToRun)
{
    std::cout << "   Testing algorithm: " << algorithmToRun << " (AI Matrices)" << std::endl;

    int max_n_level = 8; // Test (4*2^8) x (5*2^8) = 1024x1280 matrices

    for (int n = 0; n <= max_n_level; ++n)
    {
        int M = 4 * static_cast<int>(std::pow(2, n));
        int K = 5 * static_cast<int>(std::pow(2, n));
        int P = 5 * static_cast<int>(std::pow(2, n));

        std::string dim_str = "(" + std::to_string(M) + "x" + std::to_string(K) + ") * (" +
                              std::to_string(K) + "x" + std::to_string(P) + ")";
        std::cout << "\n--- Testing (Level n=" << n << "): " << dim_str << " ---" << std::endl;

        Matrix A = createMatrix(M, K, true);
        Matrix B = createMatrix(K, P, true);

        double matsBaselineKB = getPeakPrivateUsageKB();

        unsigned long long ops = 0;
        double peakMem = 0, deltaMem = 0;
        std::chrono::duration<double, std::milli> duration(0);

        Matrix C_benchmark;
        if (algorithmToRun != "iterative_ai")
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
                  << deltaMem << " KB, Correct: " << (passed ? "YES" : "NO") << std::endl;

        if (!passed)
        {
            std::cout << "CRITICAL ERROR: Algorithm did not work correctly. Stopping." << std::endl;
            break;
        }
    }
}

void writeResultsToCSV(const std::string &filename, const std::vector<BenchmarkResult> &results)
{
    std::ofstream csvfile(filename);
    if (!csvfile.is_open())
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
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
    std::cout << "\nResults saved to file: " << filename << std::endl;
}

void printUsage()
{
    std::cout << "Error: Invalid argument." << std::endl;
    std::cout << "Usage: ./benchmark.exe [ALGORITHM_TYPE]" << std::endl;
    std::cout << "Available algorithm types:" << std::endl;
    std::cout << "  iterative (for NxN test)" << std::endl;
    std::cout << "  strassen" << std::endl;
    std::cout << "  binet" << std::endl;
    std::cout << "  ai" << std::endl;
    std::cout << "  iterative_ai (for AI test)" << std::endl;
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

    if (!results.empty())
    {
        writeResultsToCSV(outputFilename, results);
    }

    std::cout << "\nBenchmark for '" << algorithmToRun << "' completed." << std::endl;
    return 0;
}