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

#include "HelperFunctionsLab2.h"
#include "RecursiveLUFactorization.h"
#include "RecursiveGaussElimination.h"
#include "RecursivelyInversingMatrix.h"

const int GAUSS_BLOCK_SIZE = 8;

struct BenchmarkResult
{
    std::string dimensions;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
};

std::vector<int> getTestSizes(int maxSize)
{
    std::vector<int> sizes;
    for (int i = 1; i <= maxSize; ++i)
    {
        sizes.push_back(i);
    }
    // for (int i = 32; i <= 128; i *= 2)
    //     sizes.push_back(i);
    // sizes.push_back(200);
    // sizes.push_back(256);
    // sizes.push_back(400);
    // sizes.push_back(512);
    // sizes.push_back(768);
    // sizes.push_back(1000);
    return sizes;
}

std::vector<double> matrixVectorMultiply(const Matrix &A, const std::vector<double> &x)
{
    int n = A.size();
    if (n == 0 || A[0].size() != x.size())
    {
        throw std::runtime_error("Incompatible dimensions in matrixVectorMultiply");
    }
    int k = x.size();
    std::vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < k; ++j)
        {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

void run_Gauss_Benchmark(std::vector<BenchmarkResult> &results, MultiplyAlgorithm algo, const std::string &algoName, int maxSize)
{
    std::cout << "\n--- Running Test: Gauss Elimination (Multiplying: " << algoName << ") ---" << std::endl;
    std::vector<int> sizes = getTestSizes(maxSize);

    for (int n : sizes)
    {
        std::string dim_str = std::to_string(n) + "x" + std::to_string(n);
        std::cout << "  Testing size: " << dim_str << "..." << std::flush;

        unsigned long long ops = 0;
        double peakMemKB = 0;
        std::chrono::duration<double, std::milli> duration(0);

        resetPeakAllocation();

        try
        {
            auto start = std::chrono::high_resolution_clock::now();

            Matrix A = createMatrix(n, n, true);
            std::vector<double> x_expected(n, 1.0);
            std::vector<double> b_expected = matrixVectorMultiply(A, x_expected);

            std::vector<double> x_result = solve_block_gauss(A, b_expected, ops, GAUSS_BLOCK_SIZE, algo);

            auto end = std::chrono::high_resolution_clock::now();
            duration = end - start;

            peakMemKB = getPeakAllocationKB();
        }
        catch (const std::exception &e)
        {
            std::cerr << "ERROR for N=" << n << ": " << e.what() << std::endl;
        }

        results.push_back({dim_str, "Gauss_" + algoName, ops, duration.count(), peakMemKB});
        std::cout << " done." << std::endl;
    }
}

void run_LU_Benchmark(std::vector<BenchmarkResult> &results, MultiplyAlgorithm algo, const std::string &algoName, int maxSize)
{
    std::cout << "\n--- Running Test: LU Factorization (Multiplying: " << algoName << ") ---" << std::endl;
    std::vector<int> sizes = getTestSizes(maxSize);

    for (int n : sizes)
    {
        std::string dim_str = std::to_string(n) + "x" + std::to_string(n);
        std::cout << "  Testing size: " << dim_str << "..." << std::flush;

        unsigned long long ops = 0;
        double peakMemKB = 0;
        std::chrono::duration<double, std::milli> duration(0);

        resetPeakAllocation();

        try
        {
            auto start = std::chrono::high_resolution_clock::now();

            Matrix A = createMatrix(n, n, true);

            LU_Result lu_result = recursive_lu_factorization(A, ops, algo);

            auto end = std::chrono::high_resolution_clock::now();
            duration = end - start;

            peakMemKB = getPeakAllocationKB();
        }
        catch (const std::exception &e)
        {
            std::cerr << "ERROR for N=" << n << ": " << e.what() << std::endl;
        }

        results.push_back({dim_str, "LU_" + algoName, ops, duration.count(), peakMemKB});
        std::cout << " done." << std::endl;
    }
}

void run_Invert_Benchmark(std::vector<BenchmarkResult> &results, MultiplyAlgorithm algo, const std::string &algoName, int maxSize)
{
    std::cout << "\n--- Running Test: Matrix Inversion (Multiplying: " << algoName << ") ---" << std::endl;
    std::vector<int> sizes = getTestSizes(maxSize);

    for (int n : sizes)
    {
        std::string dim_str = std::to_string(n) + "x" + std::to_string(n);
        std::cout << "  Testing size: " << dim_str << "..." << std::flush;

        unsigned long long ops = 0;
        double peakMemKB = 0;
        std::chrono::duration<double, std::milli> duration(0);

        resetPeakAllocation();

        try
        {
            auto start = std::chrono::high_resolution_clock::now();

            Matrix A = createMatrix(n, n, true);

            Matrix A_inv = recursive_invert(A, ops, algo);

            auto end = std::chrono::high_resolution_clock::now();
            duration = end - start;

            peakMemKB = getPeakAllocationKB();
        }
        catch (const std::exception &e)
        {
            std::cerr << "ERROR for N=" << n << ": " << e.what() << std::endl;
        }

        results.push_back({dim_str, "Invert_" + algoName, ops, duration.count(), peakMemKB});
        std::cout << " done." << std::endl;
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

    csvfile << "Dimensions,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << "\"" << res.dimensions << "\","
                << res.algorithm << ","
                << res.operations << ","
                << res.duration_ms << ","
                << res.memory_kb << "\n";
    }
    csvfile.close();
    std::cout << "\nResults saved to file: " << filename << std::endl;
}

void printUsage()
{
    std::cout << "Error: Invalid argument." << std::endl;
    std::cout << "Usage: ./Test.exe [TEST_TYPE]" << std::endl;
    std::cout << "Available test types:" << std::endl;
    std::cout << "  gauss     (Tests Gaussian Elimination with Binet and Strassen)" << std::endl;
    std::cout << "  lu        (Tests LU Factorization with Binet and Strassen)" << std::endl;
    std::cout << "  invert    (Tests Matrix Inversion with Binet and Strassen)" << std::endl;
    std::cout << "  all       (Runs all the above tests)" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printUsage();
        return 1;
    }

    std::string testToRun = argv[1];
    int maxSize = std::stoi(argv[2]);
    std::vector<BenchmarkResult> results;
    std::string outputFilename;

    if (testToRun == "gauss")
    {
        run_Gauss_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_Gauss_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        outputFilename = "benchmark_gauss.csv";
    }
    else if (testToRun == "lu")
    {
        run_LU_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_LU_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        outputFilename = "benchmark_lu.csv";
    }
    else if (testToRun == "invert")
    {
        run_Invert_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_Invert_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        outputFilename = "benchmark_invert.csv";
    }
    else if (testToRun == "all")
    {
        run_Gauss_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_Gauss_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        run_LU_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_LU_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        run_Invert_Benchmark(results, MultiplyAlgorithm::BINET, "Binet", maxSize);
        run_Invert_Benchmark(results, MultiplyAlgorithm::STRASSEN, "Strassen", maxSize);
        outputFilename = "benchmark_all.csv";
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

    std::cout << "\nTests " << testToRun << " completed." << std::endl;
    return 0;
}