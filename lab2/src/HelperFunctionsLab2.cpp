#include "HelperFunctionsLab2.h"
#include "helperFunctions.h"
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <new>
#include <cstdlib>
#include <atomic>

const double EPS = 1e-12;

const size_t c_prefixSize = sizeof(size_t);
std::atomic<size_t> g_currentAllocatedMemory(0);
std::atomic<size_t> g_peakAllocatedMemory(0);

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

void *allocateMemory(size_t size)
{
    size_t totalSize = size + c_prefixSize;
    void *block = std::malloc(totalSize);
    if (!block)
    {
        throw std::bad_alloc();
    }

    *(size_t *)block = size;

    g_currentAllocatedMemory += size;

    size_t currentPeak = g_peakAllocatedMemory.load();
    while (g_currentAllocatedMemory > currentPeak)
    {
        g_peakAllocatedMemory.compare_exchange_weak(currentPeak, g_currentAllocatedMemory.load());
    }

    return (char *)block + c_prefixSize;
}

void freeMemory(void *memory)
{
    if (!memory)
    {
        return;
    }

    void *block = (char *)memory - c_prefixSize;
    size_t size = *(size_t *)block;

    g_currentAllocatedMemory -= size;

    std::free(block);
}

void *operator new(size_t size)
{
    return allocateMemory(size);
}

void *operator new[](size_t size)
{
    return allocateMemory(size);
}

void operator delete(void *memory) noexcept
{
    freeMemory(memory);
}

void operator delete[](void *memory) noexcept
{
    freeMemory(memory);
}

void operator delete(void *memory, size_t size) noexcept
{
    (void)size;
    freeMemory(memory);
}

void operator delete[](void *memory, size_t size) noexcept
{
    (void)size;
    freeMemory(memory);
}

void resetPeakAllocation()
{
    g_peakAllocatedMemory = g_currentAllocatedMemory.load();
}

double getPeakAllocationKB()
{
    return g_peakAllocatedMemory / 1024.0;
}
