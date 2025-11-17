#include "helperFunctions.h"

// Version that works on VIEWS (offsets) of A and B, returns NEW matrix C
Matrix addMatrices(const Matrix &A, int rA, int cA,
                   const Matrix &B, int rB, int cB,
                   int rows, int cols, unsigned long long &op_count)
{
    Matrix C = createMatrix(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            C[i][j] = A[rA + i][cA + j] + B[rB + i][cB + j];
            op_count++;
        }
    }
    return C;
}

// Version that works on VIEWS (offsets) of A and B, returns NEW matrix C
Matrix subtractMatrices(const Matrix &A, int rA, int cA,
                        const Matrix &B, int rB, int cB,
                        int rows, int cols, unsigned long long &op_count)
{
    Matrix C = createMatrix(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            C[i][j] = A[rA + i][cA + j] - B[rB + i][cB + j];
            op_count++;
        }
    }
    return C;
}

// Version that works on VIEWS (offsets) of A and B, returns "in place" into C
void addMatrices_inplace(Matrix &C, int rC, int cC,
                         const Matrix &A, int rA, int cA,
                         const Matrix &B, int rB, int cB,
                         int rows, int cols, unsigned long long &op_count)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            C[rC + i][cC + j] = A[rA + i][cA + j] + B[rB + i][cB + j];
            op_count++;
        }
    }
}

// Version that works on VIEWS (offsets) of A and B, returns "in place" into C
void subtractMatrices_inplace(Matrix &C, int rC, int cC,
                              const Matrix &A, int rA, int cA,
                              const Matrix &B, int rB, int cB,
                              int rows, int cols, unsigned long long &op_count)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            C[rC + i][cC + j] = A[rA + i][cA + j] - B[rB + i][cB + j];
            op_count++;
        }
    }
}

// wrapper function for addMatrices (for whole matrices)
Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || A[0].empty())
        return B;
    if (B.empty() || B[0].empty())
        return A;
    size_t n = A.size();
    size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for addition.");
    }
    return addMatrices(A, 0, 0, B, 0, 0, n, m, op_count);
}

// wrapper function for subtractMatrices (for whole matrices)
Matrix subtractMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || A[0].empty())
    {
        return B;
    }
    if (B.empty() || B[0].empty())
        return A;
    size_t n = A.size();
    size_t m = A[0].size();
    if (n != B.size() || m != B[0].size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for subtraction.");
    }
    return subtractMatrices(A, 0, 0, B, 0, 0, n, m, op_count);
}

Matrix createMatrix(size_t rows, size_t cols, bool random)
{
    Matrix mat(rows, std::vector<double>(cols, 0.0));
    if (random && rows > 0 && cols > 0)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(0.00000001, 1.0);

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                mat[i][j] = dis(gen);
            }
        }
    }
    return mat;
}

Matrix iterativeMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    size_t m = A.size();
    size_t k = (m > 0) ? A[0].size() : 0;
    size_t p = (B.size() > 0) ? B[0].size() : 0;

    if (k != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    Matrix C = createMatrix(m, p);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < p; ++j)
        {
            double sum = 0.0;
            for (size_t l = 0; l < k; ++l)
            {
                sum += A[i][l] * B[l][j];
                op_count++;
            }
            C[i][j] = sum;
            if (k > 1)
            {
                op_count += (k - 1);
            }
        }
    }
    return C;
}

void printMatrix(const std::vector<std::vector<double>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (const auto &val : row)
        {
            printf("%0.4f ", val);
        }
        std::cout << std::endl;
    }
}

double printMemoryUsage()
{
    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        double peakMemoryKB = pmc.PeakWorkingSetSize / 1024.0;
        std::cout << "Memory Used: " << peakMemoryKB << " KB" << std::endl;
        return peakMemoryKB;
    }
    return 0;
}

double getPeakPrivateUsageKB()
{
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        return pmc.PeakPagefileUsage / 1024.0;
    }
    return 0;
}

void iterativeMultiply_inplace(Matrix &C, int rC, int cC,
                               const Matrix &A, int rA, int cA,
                               const Matrix &B, int rB, int cB,
                               int m, int k, int p, unsigned long long &op_count)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            double sum = 0.0;
            for (int l = 0; l < k; ++l)
            {
                sum += A[rA + i][cA + l] * B[rB + l][cB + j];
                op_count += 2;
            }
            C[rC + i][cC + j] = sum;
            if (k > 1)
            {
                op_count += (k - 1);
            }
        }
    }
}