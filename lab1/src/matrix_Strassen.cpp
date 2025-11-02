#include "helperFunctions.h"

// const int STRASSEN_LEAF_SIZE = 16;
const int STRASSEN_LEAF_SIZE = 0;

void strassenRecursive(Matrix &C, int rC, int cC,
                       const Matrix &A, int rA, int cA,
                       const Matrix &B, int rB, int cB,
                       int size, unsigned long long &op_count)
{

    if (size == 1)
    {
        C[rC][cC] = A[rA][cA] * B[rB][cB];
        op_count++;
        return;
    }

    if (size % 2 == 0)
    {
        int n_split = size / 2;

        Matrix s1 = subtractMatrices(B, rB, cB + n_split, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        Matrix s2 = addMatrices(A, rA, cA, A, rA, cA + n_split, n_split, n_split, op_count);
        Matrix s3 = addMatrices(A, rA + n_split, cA, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        Matrix s4 = subtractMatrices(B, rB + n_split, cB, B, rB, cB, n_split, n_split, op_count);
        Matrix s5 = addMatrices(A, rA, cA, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        Matrix s6 = addMatrices(B, rB, cB, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        Matrix s7 = subtractMatrices(A, rA, cA + n_split, A, rA + n_split, cA + n_split, n_split, n_split, op_count);
        Matrix s8 = addMatrices(B, rB + n_split, cB, B, rB + n_split, cB + n_split, n_split, n_split, op_count);
        Matrix s9 = subtractMatrices(A, rA, cA, A, rA + n_split, cA, n_split, n_split, op_count);
        Matrix s10 = addMatrices(B, rB, cB, B, rB, cB + n_split, n_split, n_split, op_count);

        Matrix p1 = createMatrix(n_split, n_split);
        Matrix p2 = createMatrix(n_split, n_split);
        Matrix p3 = createMatrix(n_split, n_split);
        Matrix p4 = createMatrix(n_split, n_split);
        Matrix p5 = createMatrix(n_split, n_split);
        Matrix p6 = createMatrix(n_split, n_split);
        Matrix p7 = createMatrix(n_split, n_split);

        strassenRecursive(p1, 0, 0, A, rA, cA, s1, 0, 0, n_split, op_count);
        strassenRecursive(p2, 0, 0, s2, 0, 0, B, rB + n_split, cB + n_split, n_split, op_count);
        strassenRecursive(p3, 0, 0, s3, 0, 0, B, rB, cB, n_split, op_count);
        strassenRecursive(p4, 0, 0, A, rA + n_split, cA + n_split, s4, 0, 0, n_split, op_count);
        strassenRecursive(p5, 0, 0, s5, 0, 0, s6, 0, 0, n_split, op_count);
        strassenRecursive(p6, 0, 0, s7, 0, 0, s8, 0, 0, n_split, op_count);
        strassenRecursive(p7, 0, 0, s9, 0, 0, s10, 0, 0, n_split, op_count);

        // C11 = P5 + P4 - P2 + P6
        addMatrices_inplace(C, rC, cC, p5, 0, 0, p4, 0, 0, n_split, n_split, op_count);
        subtractMatrices_inplace(C, rC, cC, C, rC, cC, p2, 0, 0, n_split, n_split, op_count);
        addMatrices_inplace(C, rC, cC, C, rC, cC, p6, 0, 0, n_split, n_split, op_count);

        // C12 = P1 + P2
        addMatrices_inplace(C, rC, cC + n_split, p1, 0, 0, p2, 0, 0, n_split, n_split, op_count);

        // C21 = P3 + P4
        addMatrices_inplace(C, rC + n_split, cC, p3, 0, 0, p4, 0, 0, n_split, n_split, op_count);

        // C22 = P5 + P1 - P3 - P7
        addMatrices_inplace(C, rC + n_split, cC + n_split, p5, 0, 0, p1, 0, 0, n_split, n_split, op_count);
        subtractMatrices_inplace(C, rC + n_split, cC + n_split, C, rC + n_split, cC + n_split, p3, 0, 0, n_split, n_split, op_count);
        subtractMatrices_inplace(C, rC + n_split, cC + n_split, C, rC + n_split, cC + n_split, p7, 0, 0, n_split, n_split, op_count);
    }
    else
    {
        int n1 = size - 1;

        // C11 = A11*B11 + a12*b21
        Matrix C11_p1 = createMatrix(n1, n1);
        strassenRecursive(C11_p1, 0, 0, A, rA, cA, B, rB, cB, n1, op_count);

        Matrix a12 = createMatrix(n1, 1);
        Matrix b21 = createMatrix(1, n1);
        for (int i = 0; i < n1; ++i)
            a12[i][0] = A[rA + i][cA + n1];
        for (int j = 0; j < n1; ++j)
            b21[0][j] = B[rB + n1][cB + j];

        Matrix C11_p2 = iterativeMultiply(a12, b21, op_count);
        addMatrices_inplace(C, rC, cC, C11_p1, 0, 0, C11_p2, 0, 0, n1, n1, op_count);

        // C12 = A11*b12 + a12*b22
        Matrix b12 = createMatrix(n1, 1);
        for (int i = 0; i < n1; ++i)
            b12[i][0] = B[rB + i][cB + n1];
        Matrix a22 = createMatrix(1, 1);
        a22[0][0] = A[rA + n1][cA + n1];
        Matrix b22 = createMatrix(1, 1);
        b22[0][0] = B[rB + n1][cB + n1];

        Matrix A11_sub = createMatrix(n1, n1);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n1; ++j)
                A11_sub[i][j] = A[rA + i][cA + j];

        Matrix C12_p1 = iterativeMultiply(A11_sub, b12, op_count);
        Matrix C12_p2 = iterativeMultiply(a12, b22, op_count);

        addMatrices_inplace(C, rC, cC + n1, C12_p1, 0, 0, C12_p2, 0, 0, n1, 1, op_count);

        // C21 = a21*B11 + a22*b21
        Matrix a21 = createMatrix(1, n1);
        for (int j = 0; j < n1; ++j)
            a21[0][j] = A[rA + n1][cA + j];

        Matrix B11_sub = createMatrix(n1, n1);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n1; ++j)
                B11_sub[i][j] = B[rB + i][cB + j];

        Matrix C21_p1 = iterativeMultiply(a21, B11_sub, op_count);
        Matrix C21_p2 = iterativeMultiply(a22, b21, op_count);

        addMatrices_inplace(C, rC + n1, cC, C21_p1, 0, 0, C21_p2, 0, 0, 1, n1, op_count);

        // C22 = a21*b12 + a22*b22
        Matrix C22_p1 = iterativeMultiply(a21, b12, op_count);
        Matrix C22_p2 = iterativeMultiply(a22, b22, op_count);

        C[rC + n1][cC + n1] = C22_p1[0][0] + C22_p2[0][0];
        op_count++;
    }
}

Matrix multiply_strassen_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }
    size_t n = A.size();
    if (n != A[0].size() || n != B.size() || n != B[0].size())
    {
        throw std::invalid_argument("Input matrices are not square and of the same size N.");
    }

    op_count = 0;
    Matrix C = createMatrix(n, n);

    strassenRecursive(C, 0, 0, A, 0, 0, B, 0, 0, n, op_count);

    return C;
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "Strassen's Matrix Multiplication" << std::endl;
    std::vector<int> vector_matrices = {2, 3, 5, 7, 20, 50, 100, 200};

    std::cout << " We consider following matrices sizes: " << std::endl;
    for (int n : vector_matrices)
    {
        printf("%d ", n);
    }
    std::cout << std::endl;

    unsigned long long general_op_count = 0;
    std::vector<resultsinCSV> results = {};

    double baselineMemoryKB = getPeakPrivateUsageKB();

    for (int n : vector_matrices)
    {
        Matrix A = createMatrix(n, n, true);
        Matrix B = createMatrix(n, n, true);

        auto start = std::chrono::high_resolution_clock::now();

        Matrix C = multiply_strassen_wrapper(A, B, general_op_count);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end - start;

        double peakMemoryKB = getPeakPrivateUsageKB();
        double algorithmMemoryKB = peakMemoryKB - baselineMemoryKB;

        std::cout << "Matrix size: " << n << "x" << n
                  << ", Operations count (Strassen): " << general_op_count
                  << ", Duration: " << duration.count() << " ms" << std::endl;

        printf("Memory Used (Algorithm): %.2f KB\n", algorithmMemoryKB);

        results.emplace_back(n, "Strassen", general_op_count, duration.count(), algorithmMemoryKB);

        std::cout << std::endl
                  << std::endl;
        general_op_count = 0;
    }

    std::ofstream csvfile("matrix_multiplication_results_STRASSEN.csv");
    csvfile << "Size,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.size << "," << res.algorithm << "," << res.operations << ","
                << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    return 0;
}