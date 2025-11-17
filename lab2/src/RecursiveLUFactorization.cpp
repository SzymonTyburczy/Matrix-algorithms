#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <string>

#include "helperFunctions.h"
#include "HelperFunctionsLab2.h"
#include "matrix_Strassen.h"
#include "matrix_Binet.h"
#include "RecursiveLUFactorization.h"

Matrix solve_left_unit_lower_triangular(const Matrix &L, const Matrix &B, unsigned long long &op_count)
{
    int n = L.size();
    int m = B[0].size();
    Matrix X = createMatrix(n, m);
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            for (int k = 0; k < i; ++k)
            {
                sum += L[i][k] * X[k][j];
                op_count += 2;
            }
            X[i][j] = B[i][j] - sum;
            op_count++;
        }
    }
    return X;
}

Matrix solve_right_upper_triangular(const Matrix &A, const Matrix &U, unsigned long long &op_count)
{
    int n = U.size();
    int m = A.size();
    Matrix X = createMatrix(m, n);

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < j; ++k)
            {
                sum += X[i][k] * U[k][j];
                op_count += 2;
            }

            if (fabs(U[j][j]) < EPS)
            {
                throw std::runtime_error("Zero diagonal element in solve_right_upper_triangular.");
            }

            X[i][j] = (A[i][j] - sum) / U[j][j];
            op_count += 2;
        }
    }
    return X;
}

LU_Result recursive_lu_factorization(const Matrix &A, unsigned long long &op_count,
                                     MultiplyAlgorithm algo)
{
    int n = A.size();
    if (n == 0)
    {
        return {createMatrix(0, 0), createMatrix(0, 0), 1.0};
    }
    if (n == 1)
    {
        if (fabs(A[0][0]) < EPS)
        {
            throw std::runtime_error("Pivot near zero in recursive_lu_factorization base case.");
        }
        Matrix L = {{1.0}};
        Matrix U = {{A[0][0]}};
        return {L, U, A[0][0]};
    }

    int n1 = n / 2;
    int n2 = n - n1;

    Matrix A11 = get_submatrix(A, 0, 0, n1, n1);
    Matrix A12 = get_submatrix(A, 0, n1, n1, n2);
    Matrix A21 = get_submatrix(A, n1, 0, n2, n1);
    Matrix A22 = get_submatrix(A, n1, n1, n2, n2);

    LU_Result res11 = recursive_lu_factorization(A11, op_count, algo);
    Matrix &L11 = res11.L;
    Matrix &U11 = res11.U;
    double det1 = res11.determinant;

    Matrix U12 = solve_left_unit_lower_triangular(L11, A12, op_count);
    Matrix L21 = solve_right_upper_triangular(A21, U11, op_count);

    Matrix S_temp = createMatrix(n2, n2);

    int m = L21.size();
    int k = L21[0].size();
    int p = U12[0].size();

    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(S_temp, 0, 0, L21, 0, 0, U12, 0, 0, m, k, p, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(S_temp, 0, 0, L21, 0, 0, U12, 0, 0, m, k, p, op_count);
        break;
    case MultiplyAlgorithm::ITERATIVE:
    default:
        iterativeMultiply_inplace(S_temp, 0, 0, L21, 0, 0, U12, 0, 0, m, k, p, op_count);
        break;
    }

    Matrix S = subtractMatrices(A22, S_temp, op_count);

    LU_Result res22 = recursive_lu_factorization(S, op_count, algo);
    Matrix &L22 = res22.L;
    Matrix &U22 = res22.U;
    double det2 = res22.determinant;

    double total_determinant = det1 * det2;

    Matrix L = createMatrix(n, n, false);
    Matrix U = createMatrix(n, n, false);

    copy_block_to_matrix(L, L11, 0, 0);
    copy_block_to_matrix(L, L21, n1, 0);
    copy_block_to_matrix(L, L22, n1, n1);

    copy_block_to_matrix(U, U11, 0, 0);
    copy_block_to_matrix(U, U12, 0, n1);
    copy_block_to_matrix(U, U22, n1, n1);

    return {L, U, total_determinant};
}

/*
int main()
{
    int N = 4;
    Matrix A = createMatrix(N, N, true);

    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::ITERATIVE;
    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::BINET;
    MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::STRASSEN;

    std::string algo_name;
    switch (algo_to_use)
    {
    case MultiplyAlgorithm::STRASSEN:
        algo_name = "Strassen";
        break;
    case MultiplyAlgorithm::BINET:
        algo_name = "Binet (Recursive)";
        break;
    default:
        algo_name = "Iterative";
        break;
    }

    std::cout << "Original Matrix A:" << std::endl;
    printMatrix(A);
    std::cout << "\nUsing Multiplication Algorithm: " << algo_name << "\n"
              << std::endl;

    unsigned long long flops = 0;
    try
    {
        LU_Result result = recursive_lu_factorization(A, flops, algo_to_use);

        std::cout << "Matrix L (lower triangular, 1 on diagonal):" << std::endl;
        printMatrix(result.L);
        std::cout << std::endl;

        std::cout << "Matrix U (upper triangular):" << std::endl;
        printMatrix(result.U);
        std::cout << std::endl;

        unsigned long long check_flops = 0;
        Matrix LU = iterativeMultiply(result.L, result.U, check_flops);
        std::cout << "Check: L * U =" << std::endl;
        printMatrix(LU);
        std::cout << std::endl;

        std::cout << "Determinant det(A) = " << result.determinant << std::endl;
        std::cout << "Total number of operations (FLOPS): " << flops << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "LU factorization error: " << e.what() << std::endl;
    }

    return 0;
}
*/