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

Matrix recursive_invert_internal(const Matrix &A, unsigned long long &op_count, MultiplyAlgorithm algo)
{
    int n = A.size();
    if (n == 0)
    {
        return createMatrix(0, 0);
    }

    if (n == 1)
    {
        if (fabs(A[0][0]) < EPS)
        {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }
        op_count++;
        return {{1.0 / A[0][0]}};
    }

    int n1 = n / 2;
    int n2 = n - n1;

    Matrix A11 = get_submatrix(A, 0, 0, n1, n1);
    Matrix A12 = get_submatrix(A, 0, n1, n1, n2);
    Matrix A21 = get_submatrix(A, n1, 0, n2, n1);
    Matrix A22 = get_submatrix(A, n1, n1, n2, n2);

    Matrix A11_inv = recursive_invert_internal(A11, op_count, algo);

    Matrix S_temp1 = createMatrix(n2, n1);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(S_temp1, 0, 0, A21, 0, 0, A11_inv, 0, 0, n2, n1, n1, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(S_temp1, 0, 0, A21, 0, 0, A11_inv, 0, 0, n2, n1, n1, op_count);
        break;
    default:
        iterativeMultiply_inplace(S_temp1, 0, 0, A21, 0, 0, A11_inv, 0, 0, n2, n1, n1, op_count);
        break;
    }

    Matrix S_temp2 = createMatrix(n2, n2);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(S_temp2, 0, 0, S_temp1, 0, 0, A12, 0, 0, n2, n1, n2, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(S_temp2, 0, 0, S_temp1, 0, 0, A12, 0, 0, n2, n1, n2, op_count);
        break;
    default:
        iterativeMultiply_inplace(S_temp2, 0, 0, S_temp1, 0, 0, A12, 0, 0, n2, n1, n2, op_count);
        break;
    }

    Matrix S = subtractMatrices(A22, S_temp2, op_count);

    Matrix S_inv = recursive_invert_internal(S, op_count, algo);
    Matrix B22 = S_inv;

    Matrix T1 = createMatrix(n1, n2);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(T1, 0, 0, A11_inv, 0, 0, A12, 0, 0, n1, n1, n2, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(T1, 0, 0, A11_inv, 0, 0, A12, 0, 0, n1, n1, n2, op_count);
        break;
    default:
        iterativeMultiply_inplace(T1, 0, 0, A11_inv, 0, 0, A12, 0, 0, n1, n1, n2, op_count);
        break;
    }

    Matrix B12_temp = createMatrix(n1, n2);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(B12_temp, 0, 0, T1, 0, 0, S_inv, 0, 0, n1, n2, n2, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(B12_temp, 0, 0, T1, 0, 0, S_inv, 0, 0, n1, n2, n2, op_count);
        break;
    default:
        iterativeMultiply_inplace(B12_temp, 0, 0, T1, 0, 0, S_inv, 0, 0, n1, n2, n2, op_count);
        break;
    }
    Matrix B12 = subtractMatrices(createMatrix(n1, n2), B12_temp, op_count);

    Matrix B21_temp = createMatrix(n2, n1);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(B21_temp, 0, 0, S_inv, 0, 0, S_temp1, 0, 0, n2, n2, n1, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(B21_temp, 0, 0, S_inv, 0, 0, S_temp1, 0, 0, n2, n2, n1, op_count);
        break;
    default:
        iterativeMultiply_inplace(B21_temp, 0, 0, S_inv, 0, 0, S_temp1, 0, 0, n2, n2, n1, op_count);
        break;
    }
    Matrix B21 = subtractMatrices(createMatrix(n2, n1), B21_temp, op_count);

    Matrix B11_temp = createMatrix(n1, n1);
    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        multiply_strassen_inplace(B11_temp, 0, 0, T1, 0, 0, B21, 0, 0, n1, n2, n1, op_count);
        break;
    case MultiplyAlgorithm::BINET:
        multiply_binet_inplace(B11_temp, 0, 0, T1, 0, 0, B21, 0, 0, n1, n2, n1, op_count);
        break;
    default:
        iterativeMultiply_inplace(B11_temp, 0, 0, T1, 0, 0, B21, 0, 0, n1, n2, n1, op_count);
        break;
    }

    Matrix B11 = subtractMatrices(A11_inv, B11_temp, op_count);

    Matrix B = createMatrix(n, n);
    copy_block_to_matrix(B, B11, 0, 0);
    copy_block_to_matrix(B, B12, 0, n1);
    copy_block_to_matrix(B, B21, n1, 0);
    copy_block_to_matrix(B, B22, n1, n1);

    return B;
}

Matrix recursive_invert(const Matrix &A, unsigned long long &op_count, MultiplyAlgorithm algo)
{
    if (A.size() != A[0].size() || A.empty())
    {
        throw std::invalid_argument("Matrix must be square and non-empty to be inverted.");
    }
    return recursive_invert_internal(A, op_count, algo);
}

/*
int main()
{
    int N = 4;
    Matrix A = createMatrix(N, N, true);

    MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::ITERATIVE;
    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::BINET;
    // MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::STRASSEN;

    std::string algo_name;
    switch (algo_to_use)
    {
    case MultiplyAlgorithm::STRASSEN:
        algo_name = "Strassen";
        break;
    case MultiplyAlgorithm::BINET:
        algo_name = "Binet";
        break;
    default:
        algo_name = "Iterative";
        break;
    }

    std::cout << "Original Matrix A:" << std::endl;
    printMatrix(A);
    std::cout << "\nUsing multiplication algorithm: " << algo_name << "\n"
              << std::endl;

    unsigned long long flops = 0;
    try
    {
        Matrix A_inv = recursive_invert(A, flops, algo_to_use);

        std::cout << "Inverse Matrix A_inv:" << std::endl;
        printMatrix(A_inv);
        std::cout << std::endl;

        unsigned long long check_flops = 0;
        Matrix I = iterativeMultiply(A, A_inv, check_flops);

        std::cout << "Check: A * A_inv (should be close to I):" << std::endl;
        printMatrix(I);
        std::cout << std::endl;
        std::cout << "Total number of operations (FLOPS): " << flops << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error during matrix inversion: " << e.what() << std::endl;
    }

    return 0;
}
*/