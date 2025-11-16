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

Matrix multiply_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count, MultiplyAlgorithm algo)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
    {
        throw std::runtime_error("Incompatible matrix dimensions for multiplication.");
    }
    if (A.size() == 0 || A[0].size() == 0 || B.size() == 0 || B[0].size() == 0)
    {
        return createMatrix(A.size(), B[0].size());
    }

    switch (algo)
    {
    case MultiplyAlgorithm::STRASSEN:
        return multiply_strassen_wrapper(A, B, op_count);
    case MultiplyAlgorithm::BINET:
        return multiply_recursive_wrapper(A, B, op_count);
    case MultiplyAlgorithm::ITERATIVE:
    default:
        return iterativeMultiply(A, B, op_count);
    }
}

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

    Matrix S_temp1 = multiply_wrapper(A21, A11_inv, op_count, algo);

    Matrix S_temp2 = multiply_wrapper(S_temp1, A12, op_count, algo);

    Matrix S = subtractMatrices(A22, S_temp2, op_count);

    Matrix S_inv = recursive_invert_internal(S, op_count, algo);
    Matrix B22 = S_inv;

    Matrix T1 = multiply_wrapper(A11_inv, A12, op_count, algo);

    Matrix B12_temp = multiply_wrapper(T1, S_inv, op_count, algo);
    Matrix B12 = subtractMatrices(createMatrix(n1, n2), B12_temp, op_count);

    Matrix B21_temp = multiply_wrapper(S_inv, S_temp1, op_count, algo);
    Matrix B21 = subtractMatrices(createMatrix(n2, n1), B21_temp, op_count);

    Matrix B11_temp = multiply_wrapper(T1, B21, op_count, algo);

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
        throw std::invalid_argument("Mtrix must be square and non-empty to be inverted.");
    }
    return recursive_invert_internal(A, op_count, algo);
}

int main()
{
    int N = 4;
    Matrix A = createMatrix(N, N, true);

    MultiplyAlgorithm algo_to_use = MultiplyAlgorithm::ITERATIVE;
    std::string algo_name = "Iterative";

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