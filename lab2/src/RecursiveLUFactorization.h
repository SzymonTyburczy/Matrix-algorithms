#pragma once
#include "helperFunctions.h"
#include <string>
#include "HelperFunctionsLab2.h"

struct LU_Result
{
    Matrix L;
    Matrix U;
    double determinant;
};

LU_Result recursive_lu_factorization(
    const Matrix &A,
    unsigned long long &op_count,
    MultiplyAlgorithm algo = MultiplyAlgorithm::ITERATIVE);