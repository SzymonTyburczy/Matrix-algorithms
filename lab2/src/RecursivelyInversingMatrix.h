#pragma once
#include "helperFunctions.h"
#include "RecursiveLUFactorization.h"

Matrix recursive_invert(const Matrix &A,
                        unsigned long long &op_count,
                        MultiplyAlgorithm algo);