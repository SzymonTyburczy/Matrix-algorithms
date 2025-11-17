#pragma once
#include <vector>
#include "helperFunctions.h"
#include "HelperFunctionsLab2.h"
#include "RecursiveLUFactorization.h"

std::vector<double> solve_block_gauss(Matrix A,
                                      std::vector<double> b,
                                      unsigned long long &flop_count,
                                      int block_size,
                                      MultiplyAlgorithm algo);