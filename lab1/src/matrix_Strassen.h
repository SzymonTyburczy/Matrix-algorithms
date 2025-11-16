#pragma once
#include "helperFunctions.h"
Matrix multiply_strassen_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count);
void multiply_strassen_inplace(Matrix &C, int rC, int cC,
                               const Matrix &A, int rA, int cA,
                               const Matrix &B, int rB, int cB,
                               int m, int k, int p, unsigned long long &op_count);