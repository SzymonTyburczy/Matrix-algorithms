#pragma once
#include "helperFunctions.h"
#include <vector>

extern const double EPS;

enum class MultiplyAlgorithm
{
    ITERATIVE,
    BINET,
    STRASSEN
};

Matrix get_submatrix(const Matrix &A, int r_start, int c_start, int num_rows, int num_cols);

std::vector<double> get_subvector(const std::vector<double> &b, int r_start, int num_rows);

void printVector(const std::vector<double> &v);

std::vector<double> solve_lower_triangular(const Matrix &L, const std::vector<double> &b, unsigned long long &op_count);
std::vector<double> solve_upper_triangular(const Matrix &U, const std::vector<double> &b, unsigned long long &op_count);

void copy_block_to_matrix(Matrix &Target, const Matrix &Source,
                          int r_target, int c_target);