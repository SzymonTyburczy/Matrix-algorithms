#include "helperFunctions.h"

const double EPS = 1e-12;

Matrix get_submatrix(const Matrix &A, int r_start, int c_start, int num_rows, int num_cols);

std::vector<double> get_subvector(const std::vector<double> &b, int r_start, int num_rows);

void printVector(const std::vector<double> &v);