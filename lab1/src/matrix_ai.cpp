#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <windows.h>
#include <psapi.h>
#include <iomanip>

#include "helperFunctions.h"

using Matrix = std::vector<std::vector<double>>;

// Structure to hold AI results
struct AIRresult
{
    int n_level; // 0 for 4x5, 1 for 8x10, etc.
    std::string dimensions;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;
};

void matrix_ai_inplace(Matrix &C, int rC, int cC,
                       const Matrix &A, int rA, int cA,
                       const Matrix &B, int rB, int cB,
                       unsigned long long &op_count)
{

    std::vector<double> H(76);

    H[0] = A[rA + 2][cA + 1] * (-B[rB + 1][cB + 0] - B[rB + 1][cB + 4] - B[rB + 2][cB + 0]);
    op_count += 3;
    H[1] = (A[rA + 1][cA + 1] + A[rA + 1][cA + 4] - A[rA + 2][cA + 4]) * (-B[rB + 1][cB + 4] - B[rB + 4][cB + 0]);
    op_count += 4;
    H[2] = (-A[rA + 2][cA + 0] - A[rA + 3][cA + 0] + A[rA + 3][cA + 1]) * (-B[rB + 0][cB + 0] + B[rB + 1][cB + 4]);
    op_count += 4;
    H[3] = (A[rA + 0][cA + 1] + A[rA + 0][cA + 3] + A[rA + 2][cA + 3]) * (-B[rB + 1][cB + 4] - B[rB + 3][cB + 0]);
    op_count += 4;
    H[4] = (A[rA + 0][cA + 4] + A[rA + 1][cA + 1] + A[rA + 1][cA + 4]) * (-B[rB + 1][cB + 3] + B[rB + 4][cB + 0]);
    op_count += 4;
    H[5] = (-A[rA + 1][cA + 1] - A[rA + 1][cA + 4] - A[rA + 3][cA + 4]) * (B[rB + 1][cB + 2] + B[rB + 4][cB + 0]);
    op_count += 4;
    H[6] = (-A[rA + 0][cA + 0] + A[rA + 3][cA + 0] - A[rA + 3][cA + 1]) * (B[rB + 0][cB + 0] + B[rB + 1][cB + 3]);
    op_count += 4;
    H[7] = (A[rA + 2][cA + 1] - A[rA + 2][cA + 2] - A[rA + 3][cA + 2]) * (-B[rB + 1][cB + 2] + B[rB + 2][cB + 0]);
    op_count += 4;
    H[8] = (-A[rA + 0][cA + 1] - A[rA + 0][cA + 3] + A[rA + 3][cA + 3]) * (B[rB + 1][cB + 2] + B[rB + 3][cB + 0]);
    op_count += 4;
    H[9] = (A[rA + 1][cA + 1] + A[rA + 1][cA + 4]) * B[rB + 4][cB + 0];
    op_count += 2;
    H[10] = (-A[rA + 1][cA + 0] - A[rA + 3][cA + 0] + A[rA + 3][cA + 1]) * (-B[rB + 0][cB + 0] + B[rB + 1][cB + 1]);
    op_count += 4;
    H[11] = (A[rA + 3][cA + 0] - A[rA + 3][cA + 1]) * B[rB + 0][cB + 0];
    op_count += 2;
    H[12] = (A[rA + 0][cA + 1] + A[rA + 0][cA + 3] + A[rA + 1][cA + 3]) * (B[rB + 1][cB + 1] + B[rB + 3][cB + 0]);
    op_count += 4;
    H[13] = (A[rA + 0][cA + 2] - A[rA + 2][cA + 1] + A[rA + 2][cA + 2]) * (B[rB + 1][cB + 3] + B[rB + 2][cB + 0]);
    op_count += 4;
    H[14] = (-A[rA + 0][cA + 1] - A[rA + 0][cA + 3]) * B[rB + 3][cB + 0];
    op_count += 2;
    H[15] = (-A[rA + 2][cA + 1] + A[rA + 2][cA + 2]) * B[rB + 2][cB + 0];
    op_count += 2;
    H[16] = (A[rA + 0][cA + 1] + A[rA + 0][cA + 3] - A[rA + 1][cA + 0] + A[rA + 1][cA + 1] - A[rA + 1][cA + 2] + A[rA + 1][cA + 3] - A[rA + 2][cA + 1] + A[rA + 2][cA + 2] - A[rA + 3][cA + 0] + A[rA + 3][cA + 1]) * B[rB + 1][cB + 1];
    op_count += 10;
    H[17] = A[rA + 1][cA + 0] * (B[rB + 0][cB + 0] + B[rB + 0][cB + 1] + B[rB + 4][cB + 1]);
    op_count += 3;
    H[18] = -A[rA + 1][cA + 2] * (B[rB + 2][cB + 0] + B[rB + 2][cB + 1] + B[rB + 4][cB + 1]);
    op_count += 3;
    H[19] = (-A[rA + 0][cA + 4] + A[rA + 1][cA + 0] + A[rA + 1][cA + 2] - A[rA + 1][cA + 4]) * (-B[rB + 0][cB + 0] - B[rB + 0][cB + 1] + B[rB + 0][cB + 3] - B[rB + 4][cB + 1]);
    op_count += 8;
    H[20] = (A[rA + 1][cA + 0] + A[rA + 1][cA + 2] - A[rA + 1][cA + 4]) * B[rB + 4][cB + 1];
    op_count += 3;
    H[21] = (A[rA + 0][cA + 2] - A[rA + 0][cA + 3] - A[rA + 1][cA + 3]) * (B[rB + 0][cB + 0] + B[rB + 0][cB + 1] - B[rB + 0][cB + 3] - B[rB + 2][cB + 0] - B[rB + 2][cB + 1] + B[rB + 2][cB + 3] + B[rB + 3][cB + 3]);
    op_count += 9;
    H[22] = A[rA + 0][cA + 2] * (-B[rB + 2][cB + 0] + B[rB + 2][cB + 3] + B[rB + 3][cB + 3]);
    op_count += 3;
    H[23] = A[rA + 0][cA + 4] * (-B[rB + 3][cB + 3] - B[rB + 4][cB + 0] + B[rB + 4][cB + 3]);
    op_count += 3;
    H[24] = -A[rA + 0][cA + 0] * (B[rB + 0][cB + 0] - B[rB + 0][cB + 3]);
    op_count += 2;
    H[25] = (-A[rA + 0][cA + 2] + A[rA + 0][cA + 3] + A[rA + 0][cA + 4]) * B[rB + 3][cB + 3];
    op_count += 3;
    H[26] = (A[rA + 0][cA + 2] - A[rA + 2][cA + 0] + A[rA + 2][cA + 2]) * (B[rB + 0][cB + 0] - B[rB + 0][cB + 3] + B[rB + 0][cB + 4] + B[rB + 2][cB + 4]);
    op_count += 6;
    H[27] = -A[rA + 2][cA + 3] * (-B[rB + 2][cB + 4] - B[rB + 3][cB + 0] - B[rB + 3][cB + 4]);
    op_count += 3;
    H[28] = A[rA + 2][cA + 0] * (B[rB + 0][cB + 0] + B[rB + 0][cB + 4] + B[rB + 2][cB + 4]);
    op_count += 3;
    H[29] = (A[rA + 2][cA + 0] - A[rA + 2][cA + 2] + A[rA + 2][cA + 3]) * B[rB + 2][cB + 4];
    op_count += 3;
    H[30] = (-A[rA + 0][cA + 3] - A[rA + 0][cA + 4] - A[rA + 2][cA + 3]) * (-B[rB + 3][cB + 3] - B[rB + 4][cB + 0] + B[rB + 4][cB + 3] - B[rB + 4][cB + 4]);
    op_count += 7;
    H[31] = (A[rA + 1][cA + 0] + A[rA + 3][cA + 0] + A[rA + 3][cA + 3]) * (B[rB + 0][cB + 2] - B[rB + 3][cB + 0] - B[rB + 3][cB + 1] - B[rB + 3][cB + 2]);
    op_count += 6;
    H[32] = A[rA + 3][cA + 2] * (-B[rB + 2][cB + 0] - B[rB + 2][cB + 2]);
    op_count += 2;
    H[33] = A[rA + 3][cA + 3] * (-B[rB + 0][cB + 2] + B[rB + 3][cB + 0] + B[rB + 3][cB + 2]);
    op_count += 3;
    H[34] = -A[rA + 3][cA + 4] * (B[rB + 0][cB + 2] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2]);
    op_count += 3;
    H[35] = (A[rA + 1][cA + 2] - A[rA + 1][cA + 4] - A[rA + 3][cA + 4]) * (B[rB + 2][cB + 0] + B[rB + 2][cB + 1] + B[rB + 2][cB + 2] + B[rB + 4][cB + 1]);
    op_count += 6;
    H[36] = (-A[rA + 3][cA + 0] - A[rA + 3][cA + 3] + A[rA + 3][cA + 4]) * B[rB + 0][cB + 2];
    op_count += 3;
    H[37] = (-A[rA + 1][cA + 2] - A[rA + 2][cA + 0] + A[rA + 2][cA + 2] - A[rA + 2][cA + 3]) * (B[rB + 2][cB + 4] + B[rB + 3][cB + 0] + B[rB + 3][cB + 1] + B[rB + 3][cB + 4]);
    op_count += 7;
    H[38] = (-A[rA + 2][cA + 0] - A[rA + 3][cA + 0] - A[rA + 3][cA + 3] + A[rA + 3][cA + 4]) * (B[rB + 0][cB + 2] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2] + B[rB + 4][cB + 4]);
    op_count += 7;
    H[39] = (-A[rA + 0][cA + 2] + A[rA + 0][cA + 3] + A[rA + 0][cA + 4] - A[rA + 3][cA + 3]) * (-B[rB + 2][cB + 0] - B[rB + 2][cB + 2] + B[rB + 2][cB + 3] + B[rB + 3][cB + 3]);
    op_count += 7;
    H[40] = (-A[rA + 0][cA + 0] + A[rA + 3][cA + 0] - A[rA + 3][cA + 4]) * (B[rB + 0][cB + 2] + B[rB + 2][cB + 0] + B[rB + 2][cB + 2] - B[rB + 2][cB + 3] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2] - B[rB + 4][cB + 3]);
    op_count += 9;
    H[41] = (-A[rA + 1][cA + 0] + A[rA + 1][cA + 4] - A[rA + 2][cA + 4]) * (-B[rB + 0][cB + 0] - B[rB + 0][cB + 1] - B[rB + 0][cB + 4] + B[rB + 3][cB + 0] + B[rB + 3][cB + 1] + B[rB + 3][cB + 4] - B[rB + 4][cB + 1]);
    op_count += 10;
    H[42] = A[rA + 1][cA + 3] * (B[rB + 3][cB + 0] + B[rB + 3][cB + 1]);
    op_count += 2;
    H[43] = (A[rA + 1][cA + 2] + A[rA + 2][cA + 1] - A[rA + 2][cA + 2]) * (B[rB + 1][cB + 1] - B[rB + 2][cB + 0]);
    op_count += 4;
    H[44] = (-A[rA + 2][cA + 2] + A[rA + 2][cA + 3] - A[rA + 3][cA + 2]) * (B[rB + 2][cB + 4] + B[rB + 3][cB + 0] + B[rB + 3][cB + 2] + B[rB + 3][cB + 4] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2] + B[rB + 4][cB + 4]);
    op_count += 9;
    H[45] = -A[rA + 2][cA + 4] * (-B[rB + 4][cB + 0] - B[rB + 4][cB + 4]);
    op_count += 2;
    H[46] = (A[rA + 1][cA + 0] - A[rA + 1][cA + 4] - A[rA + 2][cA + 0] + A[rA + 2][cA + 4]) * (B[rB + 0][cB + 0] + B[rB + 0][cB + 1] + B[rB + 0][cB + 4] - B[rB + 3][cB + 0] - B[rB + 3][cB + 1] - B[rB + 3][cB + 4]);
    op_count += 9;
    H[47] = (-A[rA + 1][cA + 2] + A[rA + 2][cA + 2]) * (B[rB + 1][cB + 1] + B[rB + 2][cB + 1] + B[rB + 2][cB + 4] + B[rB + 3][cB + 0] + B[rB + 3][cB + 1] + B[rB + 3][cB + 4]);
    op_count += 7;
    H[48] = (-A[rA + 0][cA + 0] - A[rA + 0][cA + 2] + A[rA + 0][cA + 3] + A[rA + 0][cA + 4] - A[rA + 1][cA + 0] - A[rA + 1][cA + 2] + A[rA + 1][cA + 3] + A[rA + 1][cA + 4]) * (-B[rB + 0][cB + 0] - B[rB + 0][cB + 1] + B[rB + 0][cB + 3]);
    op_count += 11;
    H[49] = (-A[rA + 0][cA + 3] - A[rA + 1][cA + 3]) * (B[rB + 1][cB + 1] - B[rB + 2][cB + 0] - B[rB + 2][cB + 1] + B[rB + 2][cB + 3] - B[rB + 3][cB + 1] + B[rB + 3][cB + 3]);
    op_count += 7;
    H[50] = A[rA + 1][cA + 1] * (B[rB + 1][cB + 0] + B[rB + 1][cB + 1] - B[rB + 4][cB + 0]);
    op_count += 3;
    H[51] = A[rA + 3][cA + 1] * (B[rB + 0][cB + 0] + B[rB + 1][cB + 0] + B[rB + 1][cB + 2]);
    op_count += 3;
    H[52] = -A[rA + 0][cA + 1] * (-B[rB + 1][cB + 0] + B[rB + 1][cB + 3] + B[rB + 3][cB + 0]);
    op_count += 3;
    H[53] = (A[rA + 0][cA + 1] + A[rA + 0][cA + 3] - A[rA + 1][cA + 1] - A[rA + 1][cA + 4] - A[rA + 2][cA + 1] + A[rA + 2][cA + 2] - A[rA + 3][cA + 1] + A[rA + 3][cA + 2] - A[rA + 3][cA + 3] - A[rA + 3][cA + 4]) * B[rB + 1][cB + 2];
    op_count += 10;
    H[54] = (A[rA + 0][cA + 3] - A[rA + 3][cA + 3]) * (-B[rB + 1][cB + 2] + B[rB + 2][cB + 0] + B[rB + 2][cB + 2] - B[rB + 2][cB + 3] + B[rB + 3][cB + 2] - B[rB + 3][cB + 3]);
    op_count += 7;
    H[55] = (A[rA + 0][cA + 0] - A[rA + 0][cA + 4] - A[rA + 3][cA + 0] + A[rA + 3][cA + 4]) * (B[rB + 2][cB + 0] + B[rB + 2][cB + 2] - B[rB + 2][cB + 3] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2] - B[rB + 4][cB + 3]);
    op_count += 9;
    H[56] = (-A[rA + 2][cA + 0] - A[rA + 3][cA + 0]) * (-B[rB + 0][cB + 2] - B[rB + 0][cB + 4] - B[rB + 1][cB + 4] - B[rB + 4][cB + 0] - B[rB + 4][cB + 2] - B[rB + 4][cB + 4]);
    op_count += 7;
    H[57] = (-A[rA + 0][cA + 3] - A[rA + 0][cA + 4] - A[rA + 2][cA + 3] - A[rA + 2][cA + 4]) * (-B[rB + 4][cB + 0] + B[rB + 4][cB + 3] - B[rB + 4][cB + 4]);
    op_count += 6;
    H[58] = (-A[rA + 2][cA + 2] + A[rA + 2][cA + 3] - A[rA + 3][cA + 2] + A[rA + 3][cA + 3]) * (B[rB + 3][cB + 0] + B[rB + 3][cB + 2] + B[rB + 3][cB + 4] + B[rB + 4][cB + 0] + B[rB + 4][cB + 2] + B[rB + 4][cB + 4]);
    op_count += 9;
    H[59] = (A[rA + 1][cA + 4] + A[rA + 3][cA + 4]) * (B[rB + 1][cB + 2] - B[rB + 2][cB + 0] - B[rB + 2][cB + 1] - B[rB + 2][cB + 2] - B[rB + 4][cB + 1] - B[rB + 4][cB + 2]);
    op_count += 7;
    H[60] = (A[rA + 0][cA + 3] + A[rA + 2][cA + 3]) * (B[rB + 0][cB + 0] - B[rB + 0][cB + 3] + B[rB + 0][cB + 4] - B[rB + 1][cB + 4] - B[rB + 3][cB + 3] + B[rB + 3][cB + 4] - B[rB + 4][cB + 0] + B[rB + 4][cB + 3] - B[rB + 4][cB + 4]);
    op_count += 10;
    H[61] = (A[rA + 1][cA + 0] + A[rA + 3][cA + 0]) * (B[rB + 0][cB + 1] + B[rB + 0][cB + 2] + B[rB + 1][cB + 1] - B[rB + 3][cB + 0] - B[rB + 3][cB + 1] - B[rB + 3][cB + 2]);
    op_count += 7;
    H[62] = (-A[rA + 2][cA + 2] - A[rA + 3][cA + 2]) * (-B[rB + 1][cB + 2] - B[rB + 2][cB + 2] - B[rB + 2][cB + 4] - B[rB + 3][cB + 0] - B[rB + 3][cB + 2] - B[rB + 3][cB + 4]);
    op_count += 7;
    H[63] = (A[rA + 0][cA + 0] - A[rA + 0][cA + 2] - A[rA + 0][cA + 3] + A[rA + 2][cA + 0] - A[rA + 2][cA + 2] - A[rA + 2][cA + 3]) * (B[rB + 0][cB + 0] - B[rB + 0][cB + 3] + B[rB + 0][cB + 4]);
    op_count += 8;
    H[64] = (-A[rA + 0][cA + 0] + A[rA + 3][cA + 0]) * (-B[rB + 0][cB + 2] + B[rB + 0][cB + 3] + B[rB + 1][cB + 3] - B[rB + 4][cB + 0] - B[rB + 4][cB + 2] + B[rB + 4][cB + 3]);
    op_count += 7;
    H[65] = (A[rA + 0][cA + 0] - A[rA + 0][cA + 1] + A[rA + 0][cA + 2] - A[rA + 0][cA + 4] - A[rA + 1][cA + 1] - A[rA + 1][cA + 4] - A[rA + 2][cA + 1] + A[rA + 2][cA + 2] - A[rA + 3][cA + 0] + A[rA + 3][cA + 1]) * B[rB + 1][cB + 3];
    op_count += 10;
    H[66] = (A[rA + 1][cA + 4] - A[rA + 2][cA + 4]) * (B[rB + 0][cB + 0] + B[rB + 0][cB + 1] + B[rB + 0][cB + 4] - B[rB + 1][cB + 4] - B[rB + 3][cB + 0] - B[rB + 3][cB + 1] - B[rB + 3][cB + 4] + B[rB + 4][cB + 1] + B[rB + 4][cB + 4]);
    op_count += 10;
    H[67] = (A[rA + 0][cA + 0] + A[rA + 0][cA + 2] - A[rA + 0][cA + 3] - A[rA + 0][cA + 4] - A[rA + 3][cA + 0] - A[rA + 3][cA + 2] + A[rA + 3][cA + 3] + A[rA + 3][cA + 4]) * (-B[rB + 2][cB + 0] - B[rB + 2][cB + 2] + B[rB + 2][cB + 3]);
    op_count += 10;
    H[68] = (-A[rA + 0][cA + 2] + A[rA + 0][cA + 3] - A[rA + 1][cA + 2] + A[rA + 1][cA + 3]) * (-B[rB + 1][cB + 3] - B[rB + 2][cB + 0] - B[rB + 2][cB + 1] + B[rB + 2][cB + 3] - B[rB + 4][cB + 1] + B[rB + 4][cB + 3]);
    op_count += 9;
    H[69] = (A[rA + 1][cA + 2] - A[rA + 1][cA + 4] + A[rA + 3][cA + 2] - A[rA + 3][cA + 4]) * (-B[rB + 2][cB + 0] - B[rB + 2][cB + 1] - B[rB + 2][cB + 2]);
    op_count += 6;
    H[70] = (-A[rA + 2][cA + 0] + A[rA + 2][cA + 2] - A[rA + 2][cA + 3] + A[rA + 2][cA + 4] - A[rA + 3][cA + 0] + A[rA + 3][cA + 2] - A[rA + 3][cA + 3] + A[rA + 3][cA + 4]) * (-B[rB + 4][cB + 0] - B[rB + 4][cB + 2] - B[rB + 4][cB + 4]);
    op_count += 10;
    H[71] = (-A[rA + 1][cA + 0] - A[rA + 1][cA + 3] - A[rA + 3][cA + 0] - A[rA + 3][cA + 3]) * (B[rB + 3][cB + 0] + B[rB + 3][cB + 1] + B[rB + 3][cB + 2]);
    op_count += 6;
    H[72] = (A[rA + 0][cA + 2] - A[rA + 0][cA + 3] - A[rA + 0][cA + 4] + A[rA + 1][cA + 2] - A[rA + 1][cA + 3] - A[rA + 1][cA + 4]) * (B[rB + 0][cB + 0] + B[rB + 0][cB + 1] - B[rB + 0][cB + 3] + B[rB + 1][cB + 3] + B[rB + 4][cB + 1] - B[rB + 4][cB + 3]);
    op_count += 11;
    H[73] = (A[rA + 1][cA + 0] - A[rA + 1][cA + 2] + A[rA + 1][cA + 3] - A[rA + 2][cA + 0] + A[rA + 2][cA + 2] - A[rA + 2][cA + 3]) * (B[rB + 3][cB + 0] + B[rB + 3][cB + 1] + B[rB + 3][cB + 4]);
    op_count += 8;
    H[74] = -(A[rA + 0][cA + 1] + A[rA + 0][cA + 3] - A[rA + 1][cA + 1] - A[rA + 1][cA + 4] - A[rA + 2][cA + 0] + A[rA + 2][cA + 1] + A[rA + 2][cA + 3] + A[rA + 2][cA + 4] - A[rA + 3][cA + 0] + A[rA + 3][cA + 1]) * B[rB + 1][cB + 4];
    op_count += 10;
    H[75] = (A[rA + 0][cA + 2] + A[rA + 2][cA + 2]) * (-B[rB + 0][cB + 0] + B[rB + 0][cB + 3] - B[rB + 0][cB + 4] + B[rB + 1][cB + 3] + B[rB + 2][cB + 3] - B[rB + 2][cB + 4]);
    op_count += 7;

    C[rC + 0][cC + 0] = -H[9] + H[11] + H[13] - H[14] - H[15] + H[52] + H[4] - H[65] - H[6];
    op_count += 8;
    C[rC + 1][cC + 0] = H[9] + H[10] - H[11] + H[12] + H[14] + H[15] - H[16] - H[43] + H[50];
    op_count += 8;
    C[rC + 2][cC + 0] = H[9] - H[11] + H[14] + H[15] - H[0] + H[1] + H[2] - H[3] + H[74];
    op_count += 8;
    C[rC + 3][cC + 0] = -H[9] + H[11] - H[14] - H[15] + H[51] + H[53] - H[5] - H[7] + H[8];
    op_count += 8;
    C[rC + 0][cC + 1] = H[12] + H[14] + H[19] + H[20] - H[21] + H[22] + H[24] - H[42] + H[48] + H[49];
    op_count += 9;
    C[rC + 1][cC + 1] = -H[10] + H[11] - H[12] - H[14] - H[15] + H[16] + H[17] - H[18] - H[20] + H[42] + H[43];
    op_count += 10;
    C[rC + 2][cC + 1] = -H[15] - H[18] - H[20] - H[27] - H[28] - H[37] + H[41] + H[43] - H[46] + H[47];
    op_count += 9;
    C[rC + 3][cC + 1] = H[10] - H[11] - H[17] + H[20] - H[31] + H[32] - H[33] - H[35] + H[61] - H[69];
    op_count += 9;
    C[rC + 0][cC + 2] = H[14] + H[22] + H[23] + H[33] - H[36] + H[39] - H[40] + H[54] - H[55] - H[8];
    op_count += 9;
    C[rC + 1][cC + 2] = -H[9] + H[18] + H[31] + H[34] + H[35] + H[36] - H[42] - H[59] - H[5] - H[71];
    op_count += 9;
    C[rC + 2][cC + 2] = -H[15] - H[27] + H[32] + H[36] - H[38] + H[44] - H[45] + H[62] - H[70] - H[7];
    op_count += 9;
    C[rC + 3][cC + 2] = H[9] + H[14] + H[15] - H[32] + H[33] - H[34] - H[36] - H[53] + H[5] + H[7] - H[8];
    op_count += 10;
    C[rC + 0][cC + 3] = -H[9] + H[11] + H[13] - H[15] + H[22] + H[23] + H[24] + H[25] + H[4] - H[65] - H[6];
    op_count += 10;
    C[rC + 1][cC + 3] = H[9] + H[17] - H[18] + H[19] - H[21] - H[23] - H[25] - H[4] - H[68] + H[72];
    op_count += 9;
    C[rC + 2][cC + 3] = -H[13] + H[15] - H[22] - H[25] + H[26] + H[28] + H[30] + H[45] - H[57] + H[75];
    op_count += 9;
    C[rC + 3][cC + 3] = H[11] + H[24] + H[25] - H[32] - H[34] - H[39] + H[40] + H[64] - H[67] - H[6];
    op_count += 9;
    C[rC + 0][cC + 4] = H[14] + H[23] + H[24] + H[26] - H[27] + H[29] + H[30] - H[3] + H[60] + H[63];
    op_count += 9;
    C[rC + 1][cC + 4] = -H[9] - H[17] - H[1] - H[29] - H[37] + H[41] - H[42] + H[45] + H[66] + H[73];
    op_count += 9;
    C[rC + 2][cC + 4] = -H[9] + H[11] - H[14] + H[27] + H[28] - H[1] - H[29] - H[2] + H[45] + H[3] - H[74];
    op_count += 10;
    C[rC + 3][cC + 4] = -H[11] - H[28] + H[29] - H[33] + H[34] + H[38] + H[2] - H[44] + H[56] + H[58];
    op_count += 9;
}

void multiply_ai_recursive_inplace(Matrix &C, int rC, int cC,
                                   const Matrix &A, int rA, int cA,
                                   const Matrix &B, int rB, int cB,
                                   int M, int K, int P, unsigned long long &op_count)
{
    const int M_BASE = 4;
    const int K_BASE = 5;
    const int P_BASE = 5;

    if (M == M_BASE && K == K_BASE && P == P_BASE)
    {
        matrix_ai_inplace(C, rC, cC, A, rA, cA, B, rB, cB, op_count);
        return;
    }

    if (M % 2 != 0 || K % 2 != 0 || P % 2 != 0)
    {
        std::cerr << "Error: Matrix dimensions (" << M << "x" << K << ") * ("
                  << K << "x" << P << "). "
                  << "Matrix does not have appropriate dimensions for this algorithm." << std::endl;
        throw std::invalid_argument("Matrix does not have appropriate dimensions for this algorithm.");
    }

    int m_split = M / 2;
    int k_split = K / 2;
    int p_split = P / 2;

    Matrix c11_p1 = createMatrix(m_split, p_split);
    Matrix c11_p2 = createMatrix(m_split, p_split);
    Matrix c12_p1 = createMatrix(m_split, p_split);
    Matrix c12_p2 = createMatrix(m_split, p_split);
    Matrix c21_p1 = createMatrix(m_split, p_split);
    Matrix c21_p2 = createMatrix(m_split, p_split);
    Matrix c22_p1 = createMatrix(m_split, p_split);
    Matrix c22_p2 = createMatrix(m_split, p_split);

    // A11*B11
    multiply_ai_recursive_inplace(c11_p1, 0, 0, A, rA, cA, B, rB, cB, m_split, k_split, p_split, op_count);
    // A12*B21
    multiply_ai_recursive_inplace(c11_p2, 0, 0, A, rA, cA + k_split, B, rB + k_split, cB, m_split, k_split, p_split, op_count);
    // A11*B12
    multiply_ai_recursive_inplace(c12_p1, 0, 0, A, rA, cA, B, rB, cB + p_split, m_split, k_split, p_split, op_count);
    // A12*B22
    multiply_ai_recursive_inplace(c12_p2, 0, 0, A, rA, cA + k_split, B, rB + k_split, cB + p_split, m_split, k_split, p_split, op_count);
    // A21*B11
    multiply_ai_recursive_inplace(c21_p1, 0, 0, A, rA + m_split, cA, B, rB, cB, m_split, k_split, p_split, op_count);
    // A22*B21
    multiply_ai_recursive_inplace(c21_p2, 0, 0, A, rA + m_split, cA + k_split, B, rB + k_split, cB, m_split, k_split, p_split, op_count);
    // A21*B12
    multiply_ai_recursive_inplace(c22_p1, 0, 0, A, rA + m_split, cA, B, rB, cB + p_split, m_split, k_split, p_split, op_count);
    // A22*B22
    multiply_ai_recursive_inplace(c22_p2, 0, 0, A, rA + m_split, cA + k_split, B, rB + k_split, cB + p_split, m_split, k_split, p_split, op_count);

    // C11 = c11_p1 + c11_p2
    addMatrices_inplace(C, rC, cC, c11_p1, 0, 0, c11_p2, 0, 0, m_split, p_split, op_count);
    // C12 = c12_p1 + c12_p2
    addMatrices_inplace(C, rC, cC + p_split, c12_p1, 0, 0, c12_p2, 0, 0, m_split, p_split, op_count);
    // C21 = c21_p1 + c21_p2
    addMatrices_inplace(C, rC + m_split, cC, c21_p1, 0, 0, c21_p2, 0, 0, m_split, p_split, op_count);
    // C22 = c22_p1 + c22_p2
    addMatrices_inplace(C, rC + m_split, cC + p_split, c22_p1, 0, 0, c22_p2, 0, 0, m_split, p_split, op_count);
}

Matrix multiply_ai_recursive_wrapper(const Matrix &A, const Matrix &B, unsigned long long &op_count)
{
    size_t M = A.size();
    size_t K = (M > 0) ? A[0].size() : 0;
    size_t P = (B.size() > 0) ? B[0].size() : 0;

    if (K == 0 || K != B.size())
    {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }

    op_count = 0;
    Matrix C = createMatrix(M, P, false);

    multiply_ai_recursive_inplace(C, 0, 0, A, 0, 0, B, 0, 0, M, K, P, op_count);

    return C;
}

/*
int main()
{

    std::cout << "Benchmark Matrix Multiplication (Recursive AI vs Naive)" << std::endl;

    std::vector<AIRresult> results;
    unsigned long long op_count_naive = 0;
    unsigned long long op_count_ai = 0;

    int max_n_level = 5; // Test up to size (4*2^5) x (5*2^5) = 128x160

    double baselineMemoryKB = getPeakPrivateUsageKB();

    for (int n = 0; n <= max_n_level; ++n)
    {
        int M = 4 * static_cast<int>(std::pow(2, n));
        int K = 5 * static_cast<int>(std::pow(2, n));
        int P = 5 * static_cast<int>(std::pow(2, n));

        std::string dims_str = "(" + std::to_string(M) + "x" + std::to_string(K) + ") * (" +
                               std::to_string(K) + "x" + std::to_string(P) + ")";

        std::cout << "\n--- Testing (Level n=" << n << "): " << dims_str << " ---" << std::endl;

        Matrix A = createMatrix(M, K, true);
        Matrix B = createMatrix(K, P, true);

        op_count_naive = 0;
        Matrix C_benchmark = iterativeMultiply(A, B, op_count_naive);

        op_count_ai = 0;
        Matrix C_ai = multiply_ai_recursive_wrapper(A, B, op_count_ai);

        op_count_naive = 0;
        auto start_naive = std::chrono::high_resolution_clock::now();
        iterativeMultiply(A, B, op_count_naive);
        auto end_naive = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_naive = end_naive - start_naive;

        double peakMemoryNaive = getPeakPrivateUsageKB();
        double memory_naive = peakMemoryNaive - baselineMemoryKB;
        baselineMemoryKB = peakMemoryNaive;

        std::cout << "Naive:    Operations: " << std::setw(12) << op_count_naive
                  << ", Time: " << std::fixed << std::setprecision(4) << duration_naive.count() << " ms" << std::endl;
        printf("          Memory (Delta): %.2f KB\n", memory_naive);
        results.push_back({n, dims_str, "Naive", op_count_naive, duration_naive.count(), memory_naive});

        op_count_ai = 0;
        auto start_ai = std::chrono::high_resolution_clock::now();
        multiply_ai_recursive_wrapper(A, B, op_count_ai);
        auto end_ai = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_ai = end_ai - start_ai;

        double peakMemoryAI = getPeakPrivateUsageKB();
        double memory_ai = peakMemoryAI - baselineMemoryKB;
        baselineMemoryKB = peakMemoryAI;

        std::cout << "AI (Rec): Operations: " << std::setw(12) << op_count_ai
                  << ", Time: " << std::fixed << std::setprecision(4) << duration_ai.count() << " ms" << std::endl;
        printf("          Memory (Delta): %.2f KB\n", memory_ai);
        results.push_back({n, dims_str, "AI_Recursive", op_count_ai, duration_ai.count(), memory_ai});
    }

    std::ofstream csvfile("ai_recursive_benchmark.csv");
    csvfile << "n_level,Dimensions,Algorithm,Operations,Duration_ms,Memory_kb\n";
    for (const auto &res : results)
    {
        csvfile << res.n_level << ",\"" << res.dimensions << "\"," << res.algorithm << ","
                << res.operations << "," << res.duration_ms << "," << res.memory_kb << "\n";
    }
    csvfile.close();

    std::cout << "\nBenchmark finished. Results saved to 'ai_recursive_benchmark.csv'" << std::endl;

    return 0;
}*/