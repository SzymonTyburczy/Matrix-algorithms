#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <windows.h>
#include <psapi.h>
#include <iomanip>
#include <fstream>
#include <cstdio>

using Matrix = std::vector<std::vector<double>>;

struct resultsinCSV
{
    int size;
    std::string algorithm;
    unsigned long long operations;
    double duration_ms;
    double memory_kb;

    resultsinCSV(int s, std::string algo, unsigned long long ops, double d, double mem)
        : size(s), algorithm(algo), operations(ops), duration_ms(d), memory_kb(mem) {}
};

Matrix createMatrix(size_t rows, size_t cols, bool random = false);
Matrix iterativeMultiply(const Matrix &A, const Matrix &B, unsigned long long &op_count);

Matrix addMatrices(const Matrix &A, int rA, int cA,
                   const Matrix &B, int rB, int cB,
                   int rows, int cols, unsigned long long &op_count);

Matrix subtractMatrices(const Matrix &A, int rA, int cA,
                        const Matrix &B, int rB, int cB,
                        int rows, int cols, unsigned long long &op_count);

void addMatrices_inplace(Matrix &C, int rC, int cC,
                         const Matrix &A, int rA, int cA,
                         const Matrix &B, int rB, int cB,
                         int rows, int cols, unsigned long long &op_count);

void subtractMatrices_inplace(Matrix &C, int rC, int cC,
                              const Matrix &A, int rA, int cA,
                              const Matrix &B, int rB, int cB,
                              int rows, int cols, unsigned long long &op_count);

Matrix addMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);
Matrix subtractMatrices(const Matrix &A, const Matrix &B, unsigned long long &op_count);

double printMemoryUsage();
double getPeakPrivateUsageKB();
void printMatrix(const std::vector<std::vector<double>> &matrix);