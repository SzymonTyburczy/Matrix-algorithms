#pragma once
#include "helperFunctions.h" // Dla Matrix i createMatrix
#include <string>

// Enum do wyboru algorytmu mnożenia
enum class MultiplyAlgorithm
{
    ITERATIVE,
    BINET,
    STRASSEN
};

// Struktura wynikowa dla LU
struct LU_Result
{
    Matrix L;
    Matrix U;
    double determinant;
};

// Deklaracja twojej głównej funkcji z poprzedniego pliku
LU_Result recursive_lu_factorization(
    const Matrix &A,
    unsigned long long &op_count,
    MultiplyAlgorithm algo = MultiplyAlgorithm::ITERATIVE);