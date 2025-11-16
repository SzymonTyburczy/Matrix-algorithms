#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility> // dla std::swap
#include <iomanip>
#include <string>
#include <chrono> // dla pomiaru czasu

// Załączamy nagłówek z Twoimi funkcjami
#include "helperFunctions.h"

const double EPS = 1e-12;

// =================================================================
// === KROK 1: DODATKOWE FUNKCJE POMOCNICZE (WYCINANIE BLOKÓW) ===
// =================================================================
// Te funkcje są niezbędne do kopiowania bloków A11, A12, b1
// do tymczasowych macierzy roboczych.

// Wyciąga pod-macierz z A (tworzy nową macierz)
Matrix get_submatrix(const Matrix &A, int r_start, int c_start, int num_rows, int num_cols)
{
    Matrix sub = createMatrix(num_rows, num_cols); // Używamy Twojej funkcji
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            sub[i][j] = A[r_start + i][c_start + j];
        }
    }
    return sub;
}

// Wyciąga pod-wektor (jako wektor) z b
std::vector<double> get_subvector(const std::vector<double> &b, int r_start, int num_rows)
{
    std::vector<double> sub(num_rows);
    for (int i = 0; i < num_rows; ++i)
    {
        sub[i] = b[r_start + i];
    }
    return sub;
}

// =================================================================
// === KROK 2: SOLWER PUNKTOWY (PRZYPADEK BAZOWY) ===
// =================================================================
// To jest nasz stary solwer 1x1, działający w miejscu.
// Będzie używany, gdy bloki staną się wystarczająco małe.

std::vector<double> solve_pointwise_internal(Matrix &A, std::vector<double> &b,
                                             unsigned long long &flop_count, int offset, int n_total)
{
    int current_size = n_total - offset;

    // --- Baza rekurencji (1x1) ---
    if (current_size == 1)
    {
        if (fabs(A[offset][offset]) < EPS)
            throw std::runtime_error("Zero pivot");
        flop_count++;
        return std::vector<double>{b[offset] / A[offset][offset]};
    }

    // --- Pivoting ---
    int piv = offset;
    double maxv = fabs(A[offset][offset]);
    for (int i = offset + 1; i < n_total; i++)
    {
        double av = fabs(A[i][offset]);
        if (av > maxv)
        {
            maxv = av;
            piv = i;
        }
    }
    if (maxv < EPS)
        throw std::runtime_error("Matrix singular (pointwise)");
    if (piv != offset)
    {
        std::swap(A[offset], A[piv]); // Zamiana całych wierszy
        std::swap(b[offset], b[piv]);
    }

    // --- Eliminacja ---
    for (int i = offset + 1; i < n_total; i++)
    {
        double factor = A[i][offset] / A[offset][offset];
        flop_count++;
        A[i][offset] = 0.0; // Zerowanie "w miejscu"
        for (int j = offset + 1; j < n_total; j++)
        {
            A[i][j] -= factor * A[offset][j];
            flop_count += 2;
        }
        b[i] -= factor * b[offset];
        flop_count += 2;
    }

    // --- Wywołanie rekurencyjne ---
    std::vector<double> x_sub = solve_pointwise_internal(A, b, flop_count, offset + 1, n_total);

    // --- Podstawienie wsteczne ---
    std::vector<double> x(current_size);
    double s = 0.0;
    for (int j = offset + 1; j < n_total; j++)
    {
        s += A[offset][j] * x_sub[j - (offset + 1)];
        flop_count += 2;
    }
    x[0] = (b[offset] - s) / A[offset][offset];
    flop_count += 2;
    for (int i = 1; i < current_size; i++)
        x[i] = x_sub[i - 1];
    return x;
}

// Opakowanie dla solwera punktowego (tworzy kopie)
std::vector<double> solve_pointwise_gauss(Matrix A, std::vector<double> b, unsigned long long &flop_count)
{
    return solve_pointwise_internal(A, b, flop_count, 0, A.size());
}

// =================================================================
// === KROK 3: NOWY, ZOPTAMALIZOWANY SOLWER BLOKOWY ===
// =================================================================

std::vector<double> solve_block_recursive(Matrix &A, std::vector<double> &b,
                                          unsigned long long &flop_count, int offset, int block_size)
{
    int n_total = A.size();
    int current_size = n_total - offset;

    // --- BAZA REKURENCJI: Problem jest mniejszy niż blok ---
    if (current_size <= block_size)
    {
        // Użyj solwera punktowego "w miejscu" dla reszty
        return solve_pointwise_internal(A, b, flop_count, offset, n_total);
    }

    // --- KROK BLOKOWY ---
    int b_size = block_size;
    int r_size = current_size - b_size;

    // Definicje offsetów dla bloków
    int r11 = offset;
    int c11 = offset;
    int r12 = offset;
    int c12 = offset + b_size;
    int r21 = offset + b_size;
    int c21 = offset;
    int r22 = offset + b_size;
    int c22 = offset + b_size;

    // 1. Rozwiąż A11 * X = A12   i   A11 * y = b1

    // Kopiujemy A11. Jest to KONIECZNE, ponieważ solve_pointwise_gauss
    // jest destrukcyjny (modyfikuje macierz, którą dostaje).
    Matrix A11 = get_submatrix(A, r11, c11, b_size, b_size);

    // Obliczamy X = A11_inv * A12
    // X będzie przechowywać wynik operacji
    Matrix X = createMatrix(b_size, r_size);
    for (int j = 0; j < r_size; ++j)
    {
        // Wycinamy j-tą kolumnę A12
        std::vector<double> A12_col(b_size);
        for (int i = 0; i < b_size; ++i)
            A12_col[i] = A[r12 + i][c12 + j];

        // Rozwiązujemy A11 * x_col = A12_col.
        // `solve_pointwise_gauss` tworzy WŁASNĄ KOPIĘ A11, więc
        // nasza macierz `A11` może być ponownie użyta.
        std::vector<double> x_col = solve_pointwise_gauss(A11, A12_col, flop_count);

        // Zapisz wynik do macierzy X
        for (int i = 0; i < b_size; ++i)
            X[i][j] = x_col[i];
    }

    // Obliczamy y = A11_inv * b1
    std::vector<double> b1 = get_subvector(b, r11, b_size);
    std::vector<double> y = solve_pointwise_gauss(A11, b1, flop_count);

    // 2. Aktualizuj A22 = A22 - A21 * X  (Dopełnienie Schura)
    // Tworzymy macierz tymczasową na wynik A21 * X
    Matrix A21_X = createMatrix(r_size, r_size);

    // Używamy Twojej funkcji: A21_X = A21 * X
    iterativeMultiply_inplace(A21_X, 0, 0,            // C
                              A, r21, c21,            // A (A21)
                              X, 0, 0,                // B (X)
                              r_size, b_size, r_size, // m, k, p
                              flop_count);

    // Używamy Twojej funkcji: A22 = A22 - A21_X (w miejscu)
    subtractMatrices_inplace(A, r22, c22, // C (A22)
                             A, r22, c22, // A (A22)
                             A21_X, 0, 0, // B (A21_X)
                             r_size, r_size,
                             flop_count);

    // 3. Aktualizuj b2 = b2 - A21 * y
    // Tworzymy wektor tymczasowy na wynik A21 * y
    std::vector<double> A21_y(r_size);
    // Niestety, nie masz `matrix_vector_mult_inplace`,
    // więc robimy to ręcznie (nadal wydajnie, bez kopiowania A21)
    for (int i = 0; i < r_size; ++i)
    {
        double sum = 0.0;
        for (int k = 0; k < b_size; ++k)
        {
            sum += A[r21 + i][c21 + k] * y[k];
            flop_count += 2;
        }
        A21_y[i] = sum;
    }

    // Aktualizujemy b2 "w miejscu"
    for (int i = 0; i < r_size; ++i)
    {
        b[r21 + i] -= A21_y[i];
        flop_count++;
    }

    // 4. Rozwiąż rekurencyjnie dla (A22_new, b2_new) aby znaleźć x2
    std::vector<double> x2 = solve_block_recursive(A, b, flop_count, r22, block_size);

    // 5. Podstawienie wsteczne: Oblicz x1 = y - X * x2
    // Obliczamy X * x2 (też ręcznie, bo brak helpera)
    std::vector<double> X_x2(b_size);
    for (int i = 0; i < b_size; ++i)
    {
        double sum = 0.0;
        for (int k = 0; k < r_size; ++k)
        {
            sum += X[i][k] * x2[k];
            flop_count += 2;
        }
        X_x2[i] = sum;
    }

    // Obliczamy x1 = y - X_x2
    std::vector<double> x1 = y; // Kopiujemy y
    for (int i = 0; i < b_size; ++i)
    {
        x1[i] -= X_x2[i];
        flop_count++;
    }

    // 6. Połącz x1 i x2 w jeden wektor wynikowy
    std::vector<double> x_final(current_size);
    for (int i = 0; i < b_size; ++i)
        x_final[i] = x1[i];
    for (int i = 0; i < r_size; ++i)
        x_final[b_size + i] = x2[i];

    return x_final;
}

// Opakowanie dla solwera blokowego
std::vector<double> solve_block_gauss(Matrix A, std::vector<double> b,
                                      unsigned long long &flop_count, int block_size)
{
    if (block_size <= 0)
        throw std::invalid_argument("Block size must be > 0");
    if (A.size() != b.size())
        throw std::invalid_argument("Matrix/vector size mismatch");

    // Przekazujemy A i b przez WARTOŚĆ (kopia), aby solve_block_recursive
    // mogła je bezpiecznie modyfikować "w miejscu" bez niszczenia oryginału w main.
    return solve_block_recursive(A, b, flop_count, 0, block_size);
}

// =================================================================
// === KROK 4: FUNKCJA main (TESTOWANIE) ===
// =================================================================

int main()
{
    // Ustawiamy rozmiar macierzy i bloku
    int N = 8;
    int BLOCK_SIZE = 2; // Rozmiar musi być dzielnikiem N dla tej prostej wersji

    if (N % BLOCK_SIZE != 0)
    {
        std::cerr << "Rozmiar macierzy musi byc wielokrotnoscia rozmiaru bloku." << std::endl;
        // W bardziej złożonej implementacji obsłużylibyśmy "resztki"
        // ale dla testu upraszczamy.
        return 1;
    }

    // Tworzymy macierze używając Twojej funkcji
    Matrix A = createMatrix(N, N, true);
    // Tworzymy wektor b i wypełniamy go (np. suma wierszy A)
    std::vector<double> b(N);
    for (int i = 0; i < N; ++i)
    {
        double sum = 0;
        for (int j = 0; j < N; ++j)
            sum += A[i][j];
        b[i] = sum;
        // Jeśli b = A * [1, 1, ..., 1]^T, to rozwiązaniem x powinno być [1, 1, ..., 1]^T
    }

    std::cout << "Test: Blokowa eliminacja Gaussa (rozmiar "
              << N << "x" << N << ", blok " << BLOCK_SIZE << ")" << std::endl;

    unsigned long long total_flops = 0;
    double mem_usage_kb = 0;

    // Resetujemy licznik pamięci
    getPeakPrivateUsageKB();

    // Mierzymy czas
    auto start_time = std::chrono::high_resolution_clock::now();
    try
    {
        std::vector<double> x = solve_block_gauss(A, b, total_flops, BLOCK_SIZE);

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration_ms = end_time - start_time;

        mem_usage_kb = getPeakPrivateUsageKB();

        std::cout << "\nRozwiazanie x (powinno byc bliskie 1.0):" << std::endl;
        print_vector(x); // Używamy Twojej funkcji
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Laczna liczba operacji (FLOPS): " << total_flops << std::endl;
        std::cout << "Czas wykonania: " << duration_ms.count() << " ms" << std::endl;
        std::cout << "Szczytowe zuzycie pamieci: " << mem_usage_kb << " KB" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Blad: " << e.what() << std::endl;
    }

    return 0;
}