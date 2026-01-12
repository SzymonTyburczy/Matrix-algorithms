#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <memory>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

using namespace Eigen;
using namespace std;

// Parametry globalne (domyślne)
// W analizie wrażliwości mogą być nadpisywane przez parametry funkcji, 
// ale tutaj używamy ich do metod arytmetycznych.
const double EPSILON = 1e-7; 
const int MAX_RANK = 32;     

// --- STRUKTURA MACIERZY HIERARCHICZNEJ ---
struct HNode {
    bool is_leaf;
    int rows, cols;
    int rank;
    
    // Format skompresowany: M = U * V
    MatrixXd U;
    MatrixXd V; 

    // Synowie (4 ćwiartki)
    std::vector<std::unique_ptr<HNode>> children;

    HNode(int r, int c) : rows(r), cols(c), is_leaf(false), rank(0) {}
};

// --- GENEROWANIE DANYCH ---
MatrixXd generate_3d_grid_matrix(int k) {
    int dim = pow(2, k);       
    int N = dim * dim * dim;   
    MatrixXd A = MatrixXd::Zero(N, N);

    std::random_device rd;
    std::mt19937 gen(42); 
    std::uniform_real_distribution<> dis(0.1, 1.0);

    auto get_idx = [&](int x, int y, int z) {
        return x + y * dim + z * dim * dim;
    };

    for (int z = 0; z < dim; ++z) {
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                int row = get_idx(x, y, z);
                A(row, row) = dis(gen) * 10.0; 
                int dx[] = {1, -1, 0, 0, 0, 0};
                int dy[] = {0, 0, 1, -1, 0, 0};
                int dz[] = {0, 0, 0, 0, 1, -1};

                for (int i = 0; i < 6; ++i) {
                    int nx = x + dx[i];
                    int ny = y + dy[i];
                    int nz = z + dz[i];
                    if (nx >= 0 && nx < dim && ny >= 0 && ny < dim && nz >= 0 && nz < dim) {
                        int col = get_idx(nx, ny, nz);
                        A(row, col) = dis(gen);
                    }
                }
            }
        }
    }
    return A;
}

// --- KOMPRESJA I DEKOMPRESJA ---

// Kompresja z parametrami (dla analizy wrażliwości)
std::unique_ptr<HNode> compress_matrix(const MatrixXd& block, int max_rank = MAX_RANK, double epsilon = EPSILON) {
    int rows = block.rows();
    int cols = block.cols();
    auto node = std::make_unique<HNode>(rows, cols);

    if (block.isZero(epsilon)) {
        node->is_leaf = true; node->rank = 0;
        node->U = MatrixXd::Zero(rows, 0);
        node->V = MatrixXd::Zero(0, cols);
        return node;
    }

    BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
    VectorXd singular_values = svd.singularValues();
    
    int r = 0;
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) > epsilon) r++;
    }

    bool admissible = (r <= max_rank && r <= std::min(rows, cols) / 2);

    if (admissible || rows <= max_rank * 2) { 
        node->is_leaf = true;
        node->rank = r;
        node->U = svd.matrixU().leftCols(r);
        node->V = singular_values.head(r).asDiagonal() * svd.matrixV().leftCols(r).transpose();
    } else {
        node->is_leaf = false;
        int half_r = rows / 2;
        int half_c = cols / 2;
        node->children.push_back(compress_matrix(block.block(0, 0, half_r, half_c), max_rank, epsilon));
        node->children.push_back(compress_matrix(block.block(0, half_c, half_r, cols - half_c), max_rank, epsilon));
        node->children.push_back(compress_matrix(block.block(half_r, 0, rows - half_r, half_c), max_rank, epsilon));
        node->children.push_back(compress_matrix(block.block(half_r, half_c, rows - half_r, cols - half_c), max_rank, epsilon));
    }
    return node;
}

MatrixXd decompress(const HNode* node) {
    if (node->is_leaf) {
        if (node->rank == 0) return MatrixXd::Zero(node->rows, node->cols);
        return node->U * node->V;
    } else {
        MatrixXd res(node->rows, node->cols);
        int half_r = node->children[0]->rows;
        int half_c = node->children[0]->cols;
        res.block(0, 0, half_r, half_c) = decompress(node->children[0].get());
        res.block(0, half_c, half_r, node->cols - half_c) = decompress(node->children[1].get());
        res.block(half_r, 0, node->rows - half_r, half_c) = decompress(node->children[2].get());
        res.block(half_r, half_c, node->rows - half_r, node->cols - half_c) = decompress(node->children[3].get());
        return res;
    }
}

// --- MNOŻENIE MACIERZ-WEKTOR ---
VectorXd mv_mult(const HNode* node, const VectorXd& x) {
    if (node->is_leaf) {
        if (node->rank == 0) return VectorXd::Zero(node->rows);
        return node->U * (node->V * x);
    } else {
        int half_c = node->children[0]->cols;
        VectorXd x1 = x.head(half_c);
        VectorXd x2 = x.tail(x.size() - half_c);
        VectorXd y_top = mv_mult(node->children[0].get(), x1) + mv_mult(node->children[1].get(), x2);
        VectorXd y_bottom = mv_mult(node->children[2].get(), x1) + mv_mult(node->children[3].get(), x2);
        VectorXd y(node->rows);
        y << y_top, y_bottom;
        return y;
    }
}

// --- PEŁNA ARYTMETYKA H-MACIERZY (Dodawanie i Mnożenie) ---

// Funkcja pomocnicza: Dodawanie dwóch liści z RE-KOMPRESJĄ
// C = A + B. Jeśli A i B mają rząd k, to C ma rząd max 2k.
// Wykonujemy SVD na sumie, żeby przyciąć rząd z powrotem do MAX_RANK.
std::unique_ptr<HNode> add_leaves(const HNode* A, const HNode* B) {
    auto node = std::make_unique<HNode>(A->rows, A->cols);
    node->is_leaf = true;

    // Przypadki trywialne (gdy jeden składnik jest zerowy)
    if (A->rank == 0 && B->rank == 0) { node->rank = 0; return node; }
    if (A->rank == 0) { node->rank = B->rank; node->U = B->U; node->V = B->V; return node; }
    if (B->rank == 0) { node->rank = A->rank; node->U = A->U; node->V = A->V; return node; }

    // Odtwarzamy gęste bloki (są małe lub niskiego rzędu, więc to szybkie)
    // W profesjonalnych bibliotekach używa się arytmetyki na faktorach U/V (QR), 
    // ale rekonstrukcja liścia jest w pełni akceptowalna i prostsza.
    MatrixXd DenseA = A->U * A->V;
    MatrixXd DenseB = B->U * B->V;
    MatrixXd Sum = DenseA + DenseB;

    // Re-kompresja sumy (Truncated SVD)
    BDCSVD<MatrixXd> svd(Sum, ComputeThinU | ComputeThinV);
    VectorXd sv = svd.singularValues();
    
    int r = 0;
    for (int i = 0; i < sv.size(); ++i) {
        if (sv(i) > EPSILON) r++;
    }
    if (r > MAX_RANK) r = MAX_RANK; // Przycinamy rząd, żeby nie rósł

    node->rank = r;
    if (r > 0) {
        node->U = svd.matrixU().leftCols(r);
        node->V = sv.head(r).asDiagonal() * svd.matrixV().leftCols(r).transpose();
    } else {
        node->U = MatrixXd::Zero(A->rows, 0);
        node->V = MatrixXd::Zero(0, A->cols);
    }
    
    return node;
}

// Rekurencyjne dodawanie H-macierzy: C = A + B
std::unique_ptr<HNode> h_add(const HNode* A, const HNode* B) {
    // 1. Oba są liśćmi -> używamy funkcji pomocniczej add_leaves
    if (A->is_leaf && B->is_leaf) {
        return add_leaves(A, B);
    }
    
    // 2. Jeśli struktura się zgadza (oba są węzłami), schodzimy rekurencyjnie
    if (!A->is_leaf && !B->is_leaf) {
        auto node = std::make_unique<HNode>(A->rows, A->cols);
        node->is_leaf = false;
        for (int i = 0; i < 4; ++i) {
            node->children.push_back(h_add(A->children[i].get(), B->children[i].get()));
        }
        return node;
    } 
    
    // 3. Przypadek mieszany (Liść + Węzeł). 
    // W tym zadaniu (siatka regularna) rzadko występuje, ale dla kompletności:
    // Dekompresujemy oba do gęstej i kompresujemy wynik. (Fallback)
    MatrixXd denseSum = decompress(A) + decompress(B);
    return compress_matrix(denseSum);
}

// Rekurencyjne mnożenie H-macierzy: C = A * B
// Zwraca strukturę HNode (H-macierz), a nie macierz gęstą!
std::unique_ptr<HNode> h_mult(const HNode* A, const HNode* B) {
    
    // 1. Mnożenie liści (Low-Rank * Low-Rank) -> Low-Rank
    if (A->is_leaf && B->is_leaf) {
        auto node = std::make_unique<HNode>(A->rows, B->cols);
        node->is_leaf = true;
        
        if (A->rank == 0 || B->rank == 0) {
            node->rank = 0;
            node->U = MatrixXd::Zero(A->rows, 0);
            node->V = MatrixXd::Zero(0, B->cols);
            return node;
        }

        // A * B = (Ua * Va) * (Ub * Vb) = Ua * (Va * Ub) * Vb
        // Środek (Va * Ub) jest mały (rank x rank).
        MatrixXd Cross = A->V * B->U; 
        
        // Wynik to Ua * (Cross * Vb). Żeby zachować format U*V, włączamy Cross do V.
        node->rank = A->rank; 
        node->U = A->U; 
        node->V = Cross * B->V; 
        
        return node;
    }

    // 2. Mnożenie węzłów wewnętrznych (Blokowe)
    // C11 = A11*B11 + A12*B21  <-- tu używamy h_mult i h_add
    
    // Jeśli któryś jest liściem (mieszany), dekompresujemy dla uproszczenia (fallback)
    if (A->is_leaf || B->is_leaf) {
        MatrixXd denseMult = decompress(A) * decompress(B);
        return compress_matrix(denseMult);
    }

    auto node = std::make_unique<HNode>(A->rows, B->cols);
    node->is_leaf = false;

    auto& A11 = A->children[0]; auto& A12 = A->children[1];
    auto& A21 = A->children[2]; auto& A22 = A->children[3];
    auto& B11 = B->children[0]; auto& B12 = B->children[1];
    auto& B21 = B->children[2]; auto& B22 = B->children[3];

    // C11 = A11*B11 + A12*B21
    auto P1 = h_mult(A11.get(), B11.get());
    auto P2 = h_mult(A12.get(), B21.get());
    node->children.push_back(h_add(P1.get(), P2.get()));

    // C12 = A11*B12 + A12*B22
    auto P3 = h_mult(A11.get(), B12.get());
    auto P4 = h_mult(A12.get(), B22.get());
    node->children.push_back(h_add(P3.get(), P4.get()));

    // C21 = A21*B11 + A22*B21
    auto P5 = h_mult(A21.get(), B11.get());
    auto P6 = h_mult(A22.get(), B21.get());
    node->children.push_back(h_add(P5.get(), P6.get()));

    // C22 = A21*B12 + A22*B22
    auto P7 = h_mult(A21.get(), B12.get());
    auto P8 = h_mult(A22.get(), B22.get());
    node->children.push_back(h_add(P7.get(), P8.get()));

    return node;
}

// --- EKSPORT I ANALIZY ---

void export_tree_to_file(const HNode* node, int x, int y, std::ofstream& file) {
    file << x << " " << y << " " << node->cols << " " << node->rows 
         << " " << node->rank << " " << (node->is_leaf ? 1 : 0) << "\n";
    if (!node->is_leaf) {
        int half_r = node->children[0]->rows;
        int half_c = node->children[0]->cols;
        export_tree_to_file(node->children[0].get(), x, y, file);
        export_tree_to_file(node->children[1].get(), x + half_c, y, file);
        export_tree_to_file(node->children[2].get(), x, y + half_r, file);
        export_tree_to_file(node->children[3].get(), x + half_c, y + half_r, file);
    }
}

void export_specific_case(int k, int rank, double eps) {
    int dim = pow(2, k);
    long long N = (long long)dim * dim * dim;
    std::cout << "Generowanie wizualizacji: k=" << k << " r=" << rank << " eps=" << eps << "..." << std::endl;
    MatrixXd A = generate_3d_grid_matrix(k);
    auto H = compress_matrix(A, rank, eps);
    std::ostringstream filename;
    filename << "matrix_tree_k" << k << "_r" << rank << "_eps" << eps << ".txt";
    std::ofstream outfile(filename.str());
    export_tree_to_file(H.get(), 0, 0, outfile);
    outfile.close();
}

void run_experiment(int k, std::ofstream& results_file) {
    int dim = pow(2, k); 
    long long N = (long long)dim * dim * dim;
    std::cout << "\n=== STANDARDOWY TEST k=" << k << " (N=" << N << ") ===" << std::endl;
    
    MatrixXd A_dense = generate_3d_grid_matrix(k);

    auto start_comp = std::chrono::high_resolution_clock::now();
    auto H_matrix = compress_matrix(A_dense);
    auto end_comp = std::chrono::high_resolution_clock::now();
    long long t_comp = std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count();

    if (k <= 4) { // Eksport standardowy
        std::ofstream outfile("matrix_tree_k" + std::to_string(k) + ".txt");
        export_tree_to_file(H_matrix.get(), 0, 0, outfile);
    }

    // MV
    VectorXd x = VectorXd::Random(N);
    auto start_mv = std::chrono::high_resolution_clock::now();
    VectorXd y_H = mv_mult(H_matrix.get(), x);
    auto end_mv = std::chrono::high_resolution_clock::now();
    long long t_mv = std::chrono::duration_cast<std::chrono::microseconds>(end_mv - start_mv).count();
    
    VectorXd y_dense = A_dense * x;
    double err_mv = (y_dense - y_H).squaredNorm();

    // MM (Pełna arytmetyka)
    long long t_mm = 0;
    double err_mm = -1.0;
    
    if (k <= 4) { // Teraz możemy robić MM nawet dla k=4, bo jest szybkie (H-matrix)!
        std::cout << "Mnozenie MM (H-arytmetyka)..." << std::endl;
        auto start_mm = std::chrono::high_resolution_clock::now();
        
        // Używamy nowej funkcji h_mult!
        auto H_squared = h_mult(H_matrix.get(), H_matrix.get());
        
        auto end_mm = std::chrono::high_resolution_clock::now();
        t_mm = std::chrono::duration_cast<std::chrono::milliseconds>(end_mm - start_mm).count();
        
        // Weryfikacja tylko dla małych k, bo A_dense*A_dense jest wolne
        if (k <= 3) {
            MatrixXd A2_dense = A_dense * A_dense;
            MatrixXd A2_H_rec = decompress(H_squared.get()); // Rekonstrukcja do sprawdzenia błędu
            err_mm = (A2_dense - A2_H_rec).squaredNorm();
            std::cout << "Blad MM: " << err_mm << std::endl;
        } else {
            std::cout << "Pomijam weryfikacje bledu MM dla k=4 (za duza macierz gesta)." << std::endl;
        }
    }

    results_file << k << " " << N << " " << t_comp << " " << t_mv << " " << err_mv << " " << t_mm << " " << err_mm << "\n";
    std::cout << "Zakonczono k=" << k << ". Czas MM: " << t_mm << " ms" << std::endl;
}

void run_sensitivity_analysis(std::ofstream& file) {
    // Prosta analiza dla N=512 (k=3)
    int k = 3; 
    int dim = pow(2, k); int N = dim * dim * dim;
    std::vector<int> ranks = {4, 16, 32, 64}; 
    std::vector<double> epsilons = {1e-1, 1e-3, 1e-5, 1e-7, 1e-9};

    std::cout << "\n=== ANALIZA WRAZLIWOSCI ===" << std::endl;
    MatrixXd A = generate_3d_grid_matrix(k);
    VectorXd x = VectorXd::Random(N);
    VectorXd y_ref = A * x;

    for (double eps : epsilons) {
        for (int r : ranks) {
            auto start = std::chrono::high_resolution_clock::now();
            auto H = compress_matrix(A, r, eps);
            auto end = std::chrono::high_resolution_clock::now();
            long long t_comp = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            auto start_mv = std::chrono::high_resolution_clock::now();
            VectorXd y = mv_mult(H.get(), x);
            auto end_mv = std::chrono::high_resolution_clock::now();
            long long t_mv = std::chrono::duration_cast<std::chrono::microseconds>(end_mv - start_mv).count();

            double err = (y_ref - y).squaredNorm();
            file << k << " " << r << " " << eps << " " << t_comp << " " << t_mv << " " << err << "\n";
        }
    }
}

int main() {
    try {
        std::ofstream results_file("results.txt");
        if (results_file.is_open()) {
            run_experiment(2, results_file);
            run_experiment(3, results_file);
            run_experiment(4, results_file);
            results_file.close();
        }

        std::ofstream sens_file("results_sensitivity.txt");
        if (sens_file.is_open()) {
            run_sensitivity_analysis(sens_file);
            sens_file.close();
        }

        // Wizualizacje do raportu
        export_specific_case(3, 2, 0.1);
        export_specific_case(3, 32, 1e-7);
        export_specific_case(3, 16, 1e-9);
        export_specific_case(3, 256, 1e-6);

        std::cout << "\nGotowe! Uruchom python3 main.py" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << std::endl;
    }
    return 0;
}