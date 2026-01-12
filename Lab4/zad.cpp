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

// --- STRUKTURA MACIERZY HIERARCHICZNEJ ---
struct HNode {
    bool is_leaf;
    int rows, cols;
    int rank;
    MatrixXd U;
    MatrixXd V;
    std::vector<std::unique_ptr<HNode>> children;

    HNode(int r, int c) : rows(r), cols(c), is_leaf(false), rank(0) {}
};

// Generowanie macierzy (bez zmian)
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

// ZMODYFIKOWANA KOMPRESJA: przyjmuje max_rank ORAZ epsilon
std::unique_ptr<HNode> compress_matrix(const MatrixXd& block, int max_rank, double epsilon) {
    int rows = block.rows();
    int cols = block.cols();
    auto node = std::make_unique<HNode>(rows, cols);

    if (block.isZero(epsilon)) {
        node->is_leaf = true;
        node->rank = 0;
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

    // Warunek dopuszczalności
    bool admissible = (r <= max_rank && r <= std::min(rows, cols) / 2);

    if (admissible || rows <= max_rank * 2) { 
        node->is_leaf = true;
        node->rank = r;
        node->U = svd.matrixU().leftCols(r);
        MatrixXd D = singular_values.head(r).asDiagonal();
        node->V = D * svd.matrixV().leftCols(r).transpose();
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

// (Standardowe funkcje pomocnicze - kopiowane dla kompletności)
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

MatrixXd recursive_mm_mult(const HNode* A, const HNode* B) {
     if (A->is_leaf || B->is_leaf) return decompress(A) * decompress(B);
     int half_r = A->children[0]->rows;
     int half_c = B->children[0]->cols; 
     MatrixXd C11 = recursive_mm_mult(A->children[0].get(), B->children[0].get()) + recursive_mm_mult(A->children[1].get(), B->children[2].get());
     MatrixXd C12 = recursive_mm_mult(A->children[0].get(), B->children[1].get()) + recursive_mm_mult(A->children[1].get(), B->children[3].get());
     MatrixXd C21 = recursive_mm_mult(A->children[2].get(), B->children[0].get()) + recursive_mm_mult(A->children[3].get(), B->children[2].get());
     MatrixXd C22 = recursive_mm_mult(A->children[2].get(), B->children[1].get()) + recursive_mm_mult(A->children[3].get(), B->children[3].get());
     MatrixXd res(A->rows, B->cols);
     res.block(0,0, half_r, half_c) = C11;
     res.block(0, half_c, half_r, half_c) = C12;
     res.block(half_r, 0, half_r, half_c) = C21;
     res.block(half_r, half_c, half_r, half_c) = C22;
     return res;
}

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

// --- STANDARDOWY EKSPERYMENT (Rank=32, Eps=1e-7) ---
void run_experiment(int k, std::ofstream& results_file) {
    int dim = pow(2, k); 
    long long N = (long long)dim * dim * dim;
    int default_rank = 32;
    double default_eps = 1e-7;

    std::cout << "\n=== STANDARDOWY TEST k=" << k << " (N=" << N << ") ===" << std::endl;
    MatrixXd A_dense = generate_3d_grid_matrix(k);

    auto start_comp = std::chrono::high_resolution_clock::now();
    auto H_matrix = compress_matrix(A_dense, default_rank, default_eps);
    auto end_comp = std::chrono::high_resolution_clock::now();
    long long t_comp = std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count();

    if (k <= 4) {
        std::string nazwaPliku = "matrix_tree_k" + std::to_string(k) + ".txt";
        std::ofstream outfile(nazwaPliku);
        export_tree_to_file(H_matrix.get(), 0, 0, outfile);
    }

    VectorXd x = VectorXd::Random(N);
    auto start_mv = std::chrono::high_resolution_clock::now();
    VectorXd y_H = mv_mult(H_matrix.get(), x);
    auto end_mv = std::chrono::high_resolution_clock::now();
    long long t_mv = std::chrono::duration_cast<std::chrono::microseconds>(end_mv - start_mv).count();
    
    VectorXd y_dense = A_dense * x;
    double err_mv = (y_dense - y_H).squaredNorm();

    long long t_mm = 0;
    double err_mm = -1.0;
    if (k <= 4) { 
        auto start_mm = std::chrono::high_resolution_clock::now();
        MatrixXd A2_H = recursive_mm_mult(H_matrix.get(), H_matrix.get());
        auto end_mm = std::chrono::high_resolution_clock::now();
        t_mm = std::chrono::duration_cast<std::chrono::milliseconds>(end_mm - start_mm).count();
        
        MatrixXd A2_dense = A_dense * A_dense;
        err_mm = (A2_dense - A2_H).squaredNorm();
    }

    results_file << k << " " << N << " " << t_comp << " " << t_mv << " " << err_mv << " " << t_mm << " " << err_mm << "\n";
    std::cout << "Zakonczono k=" << k << " Blad MV: " << err_mv << std::endl;
}

// --- ROZSZERZONA ANALIZA WRAŻLIWOŚCI ---
void run_sensitivity_analysis(std::ofstream& file) {
    std::vector<int> ks = {2, 3}; // k values
    std::vector<int> ranks = {4, 16, 32, 64}; // max_rank values
    std::vector<double> epsilons = {1e-1, 1e-3, 1e-5, 1e-7, 1e-9}; // epsilon values

    std::cout << "\n=== START ROZSZERZONEJ ANALIZY WRAZLIWOSCI (k, rank, epsilon) ===" << std::endl;

    for (int k : ks) {
        int dim = pow(2, k);
        int N = dim * dim * dim;
        std::cout << "Generowanie macierzy dla k=" << k << " (N=" << N << ")..." << std::endl;
        MatrixXd A_dense = generate_3d_grid_matrix(k);
        VectorXd x = VectorXd::Random(N);
        VectorXd y_dense = A_dense * x; // Referencja

        for (double eps : epsilons) {
            for (int r : ranks) {
                // Pomiar kompresji
                auto start_comp = std::chrono::high_resolution_clock::now();
                auto H_matrix = compress_matrix(A_dense, r, eps);
                auto end_comp = std::chrono::high_resolution_clock::now();
                long long t_comp = std::chrono::duration_cast<std::chrono::microseconds>(end_comp - start_comp).count();

                // Pomiar MV
                auto start_mv = std::chrono::high_resolution_clock::now();
                VectorXd y_H = mv_mult(H_matrix.get(), x);
                auto end_mv = std::chrono::high_resolution_clock::now();
                long long t_mv = std::chrono::duration_cast<std::chrono::microseconds>(end_mv - start_mv).count();

                double err_mv = (y_dense - y_H).squaredNorm();

                // Zapis: k rank epsilon t_comp(us) t_mv(us) err_mv
                file << k << " " << r << " " << eps << " " << t_comp << " " << t_mv << " " << err_mv << "\n";
                
                // Prosty progress bar w konsoli
                std::cout << "k=" << k << " eps=" << eps << " r=" << r << " -> err=" << err_mv << "\r" << std::flush;
            }
        }
        std::cout << std::endl; // Nowa linia po każdym k
    }
    std::cout << "Analiza zakonczona." << std::endl;
}

// ... (Wszystkie funkcje pomocnicze i run_experiment zostają bez zmian) ...

void export_specific_case(int k, int rank, double eps) {
    int dim = pow(2, k);
    long long N = (long long)dim * dim * dim;
    
    std::cout << "Generowanie wizualizacji: k=" << k << " r=" << rank << " eps=" << eps << "... ";
    MatrixXd A = generate_3d_grid_matrix(k);
    auto H = compress_matrix(A, rank, eps);
    
    // Tworzymy nazwę pliku z parametrami
    std::ostringstream filename;
    filename << "matrix_tree_k" << k << "_r" << rank << "_eps" << eps << ".txt";
    
    std::ofstream outfile(filename.str());
    export_tree_to_file(H.get(), 0, 0, outfile);
    outfile.close();
    std::cout << "Zapisano: " << filename.str() << std::endl;
}

int main() {
    try {
        // 1. Główne zadanie (results.txt)
        std::ofstream results_file("results.txt");
        if (results_file.is_open()) {
            run_experiment(2, results_file);
            run_experiment(3, results_file); // Standardowe k=3
            run_experiment(4, results_file);
            results_file.close();
        }

        // 2. Analiza Wrażliwości (results_sensitivity.txt)
        std::ofstream sens_file("results_sensitivity.txt");
        if (sens_file.is_open()) {
            run_sensitivity_analysis(sens_file);
            sens_file.close();
        }

        // --- 3. NOWOŚĆ: EKSPORT RÓŻNYCH STRUKTUR DO WIZUALIZACJI ---
        // Generujemy 3 przypadki dla k=3 (N=512), żeby pokazać różnice w raporcie
        
        // Przypadek A: Bardzo mocna kompresja (Mało szczegółów)
        // Rank=2, Epsilon=0.1 (duży błąd, mało czarnych pól)
        export_specific_case(3, 2, 0.1);

        // Przypadek B: Zbalansowany (Standard)
        // Rank=32, Epsilon=1e-7
        export_specific_case(3, 256, 1e-6);

        // Przypadek C: Bardzo wysoka precyzja (Prawie jak gęsta)
        // Rank=128, Epsilon=1e-9 (dużo detali)
        export_specific_case(3, 128, 1e-9);

        std::cout << "\nGotowe! Wszystkie pliki wygenerowane." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << std::endl;
    }
    return 0;
}