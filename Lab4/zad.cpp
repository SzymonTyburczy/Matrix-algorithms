#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <memory>
#include <Eigen/Dense>
#include <fstream>
// Używamy Eigen dla macierzy gęstych i SVD
using namespace Eigen;
using namespace std;

// Parametry kompresji z wykładu
const double EPSILON = 1e-7; // Próg dla wartości osobliwych
const int MAX_RANK = 32;     // Maksymalny rząd dopuszczalny dla liścia
const int MAX_RANK2 = 64;    // Maksymalny rząd dla mnożenia macierzy



// --- STRUKTURA MACIERZY HIERARCHICZNEJ ---

struct HNode {
    bool is_leaf;
    int rows, cols;
    int rank;
    
    // Dla liścia (format skompresowany U * V)
    MatrixXd U;
    MatrixXd V; // Wg slajdu 13: V zawiera już wartości osobliwe (D * V_transposed)

    // Dla węzła wewnętrznego (4 synów)
    std::vector<std::unique_ptr<HNode>> children;

    HNode(int r, int c) : rows(r), cols(c), is_leaf(false), rank(0) {}
};


MatrixXd generate_3d_grid_matrix(int k) {
    int dim = pow(2, k);       // Wymiar w jednej osi
    int N = dim * dim * dim;   // Całkowity rozmiar macierzy N x N
    MatrixXd A = MatrixXd::Zero(N, N);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.1, 1.0); // Losowe wagi

    // Mapowanie (x,y,z) -> indeks wiersza
    auto get_idx = [&](int x, int y, int z) {
        return x + y * dim + z * dim * dim;
    };

    for (int z = 0; z < dim; ++z) {
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                int row = get_idx(x, y, z);
                
                // Przekątna (wierzchołek z samym sobą)
                A(row, row) = dis(gen) * 10.0; 

                // Sąsiedzi (6 kierunków w 3D)
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


std::unique_ptr<HNode> compress_matrix(const MatrixXd& block) {
    int rows = block.rows();
    int cols = block.cols();
    auto node = std::make_unique<HNode>(rows, cols);

    // 1. Sprawdź czy blok jest zerowy (optymalizacja)
    if (block.isZero(EPSILON)) {
        node->is_leaf = true;
        node->rank = 0;
        node->U = MatrixXd::Zero(rows, 0);
        node->V = MatrixXd::Zero(0, cols);
        return node;
    }

    // 2. Oblicz SVD: A = U * D * V^T
    // Używamy BDCSVD dla dużych macierzy, JacobiSVD dla małych
    BDCSVD<MatrixXd> svd(block, ComputeThinU | ComputeThinV);
    VectorXd singular_values = svd.singularValues();
    
    // 3. Warunek dopuszczalności (Admissibility condition) [cite: 298-301]
    int r = 0;
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) > EPSILON) r++;
    }

    // Warunek: r <= max_rank ORAZ r <= size/2
    bool admissible = (r <= MAX_RANK && r <= std::min(rows, cols) / 2);

    if (admissible || rows <= MAX_RANK * 2) { // Jeśli spełniony lub macierz bardzo mała -> LIŚĆ
        node->is_leaf = true;
        node->rank = r;
        
        // v.U = U(:, 1:rank)
        node->U = svd.matrixU().leftCols(r);
        
        // v.V = D(1:rank, 1:rank) * V^T(1:rank, :)
        // Eigen zwraca V, a nie V^T w matrixV(), więc bierzemy transpozycję
        MatrixXd D = singular_values.head(r).asDiagonal();
        node->V = D * svd.matrixV().leftCols(r).transpose();
    } else {
        // Podział rekurencyjny na 4 podmacierze (ćwiartki) [cite: 304]
        node->is_leaf = false;
        int half_r = rows / 2;
        int half_c = cols / 2;

        // Top-Left, Top-Right, Bottom-Left, Bottom-Right
        node->children.push_back(compress_matrix(block.block(0, 0, half_r, half_c)));
        node->children.push_back(compress_matrix(block.block(0, half_c, half_r, cols - half_c)));
        node->children.push_back(compress_matrix(block.block(half_r, 0, rows - half_r, half_c)));
        node->children.push_back(compress_matrix(block.block(half_r, half_c, rows - half_r, cols - half_c)));
    }

    return node;
}

// --- DEKONSTRUKCJA (Rekonstrukcja macierzy gęstej)

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

// --- MNOŻENIE MACIERZ-WEKTOR (Zadanie 1)

VectorXd mv_mult(const HNode* node, const VectorXd& x) {
    if (node->is_leaf) {
        if (node->rank == 0) return VectorXd::Zero(node->rows);
        // y = U * (V * x) -> kolejność ważna dla złożoności O(Nrs)
        return node->U * (node->V * x);
    } else {
        // Podział wektora x na dwie części
        int half_c = node->children[0]->cols;
        VectorXd x1 = x.head(half_c);
        VectorXd x2 = x.tail(x.size() - half_c);

        // Wyliczenie wkładów od synów
        // Wiersze górne: Syn 1 (TL) * x1 + Syn 2 (TR) * x2
        VectorXd y_top = mv_mult(node->children[0].get(), x1) + mv_mult(node->children[1].get(), x2);
        
        // Wiersze dolne: Syn 3 (BL) * x1 + Syn 4 (BR) * x2
        VectorXd y_bottom = mv_mult(node->children[2].get(), x1) + mv_mult(node->children[3].get(), x2);

        // Złożenie wyniku
        VectorXd y(node->rows);
        y << y_top, y_bottom;
        return y;
    }
}

MatrixXd mm_mult_dense_verification(const HNode* A, const HNode* B) {
    // W celach weryfikacji i uproszczenia, w tym zadaniu często dopuszcza się
    // wykonanie operacji na poziomie liści lub dekompresję lokalną, 
    // jeśli nie implementujemy pełnej algebry macierzy H.
    // Poniżej wersja "leniwa" - rekonstrukcja -> mnożenie.
    // UWAGA: Prawdziwa wersja rekurencyjna wymagałaby procedury 'add' z re-SVD.
    
    // Zastępcza implementacja rekurencyjna (hybrydowa):
    if (A->is_leaf && B->is_leaf) {
        return (A->U * A->V) * (B->U * B->V);
    } 
    // Jeśli nie liście, musimy zejść niżej. 
    // Uproszczenie: dekompresujemy bloki do dense, mnożymy i zwracamy.
    // W pełnym rozwiązaniu tutaj nastąpiłaby kompresja wyniku.
    return decompress(A) * decompress(B);
}

// Prawdziwy szkielet rekurencyjny (bez pełnej implementacji 'add' z re-SVD)
// Zgodnie ze slajdem 23:
// C11 = A11*B11 + A12*B21
// C12 = A11*B12 + A12*B22 ... itd.
MatrixXd recursive_mm_mult(const HNode* A, const HNode* B) {
     if (A->is_leaf || B->is_leaf) {
         // Base case: mnożenie macierzy niskiej rangi
         // A * B = (Ua * Va) * (Ub * Vb) = Ua * (Va * Ub) * Vb -> wciąż niska ranga
         // Dla uproszczenia w tym kodzie zwracamy wynik gęsty
         return decompress(A) * decompress(B);
     }
     
     int half_r = A->children[0]->rows;
     int half_c = B->children[0]->cols; // Zakładamy kwadratowe podziały dla uproszczenia
     
     // Obliczamy 4 bloki wynikowe rekurencyjnie
     MatrixXd C11 = recursive_mm_mult(A->children[0].get(), B->children[0].get()) 
                  + recursive_mm_mult(A->children[1].get(), B->children[2].get());
                  
     MatrixXd C12 = recursive_mm_mult(A->children[0].get(), B->children[1].get()) 
                  + recursive_mm_mult(A->children[1].get(), B->children[3].get());
                  
     MatrixXd C21 = recursive_mm_mult(A->children[2].get(), B->children[0].get()) 
                  + recursive_mm_mult(A->children[3].get(), B->children[2].get());
                  
     MatrixXd C22 = recursive_mm_mult(A->children[2].get(), B->children[1].get()) 
                  + recursive_mm_mult(A->children[3].get(), B->children[3].get());
                  
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
// --- GŁÓWNA PĘTLA I TESTY (Zadanie + Raport) ---

void run_experiment(int k) {
    int dim = pow(2, k); 
    int N = dim * dim * dim; // N = 2^(3k)
    
    std::cout << "\n=== EKSPERYMENT DLA k=" << k << " (N=" << N << ") ===" << std::endl;

    // 1. Generowanie
    std::cout << "Generowanie macierzy..." << std::endl;
    MatrixXd A_dense = generate_3d_grid_matrix(k);

    // 2. Kompresja
    std::cout << "Kompresja macierzy..." << std::endl;
    auto start_comp = std::chrono::high_resolution_clock::now();
    auto H_matrix = compress_matrix(A_dense);
    auto end_comp = std::chrono::high_resolution_clock::now();
    std::cout << "Czas kompresji: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count() 
              << " ms" << std::endl;


    if (k == 1 || k == 2 || k == 3 || k == 4 || k == 5) { // Generujemy rysunek np. dla N=512
        std::cout << "Eksportowanie struktury drzewa do matrix_tree.txt..." << std::endl;
        std::string nazwaPliku = "matrix_tree_k" + std::to_string(k) + ".txt";
        std::ofstream outfile(nazwaPliku);
        export_tree_to_file(H_matrix.get(), 0, 0, outfile);
        outfile.close();
        std::cout << "Eksport zakonczony." << std::endl;
        
    }



    // 3. Mnożenie Macierz-Wektor (Zadanie 25/26)
    VectorXd x = VectorXd::Random(N);
    
    // a) Metoda H-Macierz
    auto start_mv = std::chrono::high_resolution_clock::now();
    VectorXd y_H = mv_mult(H_matrix.get(), x);
    auto end_mv = std::chrono::high_resolution_clock::now();
    
    // b) Metoda Gęsta (dla porównania)
    VectorXd y_dense = A_dense * x;

    // c) Weryfikacja błędu (Norma L2) 
    double diff_norm = (y_dense - y_H).squaredNorm(); // Suma kwadratów różnic
    std::cout << "Mnozenie MV czas: " 
              << std::chrono::duration_cast<std::chrono::microseconds>(end_mv - start_mv).count() 
              << " us" << std::endl;
    std::cout << "Blad MV (L2^2): " << diff_norm << std::endl;

    // 4. Mnożenie Macierz-Macierz (A^2) (Zadanie 25/27)
    // Uwaga: Dla k=4 macierz gęsta N=4096 ma 16mln elementów - mnożenie gęste może trwać długo.
    if (k <= 5) { // Ograniczamy test pełnego mnożenia MM dla dużych k żeby nie czekać wiecznie
        std::cout << "Mnozenie MM (A^2)..." << std::endl;
        
        auto start_mm = std::chrono::high_resolution_clock::now();
        MatrixXd A2_H_dense = recursive_mm_mult(H_matrix.get(), H_matrix.get());
        auto end_mm = std::chrono::high_resolution_clock::now();
        
        MatrixXd A2_dense = A_dense * A_dense;

        // Weryfikacja błędu (Norma Frobeniusa) 
        double frob_diff = (A2_dense - A2_H_dense).squaredNorm(); // squaredNorm dla macierzy to norma Frobeniusa^2
        
        std::cout << "Mnozenie MM czas: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_mm - start_mm).count() 
                  << " ms" << std::endl;
        std::cout << "Blad MM (Frobenius^2): " << frob_diff << std::endl;
    } else {
        std::cout << "Pominiecie mnozenia gęstego MM dla k=4 (zbyt duze dla demo)." << std::endl;
    }
}




int main() {
    // Zadanie prosi o wykonanie dla k = 2, 3, 4 
    // k=2 -> N=64
    // k=3 -> N=512
    // k=4 -> N=4096
    
    try {
        run_experiment(2);
        run_experiment(3);
        run_experiment(4);
        // run_experiment(5); // k=5 -> N=32768 (bardzo duża macierz, może zająć dużo czasu i pamięci)
    } catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << std::endl;
    }

    return 0;
}