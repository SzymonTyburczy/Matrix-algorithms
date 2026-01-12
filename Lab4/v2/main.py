import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import glob
import sys
import os

# --- 1. RYSOWANIE STRUKTURY (Z POPRAWIONYM TYTUŁEM) ---
def draw_hmatrix():
    print("\n--- Rysowanie struktury H-Macierzy ---")
    # Szukamy wszystkich plików zaczynających się od matrix_tree
    pliki = glob.glob("matrix_tree_*.txt")
    pliki.sort()
    
    if not pliki:
        print("Nie znaleziono plików matrix_tree_*.txt. Uruchom najpierw program C++.")
        return

    for filename in pliki:
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Błąd odczytu pliku {filename}: {e}")
            continue

        if not lines:
            continue
        
        # Pobranie wymiarów z pierwszej linii
        root_data = list(map(int, lines[0].split()))
        matrix_width = root_data[2]
        matrix_height = root_data[3]

        # --- LOGIKA TYTUŁU ---
        # Domyślne wartości ze standardowego eksperymentu w C++
        rank_str = "32"
        eps_str = "1e-7"

        # Sprawdzamy, czy nazwa pliku zawiera niestandardowe parametry
        # np. matrix_tree_k3_r2_eps0.1.txt
        if "_r" in filename and "_eps" in filename:
            try:
                parts = filename.replace(".txt", "").split("_")
                # Szukamy części zaczynających się od 'r' (ale nie root) i 'eps'
                # parts przykładowo: ['matrix', 'tree', 'k3', 'r2', 'eps0.1']
                for p in parts:
                    if p.startswith("r") and p[1:].replace('.', '', 1).isdigit():
                        rank_str = p[1:]
                    if p.startswith("eps"):
                        eps_str = p[3:]
            except:
                pass # W razie błędu parsowania zostają domyślne 32 i 1e-7

        title_suffix = f"\n(Rank={rank_str}, Epsilon={eps_str})"
        # ---------------------

        fig, ax = plt.subplots(figsize=(10, 10))
        
        for line in lines:
            data = list(map(int, line.split()))
            # Format: x, y, w, h, rank, is_leaf
            x, y, w, h, rank, is_leaf = data
            
            if is_leaf:
                if rank > 0:
                    # Blok z danymi (SVD) - czarny
                    rect = patches.Rectangle((x, y), w, h, linewidth=0.5, edgecolor='none', facecolor='black')
                    ax.add_patch(rect)
                else:
                    # Blok pusty (rank 0) - biały (nic nie rysujemy, tło jest białe)
                    pass
            else:
                # Węzły wewnętrzne - szare obramowanie
                rect = patches.Rectangle((x, y), w, h, linewidth=0.2, edgecolor='gray', facecolor='none')
                ax.add_patch(rect)

        ax.set_xlim(0, matrix_width)
        ax.set_ylim(0, matrix_height)
        ax.invert_yaxis() # (0,0) w lewym górnym rogu
        ax.set_aspect('equal')
        
        # Ustawienie tytułu z parametrami
        plt.title(f"Struktura H-Macierzy (N={matrix_width}){title_suffix}", fontsize=14)
        plt.xlabel("Kolumny")
        plt.ylabel("Wiersze")
        
        # Zapis do pliku
        output_img = f"vis_{filename.replace('.txt', '.png')}"
        plt.savefig(output_img, dpi=300)
        print(f"Wygenerowano rysunek: {output_img}")
        plt.close(fig)

# --- 2. WYKRESY WYDAJNOŚCI (Bez zmian) ---
def plot_performance():
    print("\n--- Generowanie wykresów wydajności ---")
    if not os.path.exists("results.txt"):
        print("Brak pliku results.txt.")
        return

    try:
        data = np.loadtxt("results.txt")
        if data.ndim == 1: data = data.reshape(1, -1)
        if data.shape[0] < 2: 
            print("Za mało danych w results.txt.")
            return

        N = data[:, 1]
        t_mv = data[:, 3] # us
        t_mm = data[:, 5] # ms
        
        def draw_log(N_vals, times, title, fname, unit):
            mask = times > 0
            if np.sum(mask) < 2: return
            N_f, T_f = N_vals[mask], times[mask]
            
            plt.figure(figsize=(10,6))
            plt.loglog(N_f, T_f, 'bo-', label='Dane eksperymentalne', linewidth=2)
            
            # Fit
            coeffs = np.polyfit(np.log(N_f), np.log(T_f), 1)
            beta = coeffs[0]
            alpha = np.exp(coeffs[1])
            fit_T = alpha * (N_f ** beta)
            
            plt.loglog(N_f, fit_T, 'r--', label=f'Dopasowanie: $\\beta \\approx {beta:.2f}$', linewidth=1.5)
            plt.title(title, fontsize=14)
            plt.xlabel('Rozmiar macierzy N ($2^{3k}$)', fontsize=12)
            plt.ylabel(f'Czas [{unit}]', fontsize=12)
            plt.grid(True, which="both", ls="-", alpha=0.4)
            plt.legend()
            
            for i, val in enumerate(T_f):
                plt.annotate(f"{int(val)}", (N_f[i], T_f[i]), textcoords="offset points", xytext=(0,10), ha='center')

            plt.savefig(fname, dpi=300)
            print(f"Wygenerowano: {fname} (Beta={beta:.4f})")
            plt.close()
            
        draw_log(N, t_mv, "Złożoność Mnożenia Macierz-Wektor (MV)", "wykres_mv.png", "us")
        draw_log(N, t_mm, "Złożoność Mnożenia Macierz-Macierz (MM)", "wykres_mm.png", "ms")

    except Exception as e: 
        print(f"Błąd przy results.txt: {e}")

# --- 3. ANALIZA WRAŻLIWOŚCI (Bez zmian) ---
def plot_sensitivity_full():
    print("\n--- Generowanie analizy wrażliwości ---")
    if not os.path.exists("results_sensitivity.txt"): return
    
    try:
        data = np.loadtxt("results_sensitivity.txt")
        if data.ndim == 1: data = data.reshape(1, -1)
        
        ks = np.unique(data[:, 0])
        epsilons = np.unique(data[:, 2])

        for k in ks:
            subset_k = data[data[:, 0] == k]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            fig.suptitle(f'Analiza Wrażliwości dla k={int(k)} (N={int(2**(3*k))})', fontsize=16)

            # Wykres błędu
            for eps in epsilons:
                subset_eps = subset_k[subset_k[:, 2] == eps]
                subset_eps = subset_eps[subset_eps[:, 1].argsort()]
                if len(subset_eps) > 0:
                    ax1.semilogy(subset_eps[:, 1], subset_eps[:, 5], 'o--', label=f'Eps={eps:.0e}')
            
            ax1.set_title("Błąd Aproksymacji (L2^2)")
            ax1.set_xlabel("Max Rank"); ax1.set_ylabel("Błąd")
            ax1.grid(True, which="both", alpha=0.3); ax1.legend()

            # Wykres czasu
            for eps in epsilons:
                subset_eps = subset_k[subset_k[:, 2] == eps]
                subset_eps = subset_eps[subset_eps[:, 1].argsort()]
                if len(subset_eps) > 0:
                    ax2.plot(subset_eps[:, 1], subset_eps[:, 4], 's-', label=f'Eps={eps:.0e}')

            ax2.set_title("Czas Mnożenia MV [us]")
            ax2.set_xlabel("Max Rank"); ax2.set_ylabel("Czas [us]")
            ax2.grid(True); ax2.legend()

            plt.savefig(f"wykres_sensitivity_k{int(k)}.png", dpi=300)
            print(f"Wygenerowano: wykres_sensitivity_k{int(k)}.png")
            plt.close()
    except Exception as e:
        print(f"Błąd przy results_sensitivity.txt: {e}")

if __name__ == "__main__":
    draw_hmatrix()
    plot_performance()
    plot_sensitivity_full()