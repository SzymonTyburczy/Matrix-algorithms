import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import glob
def draw_hmatrix(filenames="matrix_trees.txt"):
    print(f"Przetwarzam pliki: {filenames}")
    pliki = glob.glob("matrix_tree_k*.txt")
    pliki.sort()
    for filename in pliki:
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Błąd: Nie znaleziono pliku {filename}. Uruchom najpierw program C++.")
            return

        # Ustalenie rozmiaru macierzy na podstawie pierwszego węzła (korzenia)
        if not lines:
            return
        
        root_data = list(map(int, lines[0].split()))
        # root_data: x, y, width, height, rank, is_leaf
        matrix_width = root_data[2]
        matrix_height = root_data[3]

        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Rysowanie
        for line in lines:
            data = list(map(int, line.split()))
            x, y, w, h, rank, is_leaf = data
            
            # Logika rysowania (zgodna ze slajdami):
            # Rysujemy tylko LIŚCIE.
            # Jeśli rank > 0 -> Czarny/Szary blok (zawiera dane)
            # Jeśli rank == 0 -> Biały blok (pusty)
            
            if is_leaf:
                if rank > 0:
                    # Blok z danymi (SVD)
                    rect = patches.Rectangle((x, y), w, h, linewidth=0.5, edgecolor='none', facecolor='black')
                    ax.add_patch(rect)
                else:
                    # Blok pusty (opcjonalnie można nie rysować nic, wtedy będzie białe tło)
                    pass
            else:
                # Węzły wewnętrzne - rysujemy tylko obramowanie, żeby widzieć siatkę
                rect = patches.Rectangle((x, y), w, h, linewidth=0.2, edgecolor='gray', facecolor='none')
                ax.add_patch(rect)

        # Ustawienia wykresu
        ax.set_xlim(0, matrix_width)
        ax.set_ylim(0, matrix_height)
        ax.invert_yaxis() # Macierze mają (0,0) w lewym górnym rogu
        ax.set_aspect('equal')
        plt.title(f"Struktura H-Macierzy (N={matrix_width})")
        plt.xlabel("Kolumny")
        plt.ylabel("Wiersze")
        
        # Zapis do pliku i wyświetlenie
        output_img = f"hmatrix_visualization_{filename}.png"
        plt.savefig(output_img, dpi=300)
        print(f"Rysunek zapisano jako {output_img}")
        # plt.show() # Odkomentuj, jeśli masz środowisko graficzne

if __name__ == "__main__":
    draw_hmatrix()