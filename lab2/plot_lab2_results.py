import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

# --- Konfiguracja ---
# Lista plików do przetworzenia. 
# Klucz (str) to nazwa pliku CSV.
# Wartość (str) to prefiks, który pojawi się w tytułach wykresów i nazwach plików PNG.
FILES_TO_PLOT = {
    "benchmark_gauss.csv": "Gauss",
    "benchmark_lu.csv": "LU",
    "benchmark_invert.csv": "Invert"
}
# --------------------


def generate_plots(csv_filename, algorithm_prefix):
    """
    Wczytuje jeden plik CSV i generuje dla niego 4 wykresy.
    """
    print(f"\n--- Przetwarzanie pliku: {csv_filename} (Prefiks: {algorithm_prefix}) ---")
    
    try:
        df = pd.read_csv(csv_filename)
    except FileNotFoundError:
        print(f"BŁĄD: Nie znaleziono pliku {csv_filename}. Pomijanie.")
        return
    except pd.errors.EmptyDataError:
        print(f"BŁĄD: Plik {csv_filename} jest pusty. Pomijanie.")
        return

    # --- 1. Przygotowanie danych ---
    
    # Konwertuje "300x300" na liczbę całkowitą 300
    try:
        df['Size'] = df['Dimensions'].str.split('x').str[0].astype(int)
    except KeyError:
        print(f"BŁĄD: W pliku {csv_filename} brakuje kolumny 'Dimensions'. Sprawdź nagłówek CSV.")
        return

    # Oblicza FLOPS/s (Operations / (Duration_ms / 1000))
    # Używamy apply, aby bezpiecznie obsłużyć dzielenie przez 0
    df['FLOPS'] = df.apply(
        lambda row: row['Operations'] / (row['Duration_ms'] / 1000.0) if row['Duration_ms'] > 0 else 0,
        axis=1
    )
    
    # Filtruje dane na dwa algorytmy
    df_binet = df[df['Algorithm'] == f"{algorithm_prefix}_Binet"]
    df_strassen = df[df['Algorithm'] == f"{algorithm_prefix}_Strassen"]

    if df_binet.empty or df_strassen.empty:
        print(f"OSTRZEŻENIE: Nie znaleziono danych dla Binet lub Strassen w {csv_filename}. Plik może być niekompletny.")
        
    # --- 2. Wykres Czasu ---
    plt.figure(figsize=(10, 8))
    plt.title(f"Czas działania w zależności od wielkości macierzy ({algorithm_prefix})")
    plt.xlabel("Rozmiar macierzy (n)")
    plt.ylabel("Czas wykonania (ms) [skala log]")
    plt.plot(df_binet['Size'], df_binet['Duration_ms'], '.-', label=f"{algorithm_prefix}_Binet")
    plt.plot(df_strassen['Size'], df_strassen['Duration_ms'], '.-', label=f"{algorithm_prefix}_Strassen")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{algorithm_prefix}_duration.png")
    plt.show()

    # --- 3. Wykres Operacji ---
    plt.figure(figsize=(10, 8))
    plt.title(f"Liczba operacji w zależności od wielkości macierzy ({algorithm_prefix})")
    plt.xlabel("Rozmiar macierzy (n)")
    plt.ylabel("Liczba operacji [skala log]")
    plt.plot(df_binet['Size'], df_binet['Operations'], '.-', label=f"{algorithm_prefix}_Binet")
    plt.plot(df_strassen['Size'], df_strassen['Operations'], '.-', label=f"{algorithm_prefix}_Strassen")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{algorithm_prefix}_operations.png")
    plt.show()

    # --- 4. Wykres Pamięci ---
    plt.figure(figsize=(10, 8))
    plt.title(f"Zużycie pamięci w zależności od wielkości macierzy ({algorithm_prefix})")
    plt.xlabel("Rozmiar macierzy (n)")
    plt.ylabel("Pamięć (KB) [skala liniowa]")
    plt.plot(df_binet['Size'], df_binet['Memory_kb'], '.-', label=f"{algorithm_prefix}_Binet")
    plt.plot(df_strassen['Size'], df_strassen['Memory_kb'], '.-', label=f"{algorithm_prefix}_Strassen")
    # Usunięto yscale('log') z powodu dużej liczby zer w danych
    plt.yscale('linear')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{algorithm_prefix}_memory.png")
    plt.show()

    # --- 5. Wykres FLOPS/s ---
    plt.figure(figsize=(10, 8))
    plt.title(f"Wydajność (FLOPS/s) w zależności od wielkości macierzy ({algorithm_prefix})")
    plt.xlabel("Rozmiar macierzy (n)")
    plt.ylabel("FLOPS/s (Operacje / sek) [skala log]")
    plt.plot(df_binet['Size'], df_binet['FLOPS'], '.-', label=f"{algorithm_prefix}_Binet")
    plt.plot(df_strassen['Size'], df_strassen['FLOPS'], '.-', label=f"{algorithm_prefix}_Strassen")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{algorithm_prefix}_flops.png")
    plt.show()

    print(f"Ukończono generowanie 4 wykresów dla: {algorithm_prefix}")


# --- Główna pętla wykonawcza ---
if __name__ == "__main__":
    for filename, prefix in FILES_TO_PLOT.items():
        generate_plots(filename, prefix)
    
    print("\nWszystkie wykresy zostały wygenerowane.")