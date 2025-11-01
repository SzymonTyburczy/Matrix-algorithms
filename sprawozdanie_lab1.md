# Sprawozdanie — Lab 1: Rekurencyjne mnożenie macierzy (Binét, Strassen, AI)

Autorzy: Marek Swakoń, Szymon Tyburczy

## Cel i zakres

W ramach pierwszego zestawu zadań zaimplementowano i przeanalizowano trzy algorytmy mnożenia macierzy z losowymi wartościami z przedziału otwartego (0.00000001, 1.0):

- Rekurencyjne mnożenie macierzy metodą Binét’a
- Rekurencyjne mnożenie macierzy metodą Strassena
- Mnożenie macierzy metodą AI (na podstawie artykułu w Nature)

W trakcie eksperymentów zliczano liczbę operacji zmiennoprzecinkowych (dodawanie, odejmowanie, mnożenie, dzielenie) wykonywanych podczas mnożenia.

## Generowanie danych wejściowych

- Generator macierzy: wartości losowe z przedziału (0.00000001, 1.0), rozkład jednostajny.
- Rozmiary testowe: 1, 2, 3, …, 1000 (do maksymalnego N, który dało się policzyć na stanowisku).
- Warunki brzegowe: małe N (walidacja poprawności), duże N (pomiar wydajności i zużycia pamięci).

## Pseudokod (wysoki poziom)

### 1) Binét — rekurencyjne mnożenie

1. Wejście: A, B ∈ R^{N×N}
2. Jeśli N = 1: zwróć [A[0,0]·B[0,0]]
3. Podziel A i B na cztery bloki N/2×N/2: A11, A12, A21, A22; B11, B12, B21, B22
4. Rekurencyjnie oblicz wymagane kombinacje bloków zgodne z metodą Binét’a
5. Złącz bloki w wynik C

Warunek zakończenia: ustalony pragmatyczny próg N0, poniżej którego stosujemy metodę iteracyjną O(N^3) (optymalizacja).

### 2) Strassen — rekurencyjne mnożenie

1. Wejście: A, B ∈ R^{N×N}
2. Jeśli N = 1: zwróć [A[0,0]·B[0,0]]
3. Wyznacz 7 iloczynów pośrednich (M1…M7) zgodnie z klasycznym schematem Strassena
4. Złóż C z kombinacji M1…M7
5. Zwróć C

Uwagi: Dla nieparzystych N bez użycia paddingu należy obsłużyć „ostatni wiersz/kolumnę” oddzielnie (np. rozbicie na bloki niesymetryczne lub iteracyjne domnożenie brzegu).

### 3) Metoda AI (wg Nature)

1. Wejście: A ∈ R^{m×k}, B ∈ R^{k×p}
2. Dobierz schemat mnożenia i kolejności operacji zgodny z opisem metody AI, zliczając operacje zmiennoprzecinkowe
3. Jeśli metoda ma ograniczenia wymiarów — wypisz komunikat i przerwij (bez paddingu)
4. Zwróć C ∈ R^{m×p}

---

## Najważniejsze fragmenty kodu

- Konstrukcja macierzy losowych (C++/Python)
- Implementacja rekurencyjna Binét
- Implementacja Strassena (z obsługą brzegów lub informacją o ograniczeniach)
- Implementacja AI (oraz licznik operacji)

Fragmenty wklejamy poniżej (maks. 20–30 linii każdy) z krótkim komentarzem:

```cpp
// PRZYKŁAD: licznik operacji i dodawanie macierzy
unsigned long long op_count = 0;
Matrix add(const Matrix& A, const Matrix& B) {
	Matrix C(A.size(), std::vector<double>(A[0].size()));
	for (size_t i = 0; i < A.size(); ++i)
		for (size_t j = 0; j < A[0].size(); ++j) {
			C[i][j] = A[i][j] + B[i][j];
			op_count++; // zliczamy „+”
		}
	return C;
}
```

---

## Metodologia pomiarowa

- Czas: std::chrono::high_resolution_clock (średnia z ≥3 powtórzeń na punkt)
- Operacje: ręcznie inkrementowany licznik przy każdej operacji +, −, ·, /
- Pamięć: WinAPI (GetProcessMemoryInfo / PSAPI) lub odpowiednik; raport w KB
- Walidacja: porównanie wyników z wersją iteracyjną O(N^3) przy tolerancji ϵ = 1e−9

Format wyników (CSV):

```
Size,Algorithm,Operations,Duration_ms,Memory_kb
16,Strassen,XXXX,YY.Y,ZZZZ
...
```

Pliki CSV w repozytorium:

- `matrix_multiplication_results_BINET.csv`
- `matrix_multiplication_results_STRASSEN.csv`
- `matrix_multiplication_results_AlphaTensor.csv` lub `ai_recursive_benchmark.csv`

---

## Wyniki i wykresy

Wstaw wykresy (lub podlinkuj pliki PNG) dla zakresu N = 1…N_max:

1. Czas działania (ms) vs. rozmiar macierzy (oś X)

![Czas](Duration.png)

2. Liczba operacji zmiennoprzecinkowych vs. rozmiar

![Operacje](Operation.png)

3. Zużycie pamięci (KB) vs. rozmiar

![Pamięć](Memory.png)

Tabela przykładowa (fragment):

| N   | Algorytm | Operacje [#] | Czas [ms] | Pamięć [KB] |
| --- | -------- | ------------ | --------- | ----------- |
| 16  | Binét    | …            | …         | …           |
| 16  | Strassen | …            | …         | …           |
| 16  | AI       | …            | …         | …           |

---

## Ograniczenia, błędy i obsługa przypadków brzegowych

- Brak paddingu: jeżeli algorytm wymaga specyficznych wymiarów (np. parzyste N), program powinien zakończyć pracę i wypisać czytelny komunikat (bez sztucznego dopełniania).
- Stabilność numeryczna: porównania wyników z tolerancją ϵ.
- Zużycie pamięci: opis ewentualnych pików pamięci i ich przyczyn.

---

## Szacowanie złożoności obliczeniowej

- Binét: … (wyprowadzenie/odwołanie do literatury) — eksperymentalne dopasowanie krzywej
- Strassen: O(N^{log₂7}) ≈ O(N^{2.807}); pomiar vs. teoria
- AI: opis przewidywanej złożoności i obserwacje eksperymentalne

Metoda estymacji: regresja log–log (czas/operacje vs. N), wykres trendu i współczynnik R².

---

## Walidacja poprawności

- Testy małych rozmiarów (N=2,3,5,7,…) — porównanie z implementacją iteracyjną
- Dodatkowo testy losowe dla kilku N, weryfikacja normy błędu ||C_ref − C||\_F

---

## Instrukcja uruchomienia (Windows, PowerShell)

Kompilacja (MSYS2/MinGW, C++17):

```powershell
g++ -std=c++17 -O2 -Wall -Wextra -DNOMINMAX -o matrix_Binet.exe matrix_Binet.cpp -lpsapi
g++ -std=c++17 -O2 -Wall -Wextra -DNOMINMAX -o matrix_Strassen.exe matrix_Strassen.cpp -lpsapi
g++ -std=c++17 -O2 -Wall -Wextra -DNOMINMAX -o matrix_ai.exe matrix_ai.cpp -lpsapi
```

Uruchomienie:

```powershell
./matrix_Binet.exe
./matrix_Strassen.exe
./matrix_ai.exe
```

Skrypt z wykresami (opcjonalnie):

```powershell
py ./Binet_Strassen.py
```

---

## Dyskusja i wnioski

- Porównanie osiągów (czas, operacje, pamięć) między Binét / Strassen / AI
- Kiedy który algorytm jest korzystny i dlaczego
- Wpływ ograniczeń (brak paddingu) na implementację i wyniki
- Potencjalne kierunki optymalizacji

---

## Bibliografia

1. Volker Strassen, Gaussian Elimination is Not Optimal, 1969
2. Artykuł „Nature” dot. mnożenia macierzy metodą AI — pełny cytat
3. Inne źródła (podręczniki, wykłady, wpisy blogowe) — pełne odniesienia

---

Checklist (do odhaczenia przed oddaniem):

- [ ] Pseudokod obu algorytmów rekurencyjnych
- [ ] Fragmenty kodu (najważniejsze miejsca)
- [ ] Wykresy: czas, operacje, pamięć (1…N_max)
- [ ] CSV z wynikami dołączone do repo
- [ ] Oszacowanie złożoności (teoria/eksperyment)
- [ ] Walidacja poprawności na małych N
- [ ] Opis ograniczeń i braku paddingu
