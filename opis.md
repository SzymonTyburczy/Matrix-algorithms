# Sprawozdanie — Lab: Rekurencyjna kompresja obrazów z wykorzystaniem SVD

**Autorzy:** Marek Swakoń, Szymon Tyburczy

## 1. Cel i zakres

Celem laboratorium było zaimplementowanie i przeanalizowanie algorytmu **rekurencyjnej kompresji obrazów** (Image Compression) wykorzystującego strukturę drzewa czwórkowego (Quadtree) oraz rozkład według wartości osobliwych (SVD).

Zadanie polegało na:
1.  Wczytaniu obrazu i rozbiciu go na niezależne kanały RGB.
2.  Zaimplementowaniu funkcji rekurencyjnej dzielącej obraz na mniejsze bloki w zależności od lokalnej złożoności danych (błędu aproksymacji).
3.  Zastosowaniu **częściowego SVD (Truncated SVD)** jako metody kompresji w liściach drzewa.
4.  Wizualizacji efektów kompresji dla różnych parametrów sterujących: maksymalnego rzędu ($b$) oraz progu błędu ($\delta$).

## 2. Podstawy teoretyczne i Algorytm

Algorytm opiera się na adaptacyjnym podziale obrazu. Obszary "gładkie" (np. niebo, tło) są aproksymowane dużymi blokami o niskim rzędzie, natomiast obszary bogate w detale (krawędzie, tekstury) są rekurencyjnie dzielone na mniejsze fragmenty.

### Kryterium podziału (zgodne ze schematem blokowym)

Decyzja o podziale bloku podejmowana jest na podstawie analizy wartości osobliwych macierzy, zgodnie z warunkiem:
$$D(r+1, r+1) < \epsilon$$

Gdzie w naszej implementacji:
- $r$ (lub $b$) – maksymalny rząd (rank) kompresji.
- $\epsilon$ (lub $\delta$) – próg błędu (threshold).
- $D(r+1, r+1)$ – to $(b+1)$-sza wartość osobliwa ($\sigma_{b+1}$).

### Pseudokod: `RecursiveCompressSVD(Matrix A, b, delta)`

1.  **Oblicz SVD:** `[U, S, Vt] = SVD(A)`
2.  **Warunek Dopuszczalności:**
    * **JEŚLI** `S[b] < delta` (czyli odrzucona wartość jest mała):
        * Błąd jest akceptowalny.
        * **STOP:** Utwórz liść drzewa. Zapisz macierz używając tylko $k=b$ wartości osobliwych.
    * **W PRZECIWNYM RAZIE:**
        * Błąd jest zbyt duży (utracilibyśmy istotne informacje).
        * **REKURENCJA:** Podziel blok $A$ na 4 ćwiartki ($NW, NE, SW, SE$) i wywołaj algorytm dla każdej z nich.
3.  **Złóż wynik:** Zwróć węzeł drzewa zawierający albo dane skompresowane (liść), albo listę dzieci (węzeł).

## 3. Implementacja (Python)

Do realizacji zadania wykorzystano język Python oraz biblioteki `numpy` (algebra liniowa) i `matplotlib` (wizualizacja). Poniżej przedstawiono kluczowe fragmenty kodu.

### 3.1. Struktura Drzewa
```python
class QuadNode:
    def __init__(self, x, y, width, height, compressed_data=None, children=None):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        # compressed_data: Tuple (U, S, Vt) - tylko dla liści
        self.compressed_data = compressed_data
        # children: Lista 4 obiektów QuadNode - tylko dla węzłów
        self.children = children



def recursive_compress_svd(matrix, delta, b, x, y, min_size=4):
    h, w = matrix.shape
    
    # 1. Oblicz SVD (Truncated - full_matrices=False dla wydajności)
    try:
        U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    except np.linalg.LinAlgError:
        # Obsługa błędów numerycznych dla pustych/zdegenerowanych macierzy
        return QuadNode(x, y, w, h, compressed_data=(np.zeros((h,1)), [0], np.zeros((1,w))))

    # 2. Logika decyzyjna
    should_split = False
    
    # Jeśli naturalny rząd macierzy jest mniejszy niż 'b', nie tracimy danych
    if len(S) <= b:
        should_split = False
        k = len(S)
    else:
        # Sprawdzamy kluczowy warunek: czy (b+1)-sza wartość jest poniżej progu?
        if S[b] < delta:
            # Wartość jest mała -> Błąd jest akceptowalny -> LIŚĆ
            should_split = False
            k = b
        else:
            # Wartość jest duża -> Błąd za duży -> DZIELIMY
            should_split = True

    # Zabezpieczenie przed nieskończoną rekurencją (rozmiar minimalny)
    if w <= min_size or h <= min_size:
        should_split = False
        k = min(len(S), b)

    # 3. Wykonanie akcji
    if not should_split:
        # Tworzenie liścia z danymi skompresowanymi
        return QuadNode(x, y, w, h, compressed_data=(U[:, :k], S[:k], Vt[:k, :]))
    else:
        # Podział na 4 ćwiartki (Quadtree split)
        half_w, half_h = w // 2, h // 2
        children = [
            recursive_compress_svd(matrix[:half_h, :half_w], delta, b, x, y, min_size),
            recursive_compress_svd(matrix[:half_h, half_w:], delta, b, x + half_w, y, min_size),
            recursive_compress_svd(matrix[half_h:, :half_w], delta, b, x, y + half_h, min_size),
            recursive_compress_svd(matrix[half_h:, half_w:], delta, b, x + half_w, y + half_h, min_size)
        ]
        return QuadNode(x, y, w, h, children=children)
```


## 4. Metodologia testów

W celu weryfikacji poprawności algorytmu przeprowadzono serię eksperymentów na obrazach rzeczywistych (fotografie) oraz syntetycznych (gradienty).

**Parametry eksperymentalne:**
* **$b$ (Rank):** Testowano zakres od 1 (bardzo agresywna kompresja, każda macierz to iloczyn wektorów kolumnowego i wierszowego) do 10.
* **$\delta$ (Delta):** Próg błędu sterujący głębokością rekurencji. Testowano wartości z zakresu $10.0$ (wysoka wierność) do $1000.0$ (widoczna blokowość).
* **`min_size`:** Minimalny rozmiar bloku. Zmieniono domyślną wartość z 4px na 32px w testach wizualizacyjnych, aby wyraźnie uwidocznić strukturę bloków (efekt "szachownicy").

## 5. Wyniki eksperymentów

### 5.1. Adaptacyjność struktury Quadtree

Zaobserwowano wyraźną korelację między strukturą obrazu a głębokością podziału:
* **Obszary o niskiej częstotliwości (tło, niebo, rozmycia):** Algorytm poprawnie identyfikuje te obszary jako "łatwe" dla SVD. Warunek $\sigma_{b+1} < \delta$ jest spełniony bardzo wcześnie, co skutkuje pozostawieniem dużych bloków (np. 128x128 px).
* **Obszary o wysokiej częstotliwości (krawędzie, tekstury):** Pojedynczy rozkład SVD niskiego rzędu nie jest w stanie odwzorować ostrych krawędzi (duży błąd aproksymacji). Algorytm wymusza rekurencyjny podział, schodząc aż do poziomu `min_size`, aby zminimalizować błąd lokalny.

### 5.2. Wizualizacja efektów (Przypadek testowy)

Aby potwierdzić działanie algorytmu, przeprowadzono test "ekstremalny" z następującymi parametrami:
* `b = 1` (Maksymalna kompresja wewnątrz bloku).
* `delta = 1000.0` (Wysoki próg tolerancji błędu).
* `min_size = 32` (Wymuszenie dużych najmniejszych bloków).

**Obserwacje:**
Na obrazie wynikowym uzyskano widoczną strukturę blokową.
1.  W obszarach jednolitych widoczne są duże kwadraty (powyżej 32px), reprezentowane jako proste gradienty liniowe (cecha SVD rzędu 1).
2.  W obszarach krawędziowych (kontury obiektów) algorytm dokonał podziału na najmniejsze dozwolone bloki (32x32 px), próbując dopasować gradienty lokalnie.
3.  Dowodzi to, że mechanizm decyzyjny algorytmu funkcjonuje poprawnie.

*[Miejsce na wklejenie zrzutu ekranu z wygenerowanego obrazu]*

## 6. Wnioski

1.  **Wyższość nad globalnym SVD:** Zastosowanie struktury Quadtree eliminuje główną wadę globalnego SVD, jaką jest rozmywanie lokalnych detali na rzecz globalnego dopasowania. Dzięki rekurencji, detale są izolowane w mniejszych pod-macierzach.
2.  **Efektywność kompresji:** Metoda jest bardzo skuteczna dla obrazów zawierających duże, gładkie powierzchnie. Dla obrazów typu "szum" (np. trawa, piasek), struktura drzewa staje się bardzo gęsta, co może prowadzić do wzrostu rozmiaru danych zamiast ich redukcji (narzut na przechowywanie struktury drzewa).
3.  **Stabilność numeryczna:** Algorytm jest stabilny. Wykorzystanie SVD (nawet w wersji *truncated*) gwarantuje optymalną aproksymację macierzy w sensie normy Frobeniusa dla zadanego rzędu.
4.  **Znaczenie parametru Delta:** Parametr $\delta$ pełni kluczową rolę w balansowaniu między jakością a stopniem kompresji. Jego dynamiczne dostosowywanie (np. zależne od wariancji w bloku) mogłoby stanowić kierunek dalszego rozwoju algorytmu.
