import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import urllib.request
from io import BytesIO


# --- 1. Ładowanie obrazu (Pancerne, z User-Agent) ---
def load_and_prep_image(url=None, filename=None, size=(512, 512)):
    loaded = False
    img = None

    if url:
        try:
            req = urllib.request.Request(
                url,
                data=None,
                headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
            )
            with urllib.request.urlopen(req) as response:
                img_data = response.read()
            img = Image.open(BytesIO(img_data))
            loaded = True
        except Exception as e:
            print(f"Błąd URL: {e}")

    if not loaded and filename:
        try:
            img = Image.open(filename)
            loaded = True
        except Exception:
            pass

    if not loaded:
        print("Generuję gradient zastępczy.")
        img = Image.new('RGB', size)
        pixels = img.load()
        for i in range(img.size[0]):
            for j in range(img.size[1]):
                pixels[i, j] = (i % 255, (i + j) % 255, j % 255)

    img = img.convert('RGB').resize(size)
    img_arr = np.array(img, dtype=float)
    return img_arr.astype(np.uint8), img_arr[:, :, 0], img_arr[:, :, 1], img_arr[:, :, 2]


# --- 2. Struktura Quadtree ---
class QuadNode:
    def __init__(self, x, y, width, height, compressed_data=None, children=None):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.compressed_data = compressed_data  # (U, S, Vt)
        self.children = children


# --- 3. REKURENCYJNA KOMPRESJA (Zgodna ze slajdem) ---
def recursive_compress_svd(matrix, delta, b, x, y, min_size=4):
    """
    Algorytm zgodny ze slajdem:
    b (r): ranga kompresji
    delta (epsilon): próg wartości osobliwej
    """
    h, w = matrix.shape

    # 1. Oblicz SVD
    try:
        U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    except np.linalg.LinAlgError:
        # Pusty/zły blok -> traktujemy jako zero
        return QuadNode(x, y, w, h, compressed_data=(np.zeros((h, 1)), [0], np.zeros((1, w))))

    # 2. Logika decyzyjna (SLIDE LOGIC)
    # Warunek ze zdjęcia: if D(r+1, r+1) < epsilon
    # W Pythonie (indeksowanie od 0): sprawdzamy S[b]

    should_split = False

    # Jeśli blok jest mniejszy niż 'b', SVD jest dokładne (nie ma co odrzucać)
    if len(S) <= b:
        k = len(S)
        should_split = False
    else:
        # Tu jest kluczowy warunek ze slajdu:
        # Sprawdzamy wartość, którą byśmy odrzucili (S[b])
        if S[b] < delta:
            # Wartość jest mała -> błąd jest mały -> KOMPRESUJEMY (Liść)
            k = b
            should_split = False
        else:
            # Wartość jest duża -> błąd byłby duży -> DZIELIMY (Rekurencja)
            should_split = True

    # Zabezpieczenie przed nieskończoną rekurencją na pojedynczych pikselach
    if w <= min_size or h <= min_size:
        should_split = False
        k = min(len(S), b)

    # 3. Wykonanie akcji
    if not should_split:
        # LIŚĆ (CompressMatrix)
        U_k = U[:, :k]
        S_k = S[:k]
        Vt_k = Vt[:k, :]
        return QuadNode(x, y, w, h, compressed_data=(U_k, S_k, Vt_k))
    else:
        # WĘZEŁ (CreateTree x 4)
        half_w = w // 2
        half_h = h // 2

        children = []
        children.append(recursive_compress_svd(matrix[:half_h, :half_w], delta, b, x, y, min_size))  # TL
        children.append(recursive_compress_svd(matrix[:half_h, half_w:], delta, b, x + half_w, y, min_size))  # TR
        children.append(recursive_compress_svd(matrix[half_h:, :half_w], delta, b, x, y + half_h, min_size))  # BL
        children.append(
            recursive_compress_svd(matrix[half_h:, half_w:], delta, b, x + half_w, y + half_h, min_size))  # BR

        return QuadNode(x, y, w, h, children=children)


# --- 4. Rekonstrukcja ---
def decompress_channel(node, output_matrix):
    if node.children is None:
        U, S, Vt = node.compressed_data
        # Odtworzenie bloku: A = U * diag(S) * Vt
        if len(S) > 0:
            block = U @ np.diag(S) @ Vt
        else:
            block = np.zeros((node.height, node.width))
        output_matrix[node.y: node.y + node.height, node.x: node.x + node.width] = block
    else:
        for child in node.children:
            decompress_channel(child, output_matrix)


def plot_compressed_image(original_img, nodes_R, nodes_G, nodes_B, info):
    h, w, _ = original_img.shape
    Rec_R, Rec_G, Rec_B = np.zeros((h, w)), np.zeros((h, w)), np.zeros((h, w))

    decompress_channel(nodes_R, Rec_R)
    decompress_channel(nodes_G, Rec_G)
    decompress_channel(nodes_B, Rec_B)

    Rec_Img = np.stack([Rec_R, Rec_G, Rec_B], axis=2)
    Rec_Img = np.clip(Rec_Img, 0, 255).astype(np.uint8)

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1);
    plt.imshow(original_img);
    plt.title("Oryginał");
    plt.axis('off')
    plt.subplot(1, 2, 2);
    plt.imshow(Rec_Img);
    plt.title(f"Skompresowany\n{info}");
    plt.axis('off')
    plt.show()


# --- MAIN ---

# Link do obrazu (Astronauta - dobry do testów)
# Zachód słońca - gładkie przejścia kolorów
# Przykład dla Windows
# --- MAIN ---

# Ścieżka do Twojego pliku (używamy surowego stringa 'r' przed cudzysłowem)
file_path = r"C:\Users\szymo\Desktop\kkk.jpg"

# Poprawne wywołanie funkcji:
# 1. Przekazujemy pełną ścieżkę (file_path) do parametru filename.
# 2. Odbieramy 4 zmienne (original_img oraz R, G, B).
original_img, R, G, B = load_and_prep_image(filename=file_path, size=(512, 512))

# PARAMETRY (Zgodne z teorią)
B_RANK = 1      # r: ile wartości zachowujemy w bloku (rank)
DELTA = 600.0    # epsilon: próg błędu (wartość osobliwa)
MIN_SIZE = 32
print(f"Start kompresji algebraicznej (b={B_RANK}, delta={DELTA})...")

# Uruchomienie kompresji dla każdego kanału
root_R = recursive_compress_svd(R, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)
root_G = recursive_compress_svd(G, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)
root_B = recursive_compress_svd(B, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)

print("Gotowe. Rysowanie...")
plot_compressed_image(original_img, root_R, root_G, root_B, f"b={B_RANK}, delta={DELTA}")