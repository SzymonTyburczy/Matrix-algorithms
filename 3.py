import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os


def save_and_get_disk_size(np_array, temp_path, quality):
    try:
        img = Image.fromarray(np_array.astype(np.uint8))
        img.save(temp_path, format="JPEG", quality=quality)
        disk_size_bytes = os.path.getsize(temp_path)
        return disk_size_bytes

    except Exception as e:
        print(f"Błąd podczas zapisu/mierzenia: {e}")
        return 0

    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def load_and_prep_image(url=None, filename=None, size=(512, 512)):
    loaded = False
    img = None
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

class QuadNode:
    def __init__(self, x, y, width, height, compressed_data=None, children=None):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.compressed_data = compressed_data
        self.children = children


def recursive_compress_svd(matrix, delta, b, x, y, min_size=4):
    h, w = matrix.shape
    try:
        U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    except np.linalg.LinAlgError:
        return QuadNode(x, y, w, h, compressed_data=(np.zeros((h, 1)), [0], np.zeros((1, w))))

    should_split = False
    if len(S) <= b:
        k = len(S)
        should_split = False
    else:

        if S[b] < delta:
            k = b
            should_split = False
        else:
            should_split = True
    if w <= min_size or h <= min_size:
        should_split = False
        k = min(len(S), b)
    if not should_split:
        U_k = U[:, :k]
        S_k = S[:k]
        Vt_k = Vt[:k, :]
        return QuadNode(x, y, w, h, compressed_data=(U_k, S_k, Vt_k))
    else:
        half_w = w // 2
        half_h = h // 2

        children = []
        children.append(recursive_compress_svd(matrix[:half_h, :half_w], delta, b, x, y, min_size))  # TL
        children.append(recursive_compress_svd(matrix[:half_h, half_w:], delta, b, x + half_w, y, min_size))  # TR
        children.append(recursive_compress_svd(matrix[half_h:, :half_w], delta, b, x, y + half_h, min_size))  # BL
        children.append(
            recursive_compress_svd(matrix[half_h:, half_w:], delta, b, x + half_w, y + half_h, min_size))  # BR

        return QuadNode(x, y, w, h, children=children)


def decompress_channel(node, output_matrix):
    if node.children is None:
        U, S, Vt = node.compressed_data
        if len(S) > 0:
            block = U @ np.diag(S) @ Vt
        else:
            block = np.zeros((node.height, node.width))
        output_matrix[node.y: node.y + node.height, node.x: node.x + node.width] = block
    else:
        for child in node.children:
            decompress_channel(child, output_matrix)


def plot_compressed_image(original_path, original_img, nodes_R, nodes_G, nodes_B, info):
    h, w, _ = original_img.shape
    Rec_R, Rec_G, Rec_B = np.zeros((h, w)), np.zeros((h, w)), np.zeros((h, w))

    decompress_channel(nodes_R, Rec_R)
    decompress_channel(nodes_G, Rec_G)
    decompress_channel(nodes_B, Rec_B)

    Rec_Img = np.stack([Rec_R, Rec_G, Rec_B], axis=2)
    Rec_Img = np.clip(Rec_Img, 0, 255).astype(np.uint8)

    temp_file_path = "temp_compressed_svd.jpg"
    compressed_disk_size_mb = save_and_get_disk_size(Rec_Img, temp_file_path, quality=100) / (1024 * 1024)
    compressed_disk_text = f"Rozmiar JPG (j=95): {compressed_disk_size_mb:.2f} MB"
    try:
        original_size_bytes = os.path.getsize(original_path)
        original_size_mb = original_size_bytes / (1024 * 1024)
        original_size_text = f"Waga pliku na dysku: {original_size_mb:.2f} MB"

    except (FileNotFoundError, TypeError):
        original_size_text = "Nieznana waga (plik nie istnieje lub ścieżka jest niepoprawna)"

    compressed_size_bytes = Rec_Img.nbytes
    compressed_size_mb = compressed_size_bytes / (1024 * 1024)
    compressed_size_text = f"Rozmiar w pamięci (tablica): {compressed_disk_text} MB"

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.imshow(original_img)
    plt.title(f"Oryginał\n({original_size_text})")
    plt.axis('off')
    plt.subplot(1, 2, 2)
    plt.imshow(Rec_Img)
    plt.title(f"Skompresowany\n{info}\n({compressed_disk_text})")
    plt.axis('off')
    plt.show()



def Compress_Each_Channel(R, G, B, min_size, DELTA, B_RANK):
    root_R = recursive_compress_svd(R, delta=DELTA, b=B_RANK, x=0, y=0, min_size=min_size)
    root_G = recursive_compress_svd(G, delta=DELTA, b=B_RANK, x=0, y=0, min_size=min_size)
    root_B = recursive_compress_svd(B, delta=DELTA, b=B_RANK, x=0, y=0, min_size=min_size)
    return root_R, root_G, root_B

def do(filepath, MIN_SIZE, DELTA, B_RANK):
    original_img, R, G, B = load_and_prep_image(filename=filepath, size=(512, 512))
    print(f"Start kompresji algebraicznej (b={B_RANK}, delta={DELTA})...")
    root_R, root_G, root_B = Compress_Each_Channel(R, G, B, MIN_SIZE, DELTA, B_RANK)
    print("Gotowe. Rysowanie...")
    plot_compressed_image(filepath,original_img, root_R, root_G, root_B, f"b={B_RANK}, delta={DELTA}")



def all():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # file_path1 = os.path.join(current_dir, "1.jpg")
    # file_path2 = os.path.join(current_dir, "2.jpg")
    # file_path3 = os.path.join(current_dir, "3.jpg")
    # file_path4 = os.path.join(current_dir, "4.jpg")
    # file_path6 = os.path.join(current_dir, "6.jpg")
    # do(file_path1, MIN_SIZE)
    # do(file_path2, MIN_SIZE)
    # do(file_path3, MIN_SIZE)
    # do(file_path4, MIN_SIZE)
    # do(file_path6, MIN_SIZE)
    file_path5 = os.path.join(current_dir, "5.jpg")
    do(file_path5, 20, 10, 20)
    do(file_path5, 40, 400, 20)
    do(file_path5, 16, 1600, 20)



def plot_singular_values(S_channels, labels):
    plt.figure(figsize=(10, 6))
    colors = ['r', 'g', 'b']

    for i, (S, label, col) in enumerate(zip(S_channels, labels, colors)):
        plt.plot(S, label=f'Kanał {label}', color=col)

        plt.scatter(0, S[0], c=col, zorder=5)
        y_pos = S[0] * (0.3 ** i)
        plt.text(10, y_pos, rf' $\sigma_1 ({label})$={int(S[0])}', color=col, fontweight='bold')

    plt.yscale('log')
    plt.title("Wartości osobliwe (Singular Values) dla całej bitmapy")
    plt.xlabel("Indeks k")
    plt.ylabel("Wartość osobliwa sigma_k (log)")
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.show()

def run_report_experiments(filepath):
    img_arr, R, G, B = load_and_prep_image(filename=filepath, size=(512, 512))
    U_r, S_r, Vt_r = np.linalg.svd(R, full_matrices=False)
    U_g, S_g, Vt_g = np.linalg.svd(G, full_matrices=False)
    U_b, S_b, Vt_b = np.linalg.svd(B, full_matrices=False)

    plot_singular_values([S_r, S_g, S_b], ['R', 'G', 'B'])

    sigma_1 = S_r[0]
    sigma_last = S_r[-1]
    sigma_mid = S_r[len(S_r) // 2]

    print(f"Referencyjne Sigmy (z kanału R): s1={sigma_1:.2f}, s_mid={sigma_mid:.2f}, s_last={sigma_last:.2f}")

    cases = [
        (1, "sigma_1 (Max)", sigma_1),
        (1, "sigma_last (Min)", sigma_last),
        (1, "sigma_mid (Srodek)", sigma_mid),
        (4, "sigma_1 (Max)", sigma_1),
        (4, "sigma_last (Min)", sigma_last),
        (4, "sigma_mid (Srodek)", sigma_mid)
    ]

    for r_val, delta_name, delta_val in cases:
        print(f"\n--- Eksperyment: r={r_val}, delta={delta_name} ({delta_val:.2f}) ---")

        root_R = recursive_compress_svd(R, delta=delta_val, b=r_val, x=0, y=0, min_size=4)
        root_G = recursive_compress_svd(G, delta=delta_val, b=r_val, x=0, y=0, min_size=4)
        root_B = recursive_compress_svd(B, delta=delta_val, b=r_val, x=0, y=0, min_size=4)

        h, w = R.shape
        Rec_R, Rec_G, Rec_B = np.zeros((h, w)), np.zeros((h, w)), np.zeros((h, w))
        decompress_channel(root_R, Rec_R)
        decompress_channel(root_G, Rec_G)
        decompress_channel(root_B, Rec_B)

        Rec_Img = np.stack([Rec_R, Rec_G, Rec_B], axis=2)
        Rec_Img = np.clip(Rec_Img, 0, 255).astype(np.uint8)

        plt.figure(figsize=(15, 4))
        plt.suptitle(f"r={r_val}, delta={delta_name}")

        plt.subplot(1, 4, 1)
        plt.imshow(Rec_R, cmap='Reds')
        plt.title("Kanał R")
        plt.axis('off')

        plt.subplot(1, 4, 2)
        plt.imshow(Rec_G, cmap='Greens')
        plt.title("Kanał G")
        plt.axis('off')

        plt.subplot(1, 4, 3)
        plt.imshow(Rec_B, cmap='Blues')
        plt.title("Kanał B")
        plt.axis('off')

        plt.subplot(1, 4, 4)
        plt.imshow(Rec_Img)
        plt.title("Wynik RGB")
        plt.axis('off')

        plt.show()

    print("\n--- Najlepsza kompresja (Wybrane parametry) ---")
    my_r = 20
    my_delta = 20  # Niska delta = wysoka jakość
    print(f"Moje parametry: r={my_r}, delta={my_delta}")

    root_R = recursive_compress_svd(R, delta=my_delta, b=my_r, x=0, y=0, min_size=8)
    root_G = recursive_compress_svd(G, delta=my_delta, b=my_r, x=0, y=0, min_size=8)
    root_B = recursive_compress_svd(B, delta=my_delta, b=my_r, x=0, y=0, min_size=8)

    plot_compressed_image(filepath, img_arr, root_R, root_G, root_B, f"BEST: r={my_r}, d={my_delta}")


if __name__ == '__main__':
    all()
    run_report_experiments("5.jpg")
