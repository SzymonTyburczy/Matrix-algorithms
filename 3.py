import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import urllib.request
from io import BytesIO
import os
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

class QuadNode:
    def __init__(self, x, y, width, height, compressed_data=None, children=None):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.compressed_data = compressed_data  # (U, S, Vt)
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



def Compress_Each_Channel(R, G, B):
    root_R = recursive_compress_svd(R, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)
    root_G = recursive_compress_svd(G, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)
    root_B = recursive_compress_svd(B, delta=DELTA, b=B_RANK, x=0, y=0, min_size=32)
    return root_R, root_G, root_B

def do(filepath):
    original_img, R, G, B = load_and_prep_image(filename=filepath, size=(512, 512))
    print(f"Start kompresji algebraicznej (b={B_RANK}, delta={DELTA})...")
    root_R, root_G, root_B = Compress_Each_Channel(R, G, B)
    print("Gotowe. Rysowanie...")
    plot_compressed_image(original_img, root_R, root_G, root_B, f"b={B_RANK}, delta={DELTA}")

B_RANK = 1
DELTA = 60.0
MIN_SIZE = 2

current_dir = os.path.dirname(os.path.abspath(__file__))
file_path1 = os.path.join(current_dir, "1.jpg")
file_path2 = os.path.join(current_dir, "2.jpg")
file_path3 = os.path.join(current_dir, "3.jpg")
file_path4 = os.path.join(current_dir, "4.jpg")

do(file_path1)
do(file_path2)
do(file_path3)
do(file_path4)


