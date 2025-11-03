from __future__ import annotations

"""
Plotowanie wyników benchmarków algorytmów mnożenia macierzy.

Funkcja `plot_benchmarks` wczytuje pliki CSV (AI, Binet, Strassen) i rysuje
trzy wykresy (po jednym dla metryki): czas wykonania [ms], zużycie pamięci [KB]
oraz liczba operacji — każdy wykres zawiera wszystkie algorytmy na jednej figurze.
Zapisuje wykresy do katalogu `plots` w folderze `lab1`.

Wymagania: pandas, matplotlib
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

import re


def _ensure_packages():
	"""Leniwa kontrola obecności wymaganych paczek z przyjaznym komunikatem.

	Nie instalujemy automatycznie; jeśli brakuje, wyrzucamy czytelną informację.
	"""

	try:
		import pandas as _pd  # noqa: F401
		import matplotlib as _mpl  # noqa: F401
	except Exception as exc:  # pragma: no cover - komunikat runtime
		raise RuntimeError(
			"Brak wymaganych pakietów. Zainstaluj: pip install pandas matplotlib"
		) from exc


def _parse_size(dimensions: str) -> int:
	"""Ekstrahuje rozmiar N z pola tekstowego Dimensions.

	Obsługa formatów:
	- "NxN" (np. "100x100") → N
	- złożone, np. "(4x5) * (5x5)" → maksymalna liczba w napisie (tu: 5)
	- dowolne inne: wybiera maksymalną liczbę całkowitą obecną w stringu
	"""

	nums = [int(x) for x in re.findall(r"\d+", dimensions)]
	if not nums:
		raise ValueError(f"Nie można sparsować rozmiaru z Dimensions='{dimensions}'")
	return max(nums)


def _load_csv(path: Path, algorithm_label: str | None = None):
	import pandas as pd

	df = pd.read_csv(path)
	# Ujednolicenie nazw
	rename_map = {
		"Dimensions": "dimensions",
		"Algorithm": "algorithm",
		"Operations": "operations",
		"Duration_ms": "duration_ms",
		"Memory_kb": "memory_kb",
		"Passed": "passed",
	}
	df = df.rename(columns=rename_map)

	# Dodanie kolumny size wyliczonej z dimensions
	df["size"] = df["dimensions"].apply(_parse_size)

	# Ewentualne nadpisanie etykiety algorytmu
	if algorithm_label is not None:
		df["algorithm"] = algorithm_label

	# Wybór tylko potrzebnych kolumn
	keep = ["size", "algorithm", "operations", "duration_ms", "memory_kb"]
	df = df[keep]
	return df


def plot_benchmarks(
	ai_csv: Path | str = Path(__file__).with_name("benchmark_ai.csv"),
	binet_csv: Path | str = Path(__file__).with_name("benchmark_binet.csv"),
	strassen_csv: Path | str = Path(__file__).with_name("benchmark_strassen.csv"),
	output_dir: Path | str = Path(__file__).parent / "plots",
	show: bool = False,
) -> Tuple[Path, Path, Path]:
	"""Wczytuje wyniki i generuje po jednym wykresie dla każdej metryki.

	Zwraca krotkę ścieżek do utworzonych plików PNG:
	(mem_png, time_png, ops_png)
	"""

	_ensure_packages()
	import pandas as pd
	import matplotlib.pyplot as plt

	ai_csv = Path(ai_csv)
	binet_csv = Path(binet_csv)
	strassen_csv = Path(strassen_csv)
	output_dir = Path(output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)

	# Wczytanie i połączenie danych
	df_ai = _load_csv(ai_csv, algorithm_label="ai")
	df_binet = _load_csv(binet_csv, algorithm_label="binet")
	df_strassen = _load_csv(strassen_csv, algorithm_label="strassen")
	df = pd.concat([df_ai, df_binet, df_strassen], ignore_index=True)

	# Porządkowanie: liczby całkowite, sortowanie, filtr outlierów (opcjonalnie)
	df = df.assign(
		size=df["size"].astype(int),
		operations=df["operations"].astype(int),
		duration_ms=pd.to_numeric(df["duration_ms"], errors="coerce"),
		memory_kb=pd.to_numeric(df["memory_kb"], errors="coerce").fillna(0),
	)

	# Funkcja pomocnicza do rysowania jednego wykresu
	def _plot_metric(ax, metric: str, ylabel: str, logy: bool = False):
		for algo, g in df.groupby("algorithm"):
			g = g.sort_values("size")
			ax.plot(g["size"], g[metric], marker="o", label=algo)
		ax.set_xlabel("Rozmiar N")
		ax.set_ylabel(ylabel)
		ax.grid(True, which="both", linestyle=":", alpha=0.6)
		if logy:
			ax.set_yscale("log")
		ax.legend(title="Algorytm")

	# Osobne, wysokiej rozdzielczości wykresy dla każdej metryki
	def _single(metric: str, ylabel: str, filename: str, logy: bool = False) -> Path:
		fig2, ax2 = plt.subplots(figsize=(10, 5))
		_plot_metric(ax2, metric, ylabel, logy=logy)
		fig2.suptitle(f"{ylabel} w funkcji rozmiaru N")
		out = output_dir / filename
		fig2.savefig(out, dpi=180)
		if show:
			plt.show()
		plt.close(fig2)
		return out

	# Generujemy dokładnie trzy wykresy: pamięć, czas, operacje (FLOPs)
	out_mem = _single("memory_kb", "Pamięć [KB]", "memory_kb.png")
	out_time = _single("duration_ms", "Czas [ms]", "time_ms.png")
	out_ops = _single("operations", "Liczba operacji", "operations.png", logy=True)

	# Zwracamy tylko trzy ścieżki (bez figury 3-w-1)
	return out_mem, out_time, out_ops


if __name__ == "__main__":  # prosty CLI
	try:
		paths = plot_benchmarks(show=False)
		print("Zapisano wykresy:")
		for p in paths:
			print("-", p)
	except Exception as e:
		print("Błąd podczas generowania wykresów:", e)

