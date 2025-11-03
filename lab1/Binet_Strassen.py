import matplotlib.pyplot as plt
import pandas as pd



df1 = pd.read_csv("matrix_multiplication_results_BINET.csv")
df2 = pd.read_csv("matrix_multiplication_results_STRASSEN.csv")
sizes_Binet = df1['Size']
times_Binet = df1['Duration_ms']
sizes_STRASSEN = df2['Size']
times_STRASSEN = df2['Duration_ms']

plt.figure(figsize=(10, 8))
plt.title("Czas działania w zależności od wielkości macierzy")
plt.xlabel("Rozmiary macierzy n x n")
plt.ylabel("Czas wykonania")
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.plot(sizes_Binet, times_Binet ,marker = "o")
plt.plot(sizes_STRASSEN, times_STRASSEN ,marker = "o")
plt.legend(["BINET","STRASSEN"])
plt.savefig("Duration.png")
plt.show()




operations_Binet = df1['Operations']
operations_STRASSEN = df2['Operations']

plt.figure(figsize=(10, 8))
plt.title("liczba operacji zmienno-przecinkowych w zależności od wielkości macierzy")
plt.xlabel("Rozmiary macierzy")
plt.ylabel("liczba operacji zmienno-przecinkowych")
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.plot(sizes_Binet, operations_Binet ,marker = "o")
plt.plot(sizes_STRASSEN, operations_STRASSEN ,marker = "o")
plt.legend(["BINET","STRASSEN"])
plt.savefig("Operation.png")
plt.show()



memory_Binet = df1['Memory_kb']
memory_STRASSEN = df2['Memory_kb']

plt.figure(figsize=(10, 8))
plt.title("Zużycie pamięci w zależności od wielkości macierzy")
plt.xlabel("Rozmiary macierzy")
plt.ylabel("zużycie pamięci")
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.plot(sizes_Binet, memory_Binet ,marker = "o")
plt.plot(sizes_STRASSEN, memory_STRASSEN ,marker = "o")
plt.legend(["BINET","STRASSEN"])
plt.savefig("Memory.png")
plt.show()

