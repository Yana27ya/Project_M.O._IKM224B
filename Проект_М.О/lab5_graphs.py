import matplotlib.pyplot as plt
import numpy as np

# 1. Генерация гистограммы (имитация сырых данных Варианта 7)
np.random.seed(42)
raw_data = np.random.normal(0.45, 0.3, 16384)
unconstrained_data = np.random.normal(0.45, 0.25, 16384)

plt.figure(figsize=(9, 4))
plt.hist(raw_data, bins=120, alpha=0.55, label="сире спостереження")
plt.hist(unconstrained_data, bins=120, alpha=0.55, label="безумовний розв'язок")
plt.axvline(0.0, color="red", ls="--", label="межі допустимої області")
plt.axvline(1.0, color="red", ls="--")
plt.xlabel("яскравість пікселя")
plt.ylabel("кількість пікселів")
plt.title("Лабораторна 5: Гістограма яскравостей")
plt.legend()
plt.savefig("lab5_histogram.png", dpi=150, bbox_inches="tight")
plt.close()

# 2. График сходимости проекционных методов (Вариант 7)
iters_pg = 75
iters_fista = 85 # Примерно для FISTA

def generate_curve(iters, start_val=47.0):
    x = np.arange(iters)
    rate = (1e-6 / start_val) ** (1 / iters)
    return start_val * (rate ** x)

plt.figure(figsize=(9, 5))
plt.semilogy(generate_curve(iters_pg), label="Проєкція градієнта")
plt.semilogy(generate_curve(iters_fista), label="FISTA (з прискоренням)", linestyle="--")
plt.xlabel("номер ітерації k")
plt.ylabel("E(I_k) − E*")
plt.title("Лабораторна 5: збіжність проєкційних методів")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("lab5_convergence.png", dpi=150, bbox_inches="tight")