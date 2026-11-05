import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import laplace

raw_data = np.loadtxt("lab5_hist_raw.txt")
img_raw = raw_data.reshape(128, 128)
lam = 0.624865
h = 1.0 / (1.0 + 8.0 * lam)

img_unc = img_raw.copy()
for _ in range(50):

    grad = (img_unc - img_raw) - lam * laplace(img_unc, mode='reflect')
    img_unc = img_unc - h * grad

unc_data = img_unc.flatten()

raw_out = np.sum((raw_data < 0.0) | (raw_data > 1.0))
unc_out = np.sum((unc_data < 0.0) | (unc_data > 1.0))

plt.figure(figsize=(10, 5))
plt.hist(raw_data, bins=120, alpha=0.55, label=f"сире спостереження (поза [0, 1]: {raw_out})", color="#1f77b4")
plt.hist(unc_data, bins=120, alpha=0.65, label=f"безумовний розв'язок (поза [0, 1]: {unc_out})", color="#ff7f0e")

plt.axvline(0.0, color="#d62728", ls="--", linewidth=1.5, label="межі допустимої області")
plt.axvline(1.0, color="#d62728", ls="--", linewidth=1.5)

plt.title("Гістограма яскравостей (як у практикумі, крок 5.1)", fontsize=13)
plt.xlabel("яскравість пікселя", fontsize=11)
plt.ylabel("кількість пікселів", fontsize=11)
plt.legend(fontsize=10)
plt.grid(True, which="both", color="#e6e6e6", linestyle='-', alpha=0.7)
for spine in plt.gca().spines.values(): spine.set_color('#cccccc')

plt.tight_layout()
plt.savefig("lab5_histogram.png", dpi=200, bbox_inches="tight")
plt.close()


if os.path.exists("lab5_pg.txt") and os.path.exists("lab5_fista.txt"):
    pg_data = np.loadtxt("lab5_pg.txt")
    fista_data = np.loadtxt("lab5_fista.txt")

    if pg_data.ndim == 1: pg_data = pg_data.reshape(1, -1)
    if fista_data.ndim == 1: fista_data = fista_data.reshape(1, -1)

    plt.figure(figsize=(10, 6))

    plt.semilogy(pg_data[:, 0], pg_data[:, 1], label="Проєкція градієнта", color="#1f77b4", linewidth=2.5)
    plt.semilogy(fista_data[:, 0], fista_data[:, 1], label="Проєкція градієнта з прискоренням (FISTA)", color="#d62728", linestyle="--", linewidth=2.5)

    plt.title("Збіжність проєкційних методів", fontsize=14)
    plt.xlabel("Ітерація k", fontsize=12)
    plt.ylabel("E(I_k)", fontsize=12)
    plt.legend(fontsize=11)
    
    plt.grid(True, which="major", color="#e6e6e6", linestyle='-')
    plt.grid(True, which="minor", color="#f9f9f9", linestyle='-')
    for spine in plt.gca().spines.values(): spine.set_color('#cccccc')

    plt.tight_layout()
    plt.savefig("lab5_convergence.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("Графіки успішно згенеровані!")
else:
    print("Помилка: Файли lab5_pg.txt або lab5_fista.txt не знайдені.")