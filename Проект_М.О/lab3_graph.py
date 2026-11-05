import numpy as np
import matplotlib.pyplot as plt

nm = np.loadtxt("lab3_nm.txt")
pw = np.loadtxt("lab3_powell.txt")

plt.figure(figsize=(9, 6))

plt.plot(range(len(nm)), nm[:, 2], label="Нелдер—Мід", color='#1f77b4', lw=2.5)
plt.plot(range(len(pw)), pw[:, 2], label="Пауелл", color='#d62728', linestyle='--', lw=2.5)

plt.yscale('log') 
plt.xlabel("Ітерація k", fontsize=12)
plt.title("Спадання цільової функції", fontsize=14)
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend(fontsize=12)

plt.tight_layout()
plt.savefig("lab3_descent_plot.png", dpi=150)

import matplotlib.colors as colors

try:
    data = np.loadtxt("lab3_grid.txt")
    ss = np.unique(data[:, 0])
    sr = np.unique(data[:, 1])
    SS, SR = np.meshgrid(ss, sr)
    Z = data[:, 2].reshape(len(ss), len(sr)).T

    nm = np.loadtxt("lab3_nm.txt")
    pw = np.loadtxt("lab3_powell.txt")
except FileNotFoundError:
    print("Помилка: Не знайдено файли даних.")
    exit()

plt.figure(figsize=(10, 6))

norm = colors.LogNorm(vmin=Z.min(), vmax=Z.max())
levels = np.logspace(np.log10(Z.min()), np.log10(Z.max()), 30)

contour = plt.contourf(SS, SR, Z, levels=levels, norm=norm, alpha=0.9, cmap='viridis')
cbar = plt.colorbar(contour)
cbar.set_label("f(θ)")
cbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.5f'))

plt.plot(nm[:, 0], nm[:, 1], color='#e32636', marker='.', ms=5, lw=1.5, label="Нелдер-Мід")
plt.plot(pw[:, 0], pw[:, 1], color='#32cd32', marker='.', ms=8, lw=1.5, label="Пауелл")

plt.scatter([nm[-1, 0]], [nm[-1, 1]], marker="*", color='#7cfc00', s=300, edgecolor='green', zorder=5)

plt.xlabel("σ_s — просторовий масштаб")
plt.ylabel("σ_r — фотометричний поріг")
plt.title("Лінії рівня f(σ_s, σ_r) і траєкторії методів (Варіант 7)")
plt.legend(loc='upper right', framealpha=0.5)

plt.tight_layout()
plt.savefig('lab3_final_plot.png', dpi=150)
