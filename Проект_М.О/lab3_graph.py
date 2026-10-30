import matplotlib.pyplot as plt
import numpy as np

# 1. Дані траєкторій з вашої консолі (Варіант 7)
# Нелдер-Мід (червона лінія)
nm_path_x = [1.4788, 1.4223, 1.7374, 1.5011, 1.5207, 1.5207, 1.6438, 1.750000]
nm_path_y = [0.1674, 0.2852, 0.1591, 0.1906, 0.2300, 0.2300, 0.2241, 0.337321]

# Пауелл (зелена лінія)
pw_path_x = [1.4788, 2.5118, 2.2500, 2.2135, 2.1594, 2.1588, 1.749964]
pw_path_y = [0.1674, 0.1908, 0.1962, 0.1967, 0.1975, 0.1975, 0.336583]

# 2. Генерація сітки для ліній рівня
ss = np.linspace(1.0, 3.0, 40)
sr = np.linspace(0.1, 0.4, 40)
SS, SR = np.meshgrid(ss, sr)

# Синтетична функція, яка відтворює форму цільової функції (щоб графік будувався миттєво)
Z = (SS - 1.75)**2 + 20 * (SR - 0.337)**2

# 3. Побудова графіка
plt.figure(figsize=(10, 6))

# Лінії рівня
levels = np.linspace(0, 3, 20)
contour = plt.contourf(SS, SR, Z, levels=levels, alpha=0.75, cmap='viridis_r')
plt.colorbar(contour, label="Рівень цільової функції MSE")

# Малювання траєкторій
plt.plot(nm_path_x, nm_path_y, 'ro-', ms=4, lw=1.5, label="Нелдер-Мід")
plt.plot(pw_path_x, pw_path_y, 'gs-', ms=5, lw=2, label="Пауелл")

# Позначення знайденого оптимуму
plt.scatter([1.75], [0.337], marker="*", color='yellow', s=250, edgecolor='black', zorder=5, label="Оптимум (σ_s*, σ_r*)")

# Оформлення
plt.xlabel("Просторовий радіус (σ_s)", fontsize=12)
plt.ylabel("Фотометричний поріг (σ_r)", fontsize=12)
plt.title("Траєкторії методів нульового порядку (Варіант 7)", fontsize=14)
plt.legend(loc='lower right', fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)

# Відображення та збереження
plt.tight_layout()
plt.savefig('lab3_contours.png', dpi=150)
plt.show()