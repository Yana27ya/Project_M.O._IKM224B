import matplotlib.pyplot as plt
import numpy as np

sigma_data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
mse_data = [0.025181, 0.025180, 0.024469, 0.018683, 0.012046, 0.009173, 
            0.008427, 0.008457, 0.008770, 0.009198, 0.009680, 0.010180, 
            0.010681, 0.011186, 0.011681, 0.012167, 0.012658, 0.013135, 0.013607]

svenn_x = [0.5663, 0.7689, 1.1741]
svenn_y = [0.009799, 0.008404, 0.010050]

plt.figure(figsize=(10, 6))

plt.plot(sigma_data, mse_data, 'b-', linewidth=2.5, label='Функція MSE f(σ)')

plt.plot(svenn_x, svenn_y, 'ro', markersize=8, label='Кроки алгоритму Свенна')
labels = ['x0', 'x1', 'x2']
for i, (x, y) in enumerate(zip(svenn_x, svenn_y)):
    plt.annotate(labels[i], (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=12)

plt.axvspan(svenn_x[0], svenn_x[2], color='green', alpha=0.15, label=f'Інтервал локалізації [{svenn_x[0]}, {svenn_x[2]}]')
plt.axvline(x=svenn_x[0], color='green', linestyle='--', alpha=0.7)
plt.axvline(x=svenn_x[2], color='green', linestyle='--', alpha=0.7)

plt.title('Залежність середньоквадратичної помилки від розмиття (Варіант 7)', fontsize=14, pad=15)
plt.xlabel('Сила розмиття Гаусса (σ)', fontsize=12)
plt.ylabel('Середньоквадратична помилка (MSE)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=11)
plt.xlim(0, 2.0)
plt.ylim(0.005, 0.028)

plt.tight_layout()
plt.show()