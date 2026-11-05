import glob
import numpy as np
import matplotlib.pyplot as plt

files_map = {
    "gd_0": ("Градієнтний спуск (h = 1/L)", '#1f77b4', '-'),           
    "armijo_0": ("Градієнтний спуск (Армійо)", '#d62728', '--'),       
    "hb_17": ("Важка кулька (оптимальне β)", '#2ca02c', '-'),          
    "hb_90": ("Важка кулька β = 0.9 (завелике)", '#ff7f0e', '--'),     
    "nag_0": ("Нестеров NAG", '#9467bd', '-'),                         
    "nag_42": ("Нестеров NAG (сильно опукла)", '#8c564b', '--')        
}


loaded_data = {}
global_min_E = float('inf')

for key in files_map.keys():
    matches = glob.glob(f"lab4_{key}*.txt")
    if matches:
        data = np.loadtxt(matches[0])
        if len(data) > 0:
            if data.ndim == 1: data = data.reshape(1, -1)
            loaded_data[key] = data
            min_in_file = np.min(data[:, 1])
            if min_in_file < global_min_E:
                global_min_E = min_in_file

E_star = global_min_E

fig1, ax1 = plt.subplots(figsize=(10, 6))
fig2, ax2 = plt.subplots(figsize=(10, 6))

for key, (label, color, ls) in files_map.items():
    if key not in loaded_data:
        continue
        
    data = loaded_data[key]
    k = data[:, 0]
    E = data[:, 1]
    dist = data[:, 2] / 1e6  

    E_diff = E - E_star

    valid_idx = E_diff > 1e-13 
    
    ax1.plot(k[valid_idx], E_diff[valid_idx], label=label, color=color, linestyle=ls, linewidth=2.5)
    
    valid_dist = dist > 1e-10
    ax2.plot(k[valid_dist], dist[valid_dist], label=label, color=color, linestyle=ls, linewidth=2.5)


ax1.set_yscale('log')
ax1.set_xlabel("Ітерація k", fontsize=12)
ax1.set_ylabel("E(I_k) − f*", fontsize=12)
ax1.set_title("Збіжність за значенням функції", fontsize=14)
ax1.legend(fontsize=11)
ax1.grid(True, which="both", ls="-", alpha=0.2)

for spine in ax1.spines.values(): spine.set_color('#cccccc')

fig1.tight_layout()
fig1.savefig("lab4_convergence_f.png", dpi=150, bbox_inches='tight')

ax2.set_yscale('log')
ax2.set_xlabel("Ітерація k", fontsize=12)
ax2.set_ylabel("|| x_k − x* ||", fontsize=12)
ax2.set_title("Відстань до точного розв'язку", fontsize=14)
ax2.legend(fontsize=11)
ax2.grid(True, which="both", ls="-", alpha=0.2)
for spine in ax2.spines.values(): spine.set_color('#cccccc')

fig2.tight_layout()
fig2.savefig("lab4_convergence_dist.png", dpi=150, bbox_inches='tight')
