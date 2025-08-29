import numpy as np
import matplotlib.pyplot as plt

# Paramètres physiques
a0 = 0.547  # nm
Omega = a0**3 / 4  # nm³
b = 0.116  # nm

# Fonctions pour le rayon Rn selon la valeur de n
def rayon(n):
    R = np.zeros_like(n, dtype=float)
    
    # Pour les bulles (n <= -1)
    n_bulles = n <= -1
    R[n_bulles] = ((3 * Omega) / (4 * np.pi))**(1/3) * np.abs(n[n_bulles])**(1/3)
    
    # Pour un soluté en position interstitielle (n = 0)
    n_interstitiel = n == 0
    R[n_interstitiel] = (Omega / (np.pi * b))**0.5
    
    # Pour les boucles (n >= 1)
    n_boucles = n >= 1
    R[n_boucles] = (Omega / (np.pi * b))**0.5 * (n[n_boucles])**0.5
    
    return R, n_bulles, n_interstitiel, n_boucles

# Domaine de n
n_vals = np.arange(-50, 50)
R_vals, n_bulles, n_interstitiel, n_boucles = rayon(n_vals)
D_vals = 2 * R_vals  # diamètre 2R

# Tracé
plt.figure(figsize=(8, 6))

# Tracer chaque catégorie avec une couleur différente
plt.plot(n_vals[n_bulles], D_vals[n_bulles], 'o', color='orange', label='Bulles (n ≤ -1)', markersize=4)
plt.plot(n_vals[n_interstitiel], D_vals[n_interstitiel], 'o', color='green', label='Soluté en interstitiel (n = 0)', markersize=6)
plt.plot(n_vals[n_boucles], D_vals[n_boucles], 'o', color='magenta', label='Boucles (n ≥ 1)', markersize=4)

plt.axvline(0, linestyle='--', color='black')

# Labels
plt.xlabel("n")
plt.ylabel("2R(n) (nm)")
plt.title("Diamètre des amas en fonction de leur contenu en auto-défauts")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Affichage
plt.show()
