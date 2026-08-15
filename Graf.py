import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def graficar_matplot(f_py, etiqueta_f, valores_iniciales, x_min=-5, x_max=5):
    # Calcular raíces
    raices = list(set(np.round([fsolve(f_py, x0)[0] for x0 in valores_iniciales], 4)))
    
    # Dominio y función
    x = np.linspace(x_min, x_max, 400)
    
    plt.figure(figsize=(8, 5))
    plt.plot(x, f_py(x), label=f"f(x) = {etiqueta_f}", color="blue", linewidth=2)
    
    # Ejes y rejilla
    plt.axhline(0, color="gray", linewidth=1.2)
    plt.axvline(0, color="gray", linewidth=0.8, linestyle="--")
    plt.grid(True, linestyle=":", alpha=0.6)
    
    # Graficar raíces calculadas
    for r in raices:
        plt.plot(r, 0, 'ro', markersize=8, label=f"Raíz en x = {r}")
        
    plt.title("Inspección de raíces")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()
    plt.show()

# Uso
f_py = lambda x:  x**10 - 1
graficar_matplot(f_py, "x^10 - 1", valores_iniciales=[0.1, 2])

