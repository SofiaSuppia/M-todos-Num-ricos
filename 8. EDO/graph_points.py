# ==============================
# Comparación: solución numérica vs. solución exacta
# ==============================

import numpy as np
import matplotlib.pyplot as plt

# results.txt tiene el encabezado:
# x    y_exacta    y_numerica    error_exacto
try:
    data = np.loadtxt('results.txt', skiprows=1)
except FileNotFoundError:
    print("ERROR: El archivo 'results.txt' no se encontró.")
    print("Asegúrate de ejecutar el programa C/C++ primero para generar los datos.")
    exit()
except Exception as error:
    print(f"ERROR al leer 'results.txt': {error}")
    exit()

if data.ndim == 1:
    data = data.reshape(1, -1)

# Separar las columnas
x = data[:, 0]          # x
y_exact = data[:, 1]    # solución exacta
y_num = data[:, 2]      # solución numérica

# --- Graficar ambas curvas ---
plt.figure(figsize=(10, 6))

plt.plot(x, y_num, 'ro--', label="Solución numérica", markersize=6)
plt.plot(x, y_exact, 'b-', label="Solución exacta", linewidth=2.5)

# --- Formato del gráfico ---
plt.title("Comparación entre la solución numérica y la exacta",
          fontsize=15, fontweight='bold')
plt.xlabel("Variable independiente x", fontsize=12)
plt.ylabel("Variable dependiente y(x)", fontsize=12)
plt.legend(loc='upper right')
plt.grid(True, linestyle='--', alpha=0.7)

# Límites de los ejes extendidos ligeramente para mejor visualización
try:
    plt.xlim(min(x) - 0.2, max(x) + 0.2)
    plt.ylim(min(min(y_num), min(y_exact)) - 0.1,
             max(max(y_num), max(y_exact)) + 0.1)
except ValueError:
    pass

plt.show()

"""
Referencias para personalización:
https://matplotlib.org/stable/api/markers_api.html
https://matplotlib.org/stable/gallery/color/named_colors.html
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html
"""
