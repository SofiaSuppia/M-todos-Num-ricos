# ==============================
# Gráfico del factor de convergencia Q
# ==============================

import sys
import numpy as np
import matplotlib.pyplot as plt

# El programa C/C++ pasa el archivo como argumento:
#   python3 graph_convergence.py convergence_euler.txt
# Si no se pasa argumento, se usa convergence_rk4.txt por defecto.
FILE_NAME = sys.argv[1] if len(sys.argv) > 1 else 'convergence_rk4.txt'

try:
    # Los archivos de convergencia del código escalar tienen dos columnas: x, Q
    data = np.loadtxt(FILE_NAME)
except FileNotFoundError:
    print(f"ERROR: El archivo '{FILE_NAME}' no fue encontrado.")
    print("Asegúrate de ejecutar el programa C/C++ que genera los datos de convergencia primero.")
    exit()
except Exception as error:
    print(f"ERROR al leer '{FILE_NAME}': {error}")
    exit()

if data.ndim == 1:
    data = data.reshape(1, -1)

# Separar las columnas
x_convergence = data[:, 0]
y_convergence = data[:, 1]

# Excluir el primer punto, que suele ser Q_0 = 0 o no informativo
x_convergence = x_convergence[1:]
y_convergence = y_convergence[1:]

# Obtener un nombre legible a partir del archivo
metodo = FILE_NAME.replace('convergence_', '').replace('.txt', '')
metodo = metodo.replace('_', ' ').title()

# --- Graficar los puntos ---
plt.figure(figsize=(10, 6))
plt.plot(x_convergence, y_convergence, 'ro-', label=f"Método {metodo}", markersize=6)

plt.title(f"Factor de convergencia ({metodo})", fontsize=15, fontweight='bold')
plt.xlabel("Variable independiente x", fontsize=12)
plt.ylabel("Factor de convergencia Q[i]", fontsize=12)
plt.legend(loc='upper right')
plt.grid(True, linestyle='--', alpha=0.7)

try:
    plt.xlim(min(x_convergence) - 0.2, max(x_convergence) + 0.2)
    plt.ylim(min(y_convergence) - 0.2, max(y_convergence) + 0.2)
except ValueError:
    print("Advertencia: no hay suficientes puntos para establecer los límites del gráfico.")

plt.show()

"""
Referencias para personalización:
https://matplotlib.org/stable/api/markers_api.html
https://matplotlib.org/stable/gallery/color/named_colors.html
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html
"""
