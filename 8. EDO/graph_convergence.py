# ==============================
# 📘 Gráfico del Factor de Convergencia (Método de Euler)
# ==============================

import numpy as np
import matplotlib.pyplot as plt

# --- Leer datos de convergence_rk4.txt (generado por el código C) ---
# Cada fila se asume que contiene: X[i]  Factor_Q[i]
FILE_NAME = 'convergence_rk4.txt' 

try:
    data = np.loadtxt(FILE_NAME) 
except FileNotFoundError:
    print(f"ERROR: El archivo '{FILE_NAME}' no fue encontrado.")
    print("Asegúrate de ejecutar el programa C/C++ que genera los datos de convergencia primero.")
    exit()

# Separar las columnas en arrays
x_convergence = data[:, 0]  # Primera columna (valores de x)
y_convergence = data[:, 1]  # Segunda columna (Factor de Convergencia Q)

# --- Excluir el primer punto (índice 0) ---
# A menudo, el primer punto en el factor de convergencia (Q_0) no es útil o es 0.
x_convergence = x_convergence[1:]
y_convergence = y_convergence[1:]

# --- Graficar los puntos ---
plt.figure(figsize=(10, 6))

# Puntos del factor de convergencia (Rojo, círculos, línea discontinua)
# 'ro--' -> Red (rojo), Marker Circle (o), Line Dashed (--)
plt.plot(x_convergence, y_convergence, 'ro-', 
         label="Método de Euler", markersize=6)

# --- Formato del Gráfico ---
plt.title("Factor de Convergencia (Método de Euler)", 
          fontsize=15, fontweight='bold')
plt.xlabel("Variable Independiente x", fontsize=12)
plt.ylabel("Factor de Convergencia Q[i]", fontsize=12)
plt.legend(loc='upper right')
plt.grid(True, linestyle='--', alpha=0.7)

# Límites de los ejes extendidos ligeramente para mejor visualización
try:
    plt.xlim(min(x_convergence) - 0.2, max(x_convergence) + 0.2)
    plt.ylim(min(y_convergence) - 0.2, max(y_convergence) + 0.2)
except ValueError:
    print("Advertencia: No hay suficientes puntos para establecer los límites del gráfico.")

# Mostrar gráfico
plt.show()

"""
Referencias para personalización:
https://matplotlib.org/stable/api/markers_api.html ---> Diferentes tipos de puntos
https://matplotlib.org/stable/gallery/color/named_colors.html ----> Nombres de colores
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle --> Estilos de línea
"""