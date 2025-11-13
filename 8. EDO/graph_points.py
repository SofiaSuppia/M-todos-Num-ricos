# ==============================
# Comparación: Euler vs. Solución Exacta
# ==============================

import numpy as np
import matplotlib.pyplot as plt

# --- Leer datos de results.txt (generado por el código C) ---
# Se asume que cada fila tiene la estructura: X[i]  Y[i]
try:
    data = np.loadtxt('results.txt')
except FileNotFoundError:
    print("ERROR: El archivo 'results.txt' no se encontró.")
    print("Asegúrate de ejecutar el programa C/C++ primero para generar los datos.")
    exit()

# Separar las columnas en arrays
x_euler = data[:, 0]  # Primera columna (valores de x del método de Euler)
y_euler = data[:, 1]  # Segunda columna (valores de y del método de Euler)

# --- Calcular la solución analítica exacta y(x) = 1 / (x² + 1) ---
# Generamos 200 puntos entre el mínimo y el máximo de x_euler para una curva suave
x_exact = np.linspace(min(x_euler), max(x_euler), 200)
# La solución exacta se calcula para cada punto de x_exact
y_exact = 1.0 / (x_exact * x_exact + 1.0)

# --- Graficar ambas curvas ---
plt.figure(figsize=(10, 6))

# Puntos del método de Euler (Rojo, círculos, línea discontinua)
plt.plot(x_euler, y_euler, 'ro--', label="Método de Euler", markersize=6)

# Curva de la solución exacta (Azul, línea sólida)
plt.plot(x_exact, y_exact, 'b-', label="Solución Exacta", linewidth=2.5)

# --- Formato del Gráfico ---
plt.title("Comparación entre el Método de Euler y la Solución Exacta", fontsize=15, fontweight='bold')
plt.xlabel("Variable Independiente x", fontsize=12)
plt.ylabel("Variable Dependiente y(x)", fontsize=12)
plt.legend(loc='upper right') # Muestra la leyenda en la esquina superior derecha
plt.grid(True, linestyle='--', alpha=0.7)

# Límites de los ejes extendidos ligeramente para mejor visualización
try:
    plt.xlim(min(x_euler) - 0.2, max(x_euler) + 0.2)
    plt.ylim(min(min(y_euler), min(y_exact)) - 0.1, max(max(y_euler), max(y_exact)) + 0.1)
except ValueError:
    # Esto maneja el caso donde los arrays están vacíos (si el archivo no se pudo leer)
    pass


# Mostrar gráfico
plt.show()

"""
Referencias para personalización:
https://matplotlib.org/stable/api/markers_api.html ---> Diferentes tipos de puntos
https://matplotlib.org/stable/gallery/color/named_colors.html ----> Nombres de colores
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle --> Estilos de línea
"""