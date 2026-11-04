# ============================================================================
# PROBLEMA 2 - SEGUNDO PARCIAL 2024
# ============================================================================
# Lee resultados_problema2.csv generado por problema2.c++:
#
#   x, y_exacta, y_ab2, y_rk4, error_ab2, error_rk4
#
# Guarda la figura como grafico_problema2.png.
# ============================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    datos = np.loadtxt("resultados_problema2.csv", delimiter=",", skiprows=1)
except FileNotFoundError:
    print("ERROR: no se encontro resultados_problema2.csv.")
    print("Ejecute primero el programa problema2.")
    raise SystemExit(1)

x = datos[:, 0]
y_exacta = datos[:, 1]
y_ab2 = datos[:, 2]
y_rk4 = datos[:, 3]

plt.figure(figsize=(9, 6))
plt.plot(x, y_exacta, "b-", linewidth=2.5, label="Solucion exacta")
plt.plot(x, y_ab2, "ro--", markersize=5, label="Adams-Bashforth 2 pasos")
plt.plot(x, y_rk4, "gs-.", markersize=5, label="RK4")

plt.title("Problema 2: solucion exacta, AB2 y RK4")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.7)
plt.xlim(min(x) - 0.02, max(x) + 0.02)
plt.ylim(min(y_exacta) - 0.03, max(y_exacta) + 0.03)

plt.tight_layout()
plt.savefig("grafico_problema2.png", dpi=300)
print("Figura guardada en grafico_problema2.png")
