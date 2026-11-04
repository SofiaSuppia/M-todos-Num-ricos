# ============================================================================
# PROBLEMA 3 - Comparacion entre la solucion exacta y RK4 3/8
# ============================================================================
# Lee resultados_rk4.csv generado por problema3.c++:
#
#   x, y_numerica, y_exacta, error
#
# Guarda la figura como grafico_rk4.png.
# ============================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    datos = np.loadtxt("resultados_rk4.csv", delimiter=",", skiprows=1)
except FileNotFoundError:
    print("ERROR: no se encontro resultados_rk4.csv.")
    print("Ejecute primero el programa problema3.")
    raise SystemExit(1)

x = datos[:, 0]
y_numerica = datos[:, 1]
y_exacta = datos[:, 2]

plt.figure(figsize=(9, 6))
plt.plot(x, y_exacta, "b-", linewidth=2.5, label="Solucion exacta")
plt.plot(x, y_numerica, "ro--", markersize=5, label="RK4 3/8, h = 0.1")

plt.title("Problema 3: solucion exacta vs. RK4 3/8")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.7)
plt.xlim(min(x) - 0.05, max(x) + 0.05)
plt.ylim(min(y_exacta) - 0.3, max(y_exacta) + 0.3)

plt.tight_layout()
plt.savefig("grafico_rk4.png", dpi=300)
print("Figura guardada en grafico_rk4.png")
