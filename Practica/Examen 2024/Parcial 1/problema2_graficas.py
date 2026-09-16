# ============================================================
# PARCIAL 1 - 2024 - PROBLEMA 2
# Graficos:
#  1) E(x) = |P4(x) - f(x)|  (parte c)   -> problema2_error(x).png
#  2) Datos + regresion lineal y cuadratica (d y e)
#                                          -> problema2_ajustes.png
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

# Datos de la tabla
X = np.array([1.0, 1.2, 1.5, 1.75, 2.0])
Y = np.array([-0.148, -0.040, 0.181, 0.419, 0.700])

# Funcion real del problema 1
def f(x):
    return np.log(x * x + 1) - np.sin(x)

# Polinomio interpolador de grado 4 (Lagrange) por los 5 nodos
coef = np.polyfit(X, Y, 4)          # coef [a4,...,a0]
P = np.poly1d(coef)

# ------------------------------------------------------------
# PARTE c) E(x) = |P(x) - f(x)|
# ------------------------------------------------------------
x = np.linspace(0.85, 2.15, 600)
E = np.abs(P(x) - f(x))

plt.figure(figsize=(8, 5))
plt.plot(x, E, "-", color="blue", linewidth=2,
         label="E(x) = |P4(x) - f(x)|")
plt.plot(X, np.abs(P(X) - f(X)), 'o', color="red", markersize=7,
         label="En nodos E(xi) = 0")
plt.axhline(0, color="gray", linewidth=1.2)

plt.grid(True, linestyle=":", alpha=0.6)
plt.title("Funcion error: E(x) = |P4(x) - f(x)|")
plt.xlabel("x")
plt.ylabel("E(x)")
plt.legend()
plt.savefig("problema2_error(x).png", dpi=150)
print("Grafico guardado en: problema2_error(x).png")

# ------------------------------------------------------------
# PARTES d) y e) Ajustes lineal y cuadratico
# ------------------------------------------------------------
lin = np.poly1d(np.polyfit(X, Y, 1))      # regresion lineal
cuad = np.poly1d(np.polyfit(X, Y, 2))     # regresion cuadratica

plt.figure(figsize=(8, 5))
plt.plot(X, Y, 'o', color="black", markersize=9, label="Datos de la tabla")
plt.plot(x, lin(x), "--", color="green", linewidth=2,
         label="Regresion lineal  (r = 0.9922)")
plt.plot(x, cuad(x), "-", color="purple", linewidth=2,
         label="Regresion cuadratica  (r = 0.999996)")
plt.plot(x, f(x), ":", color="blue", linewidth=2,
         label="f(x) real")
plt.axhline(0, color="gray", linewidth=1.2)

plt.grid(True, linestyle=":", alpha=0.6)
plt.title("Ajuste de los datos por regresion lineal y cuadratica")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("problema2_ajustes.png", dpi=150)
print("Grafico guardado en: problema2_ajustes.png")

if plt.get_backend().lower() != "agg":
    plt.show()