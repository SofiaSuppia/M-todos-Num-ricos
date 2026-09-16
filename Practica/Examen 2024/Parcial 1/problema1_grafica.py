# ============================================================
# PARCIAL 1 - 2024 - PROBLEMA 1
# Plot de f(x) = ln(x^2 + 1) - sin(x)  (parte a)
# Guarda el grafico en "problema1_f(x).png"
# ============================================================
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.log(x * x + 1) - np.sin(x)

alpha = 1.26203265474

x = np.linspace(-2.5, 2.5, 800)

plt.figure(figsize=(8, 5))
plt.plot(x, f(x), label="f(x) = ln(x^2 + 1) - sin(x)", color="blue", linewidth=2)
plt.plot(alpha, 0, 'o', color="red", markersize=8,
         label=f"Raiz en x = {alpha:.4f}")

plt.axhline(0, color="gray", linewidth=1.2)
plt.axvline(alpha, color="red", linewidth=0.8, linestyle="--", alpha=0.6)

# Marcar el intervalo [1, 2]
plt.axvspan(1.0, 2.0, color="orange", alpha=0.15,
            label="Intervalo [1.0, 2.0]")

plt.grid(True, linestyle=":", alpha=0.6)
plt.title("f(x) = ln(x^2 + 1) - sin(x)")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()

# Guardar el grafico (parte a: adjuntar por aula virtual)
plt.savefig("problema1_f(x).png", dpi=150)
print("Grafico guardado en: problema1_f(x).png")

if plt.get_backend().lower() != "agg":
    plt.show()