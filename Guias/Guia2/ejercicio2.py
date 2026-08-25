import numpy as np
import matplotlib.pyplot as plt

def punto_fijo(g, x0, tol=1e-12, max_iter=1000):
    """
    Aproxima el punto fijo de g(x) dado un valor inicial x0.
    """
    x = x0
    for i in range(1, max_iter + 1):
        try:
            x_new = g(x)
            if abs(x_new - x) < tol:
                return x_new, i
            x = x_new
        except OverflowError:
            return None, i
    return None, max_iter

# Definimos las funciones
def g1(x): return x**5 - 3*x**3 - 2*x**2 + 2
def g2(x): return np.cos(np.sin(x))
def g3(n): return n - np.cos(n)

# Evaluamos cada función con una semilla (x0) adecuada
p1, iter1 = punto_fijo(g1, 0.0) # Semilla x0 = 0
p2, iter2 = punto_fijo(g2, 0.0) # Semilla x0 = 0
p3, iter3 = punto_fijo(g3, -1.0) # Semilla n0 = -1

print(f"Punto fijo de g1: {p1:.12f} (encontrado en {iter1} iteraciones)")
print(f"Punto fijo de g2: {p2:.12f} (encontrado en {iter2} iteraciones)")
print(f"Punto fijo de g3: {p3:.12f} (encontrado en {iter3} iteraciones)")

# --- SECCIÓN DE GRÁFICAS ---
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Gráfica 1
x1 = np.linspace(-2, 2.2, 400)
axs[0].plot(x1, g1(x1), label='g(x) = x^5 - 3x^3 - 2x^2 + 2', color='blue')
axs[0].plot(x1, x1, 'k--', label='y = x')
axs[0].scatter([p1], [p1], color='red', zorder=5)
axs[0].set_title('Función 1')
axs[0].legend()
axs[0].grid(True)

# Gráfica 2
x2 = np.linspace(0, 1.5, 400)
axs[1].plot(x2, g2(x2), label='g(x) = cos(sin(x))', color='green')
axs[1].plot(x2, x2, 'k--', label='y = x')
axs[1].scatter([p2], [p2], color='red', zorder=5)
axs[1].set_title('Función 2')
axs[1].legend()
axs[1].grid(True)

# Gráfica 3
x3 = np.linspace(-3, 6, 400)
axs[2].plot(x3, g3(x3), label='g(n) = n - cos(n)', color='purple')
axs[2].plot(x3, x3, 'k--', label='y = x')
axs[2].scatter([p3], [p3], color='red', zorder=5)
axs[2].set_title('Función 3')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.show()