import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Ejercicio: Polinomio de interpolacion de Lagrange cuadratico
# Nodos: x0=1, x1=2, x2=2.5
# f(x) = x + 2/x
# ============================================================

# Nodos y sus imagenes
x_nodos = np.array([1.0, 2.0, 2.5])
f = lambda x: x + 2.0 / x          # funcion real
y_nodos = f(x_nodos)               # f en los nodos
N = len(x_nodos) - 1               # grado del polinomio (N=2)

# --- Funcion base de Lagrange C(N,k,x) ---
#   C(N,k,x) = PROD(i!=k) (x - x[i]) / (x[k] - x[i])
def CNk(k, x):
    prod = np.ones_like(x)
    for i in range(N + 1):
        if i != k:
            prod = prod * (x - x_nodos[i]) / (x_nodos[k] - x_nodos[i])
    return prod

# --- Polinomio interpolante Pn(x) ---
#   Pn(x) = SUMA(k=0..N) f(x[k]) * C(N,k,x)
def Pn(x):
    suma = np.zeros_like(x)
    for k in range(N + 1):
        suma = suma + y_nodos[k] * CNk(k, x)
    return suma

# Dominio de grafico (cerca de los nodos)
x = np.linspace(0.5, 3.5, 500)

# ============================================================
# GRAFICO 1: Coeficientes de Lagrange CN,k(x)
# Se verifica que CN,k(xk) = 1  y  CN,k(xj) = 0 para j != k
# (en cada nodo pasa exactamente por 1 para su k y por 0 para las demas)
# ============================================================
plt.figure(figsize=(9, 5))
for k in range(N + 1):
    plt.plot(x, CNk(k, x), label=f"C$_{N},{k}$(x) = L${k}$(x)", linewidth=2)

# Se marcan los nodos
plt.plot(x_nodos, np.ones_like(x_nodos), 'ko', markersize=6)

for k in range(N + 1):
    for j in range(N + 1):
        # CNk(xj): en el nodo j, C vale 1 si j==k y 0 si j!=k
        plt.plot(x_nodos[j], CNk(k, np.array([x_nodos[j]]))[0], 'o', markersize=6,
                 color=f'C{k}')

plt.axhline(0, color="gray", linewidth=1)
plt.axhline(1, color="gray", linewidth=1, linestyle="--")
plt.grid(True, linestyle=":", alpha=0.6)
plt.title("Coeficientes de Lagrange C$_{2,k}$(x)  |  verificar C$_{2,k}$(xk)=1, C$_{2,k}$(xj)=0")
plt.xlabel("x")
plt.ylabel("C$_{2,k}$(x)")
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# GRAFICO 2: f(x) vs Pn(x)
# El polinomio pasa exactamente por los nodos (e(xk)=0)
# ============================================================
plt.figure(figsize=(9, 5))
plt.plot(x, f(x), label="f(x) = x + 2/x", color="blue", linewidth=2)
plt.plot(x, Pn(x), label="P$_2$(x) (Lagrange)", color="red", linewidth=2, linestyle="--")
plt.plot(x_nodos, y_nodos, 'ko', markersize=8, label="Nodos (x$_k$, f(x$_k$))")
plt.axhline(0, color="gray", linewidth=1)
plt.grid(True, linestyle=":", alpha=0.6)
plt.title("f(x) vs Polinomio interpolante P$_2$(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# GRAFICO 3: Error e(x) = |f(x) - Pn(x)|
# Se nota que e(xk) = 0 en todos los nodos
# ============================================================
e = np.abs(f(x) - Pn(x))
plt.figure(figsize=(9, 5))
plt.plot(x, e, label="e(x) = |f(x) - P$_2$(x)|", color="green", linewidth=2)
plt.plot(x_nodos, np.zeros_like(x_nodos), 'ro', markersize=8, label="e(x$_k$) = 0")
plt.axvline(1.5, color="gray", linestyle="--", alpha=0.7)
plt.axvline(1.2, color="gray", linestyle="--", alpha=0.7)
plt.grid(True, linestyle=":", alpha=0.6)
plt.title("Error de interpolacion e(x) = |f(x) - P$_2$(x)|")
plt.xlabel("x")
plt.ylabel("|e(x)|")
plt.legend()
plt.tight_layout()
plt.show()

# Valores pedidos por el ejercicio
for xh in [1.5, 1.2]:
    coefs = [CNk(k, np.array([xh]))[0] for k in range(N + 1)]
    print(f"x = {xh}")
    print(f"  C2,0 = {coefs[0]:.6f}  C2,1 = {coefs[1]:.6f}  C2,2 = {coefs[2]:.6f}")
    print(f"  Suma de coeficientes = {sum(coefs):.6f}")
    print(f"  P2(x) = {Pn(np.array([xh]))[0]:.6f}   f(x) = {f(xh):.6f}   error = {abs(f(xh) - Pn(np.array([xh]))[0]):.6f}")
    print()