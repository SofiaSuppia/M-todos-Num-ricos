# ============================================================
# GRAFICADOR DE FUNCIONES (inspeccion de raices)
# Permite graficar una o varias funciones f(x) y marcar sus raices.
#
# Uso (una sola funcion, igual que el original):
#   graficar_matplot(f, "cos(sin(x))", valores_iniciales=[0.1, 2])
#
# Uso (varias funciones):
#   graficar_matplot([
#       {"funcion": f1, "etiqueta": "cos(sin(x))", "color": "blue",
#        "valores_iniciales": [0.1, 2]},
#       {"funcion": f2, "etiqueta": "sin(x)", "color": "red",
#        "valores_iniciales": [0, 3.1]}
#   ])
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def graficar_matplot(funciones, etiqueta_f="", valores_iniciales=None, x_min=-5, x_max=5):
    # Si me pasan una sola funcion (forma clasica), se convierte en lista
    if callable(funciones):
        funciones = [{"funcion": funciones,
                      "etiqueta": etiqueta_f,
                      "color": "blue",
                      "valores_iniciales": valores_iniciales}]

    # Dominio comun para todas las curvas
    x = np.linspace(x_min, x_max, 400)

    plt.figure(figsize=(8, 5))

    for fun in funciones:
        f_py = fun["funcion"]
        etiqueta = fun["etiqueta"]
        color = fun.get("color", "blue")
        iniciales = fun.get("valores_iniciales", [])

        # Calcular raices con los valores iniciales dados
        raices = list(set(np.round([fsolve(f_py, x0)[0] for x0 in iniciales], 4)))

        # Graficar la curva
        plt.plot(x, f_py(x), label=f"f(x) = {etiqueta}", color=color, linewidth=2)

        # Marcar las raices de esta curva
        for r in raices:
            plt.plot(r, 0, 'o', color=color, markersize=8, label=f"Raiz de {etiqueta} en x = {r}")

    # Ejes y rejilla
    plt.axhline(0, color="gray", linewidth=1.2)
    plt.axvline(0, color="gray", linewidth=0.8, linestyle="--")
    plt.grid(True, linestyle=":", alpha=0.6)

    plt.title("Inspección de raíces")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()
    plt.show()

# Uso con una sola funcion (igual que el original)
f_py = lambda x: np.cos(np.sin(x))
graficar_matplot(f_py, "cos(sin(x))", valores_iniciales=[0.1, 2])

# Uso con varias funciones (descomentar para probar)
# f1 = lambda x: np.cos(np.sin(x))
# f2 = lambda x: np.sin(x)
# graficar_matplot([
#     {"funcion": f1, "etiqueta": "cos(sin(x))", "color": "blue", "valores_iniciales": [0.1, 2]},
#     {"funcion": f2, "etiqueta": "sin(x)", "color": "red", "valores_iniciales": [0, 3.1]}
# ])