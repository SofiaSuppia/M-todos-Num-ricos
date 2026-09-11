#Paquetes
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# SECCION: LECTURA DE ARCHIVO
# Se lee una tabla de datos desde un archivo de texto
# (mismo criterio que los programas en C++: el archivo debe
# estar en la misma carpeta que este script).
# Formato esperado: una linea por punto ->  xi yi
# La lectura es identica a la original: np.loadtxt y columnas 0 y 1.
# ============================================================
NOMBRE_ARCHIVO = 'data.txt'

def leer_tabla(nombre_archivo=NOMBRE_ARCHIVO):
    data = np.loadtxt(nombre_archivo)  # leo el archivo
    x = data[:, 0]  # primer columna
    y = data[:, 1]  # segunda columna
    return x, y

# ============================================================
# SECCION: GRAFICO DE LA TABLA
# ============================================================
def graficar_tabla(x, y, titulo='plot', xlabel='x', ylabel='y',
                   xmin=0, xmax=20, ymin=0, ymax=20,
                   color='r', marcador='o'):
    plt.plot(x, y, marker=marcador, color=color, linestyle='-')  # puntos y linea
    plt.xlim(xmin, xmax)  # limites del eje x
    plt.ylim(ymin, ymax)  # limites del eje y
    plt.title(titulo)  # titulo del grafico
    plt.xlabel(xlabel)  # nombre del eje x
    plt.ylabel(ylabel)  # nombre del eje y
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.show()

# Uso
x, y = leer_tabla()
graficar_tabla(x, y)

"""
https://matplotlib.org/stable/api/markers_api.html ---> distintos tipos de puntos
https://matplotlib.org/stable/gallery/color/named_colors.html ----> colores
https://matplotlib.org/stable/api/_gen/as matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle --> tipo de linea
"""