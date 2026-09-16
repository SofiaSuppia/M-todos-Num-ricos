#Paquetes
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# SECCION: LECTURA DE ARCHIVO
# Se lee una tabla de datos desde un archivo de texto.
# Formato esperado: una linea por punto ->  xi yi
#
# Para que se pueda ejecutar desde CUALQUIER carpeta, la busqueda
# es así:
#   1) Si se pasa un argumento (ruta del archivo) se usa ese:
#        python3 GrafDatosTabla.py /ruta/archivo.txt
#   2) Si no se pasa nada, se busca 'data.txt' en la carpeta donde
#      esta este script (asi es portable: se copia script + data.txt
#      a la carpeta nueva y funciona igual).
# ============================================================
NOMBRE_ARCHIVO = 'data.txt'

def ruta_data():
    if len(sys.argv) > 1:
        return sys.argv[1]  # ruta pasada por el usuario
    # carpeta donde esta el script, no la carpeta desde donde se ejecuta
    carpeta_script = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(carpeta_script, NOMBRE_ARCHIVO)

def leer_tabla(nombre_archivo=None):
    if nombre_archivo is None:
        nombre_archivo = ruta_data()
    data = np.loadtxt(nombre_archivo)  # leo el archivo
    x = data[:, 0]  # primera columna
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