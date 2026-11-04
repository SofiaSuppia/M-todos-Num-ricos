# Packages
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Read data
BASE_DIR = Path(__file__).resolve().parent
RUTA_DATOS = BASE_DIR / "output" / "results.txt"

# Si el programa que genera los datos los dejo en la misma carpeta, tambien
# se acepta ese caso.
if not RUTA_DATOS.exists():
    RUTA_DATOS = BASE_DIR / "results.txt"

data = np.loadtxt(RUTA_DATOS) # Read txt file
x = data[:,0] # First column
y = data[:,1] # Second column

# Plot
plt.plot(x,y, marker='o', color='r') # Plot the graph
plt.xlim(min(x)-0.5, max(x)+0.5) # X-axis limits based on data
plt.ylim(min(y)-0.5, max(y)+0.5) # Y-axis limits based on data
plt.title('Plot') # Title of the graph
plt.xlabel('X') # X-axis label
plt.ylabel('Y') # Y-axis label
ruta_salida = Path(__file__).resolve().parent / "grafico_integracion.png"

plt.savefig(ruta_salida, dpi=300, bbox_inches="tight")
print(f"Gráfico guardado en: {ruta_salida}")

plt.show()

"""
https://matplotlib.org/stable/api/markers_api.html ---> Different types of points
https://matplotlib.org/stable/gallery/color/named_colors.html ----> Colors
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle --> Line styles
"""