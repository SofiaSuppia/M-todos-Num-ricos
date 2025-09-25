# Genera data.dat para el sistema de banda en el Problema 4
n = 50  # Tamaño del sistema (número de ecuaciones)
# Abre el archivo data.dat en modo escritura ('w')
with open("data.dat", "w") as f:
    # Itera a través de cada fila de la matriz, desde i = 0 hasta 49
    for i in range(n):
        # Crea una lista de 50 ceros para representar la fila actual de coeficientes
        fila = [0]*n
        
        # Coeficientes de la banda (estableciendo valores no nulos)
        
        # Coeficiente de la diagonal i-2 (dos posiciones a la izquierda de la diagonal principal)
        # Solo se aplica si i-2 >= 0 (es decir, a partir de la fila 2)
        if i-2 >= 0:  
            fila[i-2] = 1
            
        # Coeficiente de la diagonal i-1 (una posición a la izquierda)
        # Solo se aplica si i-1 >= 0 (es decir, a partir de la fila 1)
        if i-1 >= 0:  
            fila[i-1] = -2
            
        # Coeficiente de la diagonal principal (posición i)
        fila[i] = 12
        
        # Coeficiente de la diagonal i+1 (una posición a la derecha)
        # Solo se aplica si i+1 < n (es decir, hasta la fila 48)
        if i+1 < n:  
            fila[i+1] = -2
            
        # Coeficiente de la diagonal i+2 (dos posiciones a la derecha)
        # Solo se aplica si i+2 < n (es decir, hasta la fila 47)
        if i+2 < n:  
            fila[i+2] = 1
            
        # Escribe la fila de coeficientes separada por espacios, y añade el término independiente (5)
        f.write(" ".join(f"{valor:3d}" for valor in fila) + " 5\n")