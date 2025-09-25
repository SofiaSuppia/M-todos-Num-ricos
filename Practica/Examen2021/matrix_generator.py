n = 100  # System size
"""
--> El índice i del for representa la fila actual que estás construyendo.
--> El [] dentro del row indica las columnas de la fila
--> Esto sucede porque primero definis una fila, luego se escribe esa fila, y vuelve a empezar el bucle for escribiendo otra fila debajo
--> La lista row representa el contenido de esa fila.
--> El write con \n asegura que cada fila esté en una línea diferente dentro del archivo. """
with open("data.dat", "w") as f:
    for i in range(n):
        row = [0]*(n+1)            
                                      

        # Diagonal principal
        row[i] = -2

        # Diagonal superior
        if i + 1 < n:
            row[i+1] = 1
        
        # Diagonal inferior
        if i - 1 >= 0:
            row[i-1] = 1
        
        """ # Condiciones de frontera
        if i == 0:
            row = [0]*(n+1)
            row[0] = 1
        if i == n-1:
            row = [0]*(n+1)
            row[n-1] = 1

        # Términos independientes
        if i == 0 or i == n-1:
            row[n] = 1
        else:
            row[n] = 4  """
            
            


        # Escribir fila con término independiente
        f.write(" ".join(f"{val:3d}" for val in row) + "\n")

    """ for i in range(n+1):
        if i == 0:
            row[n+1] = 1 """