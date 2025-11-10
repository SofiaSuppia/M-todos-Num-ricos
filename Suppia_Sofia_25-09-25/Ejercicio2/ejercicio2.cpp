#include <cstdio>  
#include <stdlib.h> 
#include <cmath>   

#define MAX_TAMANO 50  

/**
 * @param nombre_archivo Nombre del archivo a leer
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param n Puntero al tamaño del sistema (número de ecuaciones)
 * @return true si el archivo fue leído exitosamente, false en caso contrario
 * */
bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_TAMANO+1], double b[], int* n);

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma_temp, aux;
    int se_uso_pivoteo_parcial = 0;
    
    // Definición de arreglos usando el MAX_TAMANO global
    double a[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], X[MAX_TAMANO+1]; // X es el vector solución

    // Leer arreglo desde el archivo usando la función. '
    if(!leer_archivo_matriz("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Suppia_Sofia_25-09-25\\Ejercicio2\\data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Sistema de ecuaciones original:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);      // 10 espacios totales, 6 decimales
    }
    printf("\n");
    
    // Recorrer las filas de la matriz (Eliminación Gaussiana)
    for(int i = 1; i <= n-1; i++) {
        
        printf("Paso %d: Analizando elemento diagonal a[%d][%d] = %.10lf\n", i, i, i, a[i][i]);
        
        // CRITERIO PARA PIVOTEO: Verificar si el elemento diagonal es muy pequeño
        if(fabs(a[i][i]) < 1e-5) {
            printf(">>> NECESITA PIVOTEO: |a[%d][%d]| = %.2e < 1e-5 (muy pequeño)\n", i, i, fabs(a[i][i]));
            printf("    Razones:\n");
            printf("    - División por número muy pequeño causa errores numéricos\n");
            printf("    - Puede causar overflow en los cálculos\n");
            printf("    - Reduce la estabilidad numérica del algoritmo\n");
            se_uso_pivoteo_parcial = 1;
            
            // Buscar el mejor pivote (elemento con mayor valor absoluto)
            p = i;
            printf("    Buscando mejor pivote en columna %d:\n", i);
            for(int l = i+1; l <= n; l++) {
                printf("      a[%d][%d] = %.6lf (|valor| = %.6lf)\n", 
                       l, i, a[l][i], fabs(a[l][i]));
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l; // Encontramos la fila con el elemento más grande en la columna i
                    printf("      >>> Nuevo mejor pivote encontrado en fila %d\n", l);
                }
            }
            printf("    Pivote seleccionado: a[%d][%d] = %.6lf\n", p, i, a[p][i]);
            
            // Si la fila con el pivote más grande es diferente, intercambiamos filas
            if (p != i) { 
                printf("    >>> INTERCAMBIANDO: Fila %d ↔ Fila %d\n", i, p);
                for(int m = i; m <= n; m++) {
                    aux = a[p][m];
                    a[p][m] = a[i][m];
                    a[i][m] = aux; // Intercambiar las filas p e i en la matriz A
                }
                aux = b[p];
                b[p] = b[i];
                b[i] = aux; // Intercambiar el término independiente
                printf("    Nuevo pivote después del intercambio: a[%d][%d] = %.6lf\n", i, i, a[i][i]);
            } else if (fabs(a[i][i]) < 1e-5) {
                // Si incluso después del pivoteo el pivote es cero, el sistema no tiene solución única
                printf("    ERROR: Incluso después del pivoteo, no hay pivote válido\n");
                printf("    El sistema es singular (no tiene solución única)\n");
                exit(1);
            }
        } else {
            printf(">>> NO NECESITA PIVOTEO: |a[%d][%d]| = %.6lf >= 1e-5 (suficientemente grande)\n", 
                   i, i, fabs(a[i][i]));
            printf("    El elemento diagonal es estable para la eliminación\n");
            p = i; // No hay cambio de fila
        }

        // Hacer cero los elementos debajo de la diagonal en la columna actual
        for(int j = i+1; j <= n; j++) {
            factor = a[j][i] / a[i][i]; // Calculamos el factor multiplicador (sin el signo negativo)
            
            // Recorrer las columnas en la fila j para la eliminación
            for(int k = i; k <= n; k++) {
                // Fila j = Fila j - factor * Fila i
                a[j][k] = a[j][k] - factor * a[i][k];
            }
            // Actualizar el vector de términos independientes
            b[j] = b[j] - factor * b[i];
        }
    }

    // RESUMEN DE ANÁLISIS DE PIVOTEO
    printf("\n=============== ANÁLISIS DE PIVOTEO ===============\n");
    if(se_uso_pivoteo_parcial == 1) {
        printf("✅ SE UTILIZÓ PIVOTEO PARCIAL\n");
        printf("Justificación:\n");
        printf("- Se encontraron elementos diagonales menores a 1e-5\n");
        printf("- El pivoteo mejora la estabilidad numérica\n");
        printf("- Evita errores por división entre números muy pequeños\n");
        printf("- Reduce la propagación de errores de redondeo\n");
    } else {
        printf("❌ NO SE NECESITÓ PIVOTEO\n");
        printf("Justificación:\n");
        printf("- Todos los elementos diagonales fueron >= 1e-5\n");
        printf("- La matriz tiene buena estabilidad numérica\n");
        printf("- No hay riesgo de división por números pequeños\n");
    }
    printf("==================================================\n\n");

    // Imprimimos la matriz después de la eliminación Gaussiana (Matriz Triangular Superior)
    printf("La matriz despues de la eliminacion Gaussiana es:\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);
        }
        printf("| %10.6lf\n", b[i]);
    }
    printf("\n");


    // Verificamos el determinante de la matriz
    producto = 1.0;
    for(int i = 1; i <= n; i++) {
        producto = producto * a[i][i]; // El determinante es el producto de los elementos de la diagonal
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", producto);

    // Si el determinante es cero (o muy cercano a cero), el sistema no tiene solución única.
    if(fabs(producto) < 1e-12) {
        printf("El determinante de la matriz es cero, el sistema no tiene solución unica.\n");
        exit(0); // Salir con código de éxito (o 1 si se prefiere error)
    }

    // Realizamos la sustitución regresiva para encontrar la solución
    // Empezamos por la última incógnita
    X[n] = b[n] / a[n][n];

    // Sustitución regresiva: i = n-1 hasta 1
    for(int i = n-1; i >= 1; i--) {
        suma_temp = b[i]; // Inicializar la suma con el término independiente
        
        // Restar los términos conocidos (a[i][j] * X[j] para j > i)
        for(int j = i+1; j <= n; j++) {
            suma_temp = suma_temp - a[i][j] * X[j];
        }
        // Despejar la incógnita actual X[i]
        X[i] = suma_temp / a[i][i];
    }
    
    printf("------------------SOLUCION------------------\n");
    printf("La solucion del sistema es:\n");
    for(int i = 1; i <= n; i++) {
        printf("X[%d] = %lf\n", i, X[i]);
    }
    

    return 0;
}

// Implementación de la función de lectura
bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_TAMANO+1], double b[], int* n) {
    FILE *p_archivo;
    
    // Abrir archivo de datos
    p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }
    
    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar filas en el archivo (mejor alternativa)
    int filas = 0;
    while(!feof(p_archivo)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), p_archivo) != NULL) {
            filas++;
        }
    }
    
    // El sistema debe ser cuadrado, entonces n = filas
    *n = filas;
    printf("Tamanio del sistema: %d x %d\n", *n, *n);

    // Cerrar y reabrir el archivo para reiniciar el puntero
    fclose(p_archivo);
    p_archivo = fopen(nombre_archivo, "r");
    
    // Verificar el tamaño máximo usando MAX_TAMANO global
    if(*n > MAX_TAMANO) {
        printf("Error: Sistema demasiado grande (%d). Maximo permitido: %d\n", *n, MAX_TAMANO);
        fclose(p_archivo);
        return false;
    }

    // Leer la matriz aumentada del archivo
    // Formato esperado: cada fila contiene n coeficientes + 1 término independiente
    
    int i, j;
    for(i = 1; i <= *n; i++) {
        // Leer los coeficientes de la matriz
        for(j = 1; j <= *n; j++) {
            if(fscanf(p_archivo, "%lf", &a[i][j]) != 1) {
                printf("Error al leer el elemento a[%d][%d]\n", i, j);
                fclose(p_archivo);
                return false;
            }
        }
        // Leer el término independiente
        if(fscanf(p_archivo, "%lf", &b[i]) != 1) {
            printf("Error al leer el término independiente b[%d]\n", i);
            fclose(p_archivo);
            return false;
        }
    }
    
    fclose(p_archivo);
    printf("Matriz leída exitosamente del archivo\n\n");
    
    return true;
}
