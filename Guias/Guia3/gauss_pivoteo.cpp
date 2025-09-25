#include <cstdio>  
#include <stdlib.h> 
#include <cmath>   

// Ahora puedes manejar matrices de hasta 50x50
#define MAX_TAMANO 50  // Tamaño máximo de la matriz y vectores

/**
 * Lee una matriz aumentada desde un archivo de texto.
 * El formato esperado es:
 * a11 a12 ... a1n b1
 * a21 a22 ... a2n b2
 * ...
 * an1 an2 ... ann bn
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
    // Bandera para saber si se usó pivoteo parcial
    int se_uso_pivoteo_parcial = 0;
    
    // Definición de arreglos usando el MAX_TAMANO global
    double a[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], X[MAX_TAMANO+1]; // X es el vector solución

    // Leer arreglo desde el archivo usando la función. Se espera 'data2.dat'
    if(!leer_archivo_matriz("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia3\\data.dat", a, b, &n)) {
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

        // ¿Qué sucede si a[i][i] está cerca de cero?
        // Usaremos pivoteo parcial para evitar la división por cero o la inestabilidad numérica
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            se_uso_pivoteo_parcial = 1;
            for(int l = i+1; l <= n; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l; // Encontramos la fila con el elemento más grande en la columna i
                }
            }
            
            // Si la fila con el pivote más grande es diferente, intercambiamos filas
            if (p != i) { 
                for(int m = i; m <= n; m++) {
                    aux = a[p][m];
                    a[p][m] = a[i][m];
                    a[i][m] = aux; // Intercambiar las filas p e i en la matriz A
                }
                aux = b[p];
                b[p] = b[i];
                b[i] = aux; // Intercambiar el término independiente
            } else if (fabs(a[i][i]) < 1e-5) {
                // Si incluso después del pivoteo el pivote es cero, el sistema no tiene solución única
                printf("Error: Pivote cero o casi cero encontrado en la matriz. El sistema no tiene solución única.\n");
                exit(1);
            }
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

    // Imprimir si se usó pivoteo parcial en la resolución de la matriz
    if(se_uso_pivoteo_parcial == 1) {
        printf("Se utilizó pivoteo parcial.\n");
    }

    // Imprimimos la matriz después de la eliminación Gaussiana (Matriz Triangular Superior)
    printf("La matriz después de la eliminación Gaussiana es:\n");
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
        printf("El determinante de la matriz es cero, el sistema no tiene solución única.\n");
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
    
    printf("------------------SOLUCIÓN------------------\n");
    printf("La solución del sistema es:\n");
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
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
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
    printf("📊 Tamaño del sistema: %d x %d\n", *n, *n);

    // Cerrar y reabrir el archivo para reiniciar el puntero
    fclose(p_archivo);
    p_archivo = fopen(nombre_archivo, "r");
    
    // Verificar el tamaño máximo usando MAX_TAMANO global
    if(*n > MAX_TAMANO) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_TAMANO);
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
                printf("❌ Error al leer el elemento a[%d][%d]\n", i, j);
                fclose(p_archivo);
                return false;
            }
        }
        // Leer el término independiente
        if(fscanf(p_archivo, "%lf", &b[i]) != 1) {
            printf("❌ Error al leer el término independiente b[%d]\n", i);
            fclose(p_archivo);
            return false;
        }
    }
    
    fclose(p_archivo);
    printf("✅ Matriz leída exitosamente del archivo\n\n");
    
    return true;
}