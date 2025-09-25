#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_TAMANO 100 // Tamaño máximo de la matriz

/**
 * @brief Resuelve un sistema de ecuaciones lineales A*X = B mediante eliminación Gaussiana con pivoteo parcial.
 * * @param n El tamaño del sistema (número de ecuaciones/variables).
 * @param a La matriz de coeficientes (se modifica durante la eliminación).
 * @param b El vector de términos independientes (se modifica durante la eliminación).
 * @param X El vector para almacenar la solución final.
 */
void eliminacion_gaussiana(int n, double a[MAX_TAMANO+1][MAX_TAMANO+1], double b[MAX_TAMANO+1], double X[MAX_TAMANO+1]) {
    int p_pivote; // Índice para el pivote
    double factor, producto, suma, aux;

    // Recorrer las filas de la matriz (Eliminación Gaussiana)
    for(int i = 0; i <= n-2; i++) {
        // Pivoteo parcial
        p_pivote = i;
        if(fabs(a[i][i]) < 1e-5) {
            // Buscar el elemento con el mayor valor absoluto en la columna i
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(a[l][i]) > fabs(a[p_pivote][i])) {
                    p_pivote = l;
                }
            }
            // Intercambiar la fila actual (i) con la fila del pivote (p_pivote) en la matriz a
            for(int m = i; m <= n-1; m++) {
                aux = a[p_pivote][m];
                a[p_pivote][m] = a[i][m];
                a[i][m] = aux;
            }
            // Intercambiar los elementos correspondientes en el vector b
            aux = b[p_pivote];
            b[p_pivote] = b[i];
            b[i] = aux;
        }

        // Hacer cero los elementos debajo de la diagonal (reducción por filas)
        for(int j = i+1; j <= n-1; j++) {
            factor = a[j][i] / a[i][i];
            for(int k = i; k <= n-1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Imprimimos la matriz después de la eliminación Gaussiana (matriz triangular superior), es más conveniente
    printf("La matriz después de la eliminación Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("El vector b después de la eliminación Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // Verificación del determinante (el producto de los elementos de la diagonal principal de la matriz triangular)
    producto = 1.0;
    for(int i = 0; i < n; i++) {
        producto *= a[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", producto);

    if(fabs(producto) < 1e-10) {
        printf("Error: el determinante es cero, no hay solución única.\n");
        exit(1);
    }

    // Sustitución hacia atrás
    // Calcular X[n-1]
    X[n-1] = b[n-1] / a[n-1][n-1];
    // Calcular el resto de las variables (X[i] para i = n-2 hasta 0)
    for(int i = n-2; i >= 0; i--) {
        suma = b[i];
        for(int j = i+1; j < n; j++) {
            suma -= a[i][j] * X[j];
        }
        X[i] = suma / a[i][i];
    }
}

#endif