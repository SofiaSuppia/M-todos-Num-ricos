#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_TAMANO 100 // Ajusta si lo necesitas

/**
 * Función para resolver un sistema de ecuaciones lineales mediante Eliminación Gaussiana
 * con pivoteo parcial y sustitución hacia atrás.
 * * @param n El tamaño del sistema (número de ecuaciones/variables)
 * @param matriz_a La matriz de coeficientes del sistema
 * @param vector_b El vector de términos independientes
 * @param vector_solucion El arreglo para almacenar la solución (X)
 */
void eliminacion_gaussiana(int n, double matriz_a[MAX_TAMANO+1][MAX_TAMANO+1], double vector_b[MAX_TAMANO+1], double vector_solucion[MAX_TAMANO+1]) {
    int p_pivote;
    double factor, producto, suma, aux;

    // Recorrer las filas de la matriz (Eliminación Gaussiana)
    for(int i = 0; i <= n-2; i++) {
        // Pivoteo parcial
        p_pivote = i;
        if(fabs(matriz_a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(matriz_a[l][i]) > fabs(matriz_a[p_pivote][i])) {
                    p_pivote = l;
                }
            }
            // Intercambio de filas en la matriz
            for(int m = i; m <= n-1; m++) {
                aux = matriz_a[p_pivote][m];
                matriz_a[p_pivote][m] = matriz_a[i][m];
                matriz_a[i][m] = aux;
            }
            // Intercambio de elementos en el vector b
            aux = vector_b[p_pivote];
            vector_b[p_pivote] = vector_b[i];
            vector_b[i] = aux;
        }

        // Hacer cero los elementos debajo de la diagonal
        for(int j = i+1; j <= n-1; j++) {
            factor = matriz_a[j][i] / matriz_a[i][i];
            for(int k = i; k <= n-1; k++) {
                matriz_a[j][k] -= factor * matriz_a[i][k];
            }
            vector_b[j] -= factor * vector_b[i];
        }
    }

    // Imprimimos la matriz después de la eliminación Gaussiana, es más conveniente
    printf("La matriz después de la eliminación Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", matriz_a[i][j]);
        }
        printf("\n");
    }
    printf("El vector b después de la eliminación Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", vector_b[i]);
    }
    printf("\n\n");

    // Verificación del determinante (el producto de los elementos de la diagonal principal)
    producto = 1.0;
    for(int i = 0; i < n; i++) {
        producto *= matriz_a[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", producto);

    if(fabs(producto) < 1e-10) {
        printf("Error: el determinante es cero, no hay solución única.\n");
        exit(1);
    }

    // Sustitución hacia atrás
    vector_solucion[n-1] = vector_b[n-1] / matriz_a[n-1][n-1];
    for(int i = n-2; i >= 0; i--) {
        suma = vector_b[i];
        for(int j = i+1; j < n; j++) {
            suma -= matriz_a[i][j] * vector_solucion[j];
        }
        vector_solucion[i] = suma / matriz_a[i][i];
    }
}

#endif