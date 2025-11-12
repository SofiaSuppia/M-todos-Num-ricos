#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definimos el tamaño máximo de la matriz para facilitar la declaración de arrays.
#define MAX_SIZE 100

/**
 * @brief Resuelve un sistema de ecuaciones lineales A * X = B utilizando
 * la Eliminación Gaussiana con Pivoteo Parcial.
 *
 * @param n Tamaño de la matriz cuadrada A (número de filas/columnas).
 * @param a Matriz de coeficientes [n x n].
 * @param b Vector de términos independientes [n].
 * @param X Vector donde se almacenarán las soluciones [n].
 */
void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int p; // Índice para el pivote
    double factor, producto, suma, aux;

    // --- FASE 1: ELIMINACIÓN HACIA ADELANTE (Triangularización) ---
    // Recorremos las filas de la matriz hasta la penúltima (i = 0 a n-2)
    for(int i = 0; i <= n-2; i++) {
        
        // 1. Pivoteo Parcial: Buscamos el elemento con mayor magnitud en la columna i,
        // desde la fila actual i hacia abajo.
        p = i;
        // Si el pivote actual es muy cercano a cero (problema potencial de división por cero o error)
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i + 1; l <= n - 1; l++) {
                // Si encontramos un valor absoluto más grande, actualizamos el índice del pivote.
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l;
                }
            }

            // Intercambiamos la fila actual (i) con la fila pivote (p)
            for(int m = i; m <= n - 1; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }
            // Intercambiamos también los elementos correspondientes en el vector b
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        // 2. Eliminación: Convertimos en cero los elementos debajo de la diagonal (a[j][i])
        for(int j = i + 1; j <= n - 1; j++) {
            // Calculamos el factor de multiplicación para la eliminación (m_ji)
            factor = a[j][i] / a[i][i];

            // Aplicamos la operación de fila: Fila_j = Fila_j - factor * Fila_i
            for(int k = i; k <= n - 1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            // Aplicamos la misma operación al vector b
            b[j] -= factor * b[i];
        }
    }

    // --- RESULTADOS DE LA FASE 1 ---
    
    // Imprimimos la matriz después de la Eliminación Gaussiana (ahora triangular superior)
    printf("La matriz después de la Eliminación Gaussiana (Triangular Superior) es:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("El vector b después de la Eliminación Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // Cálculo del Determinante: Es el producto de los elementos de la diagonal principal
    producto = 1.0;
    for(int i = 0; i < n; i++) {
        producto *= a[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", producto);

    // Verificación de existencia de solución única
    if(fabs(producto) < 1e-10) {
        printf("Error: el determinante es cero. No hay solución única.\n");
        // Terminamos el programa si no se puede resolver el sistema de forma única
        exit(1); 
    }

    // --- FASE 2: SUSTITUCIÓN HACIA ATRÁS ---
    
    // 1. Calculamos la última incógnita (X[n-1])
    // X[n-1] = b[n-1] / a[n-1][n-1]
    X[n-1] = b[n-1] / a[n-1][n-1];
    
    // 2. Calculamos las incógnitas restantes de forma ascendente (i = n-2 hasta 0)
    for(int i = n - 2; i >= 0; i--) {
        // Inicializamos la suma con el elemento b[i]
        suma = b[i];
        
        // Restamos los términos conocidos que ya hemos calculado (a[i][j] * X[j])
        for(int j = i + 1; j < n; j++) {
            suma -= a[i][j] * X[j];
        }
        
        // Despejamos la incógnita actual X[i]
        X[i] = suma / a[i][i];
    }
}

#endif // GAUSS_H