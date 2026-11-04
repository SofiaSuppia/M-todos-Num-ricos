#ifndef ELIMINACION_GAUSSIANA_H
#define ELIMINACION_GAUSSIANA_H

/*
 * ============================================================================
 * ELIMINACION GAUSSIANA CON SUSTITUCION HACIA ATRAS
 * ============================================================================
 *
 * MAPA DE CALCULOS
 * ----------------
 * - DERIVADAS: este archivo no calcula derivadas.
 * - INTEGRALES: este archivo no calcula integrales.
 * - SUMATORIAS:
 *     La unica sumatoria esta en la sustitucion hacia atras. Para cada fila i
 *     se acumula:
 *         sum = b[i] - a[i][i+1]*X[i+1] - ... - a[i][n-1]*X[n-1]
 *     usando el bucle "for (j = i+1; j < n; j++) sum -= a[i][j] * X[j];"
 * - PRODUCTO:
 *     El determinante se calcula multiplicando la diagonal de la matriz
 *     triangular: product *= a[i][i].
 * ============================================================================
 */

// Incluye las librerías estándar necesarias
#include <stdio.h> // Para entrada/salida (printf)
#include <stdlib.h> // Para funciones generales (exit)
#include <math.h> // Para funciones matemáticas (fabs)

/*
 * Tamano maximo de la matriz.
 * Si el archivo que incluye gauss.h ya definio MAX_SIZE, se respeta ese valor.
 * Por defecto se usa 200 para que el spline cubico pueda manejar hasta 50
 * puntos: 4*(50-1) = 196 incognitas.
 */
#ifndef MAX_SIZE
#define MAX_SIZE 200
#endif

/**
 * @brief Resuelve un sistema de ecuaciones lineales Ax = b utilizando la eliminación de Gauss con pivoteo parcial.
 * * @param n El tamaño de la matriz (n x n).
 * @param a La matriz de coeficientes aumentada (se modifica in-place a forma triangular superior).
 * @param b El vector de resultados (se modifica in-place).
 * @param X El vector solución donde se almacenarán los resultados.
 */
void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int p; // Índice para almacenar la fila del pivote más grande
    double factor, product, sum, aux; // Variables auxiliares para los cálculos

    // --- FASE 1: ELIMINACIÓN GAUSSIANA (Conversión a Matriz Triangular Superior) ---
    // Recorre las filas de la matriz para la eliminación
    for(int i = 0; i <= n-2; i++) {
        
        // --- Pivoteo Parcial ---
        p = i;
        // Comprueba si el pivote actual es demasiado pequeño (cercano a cero)
        if(fabs(a[i][i]) < 1e-5) {
            // Busca la fila 'p' con el elemento más grande en la columna 'i' debajo del pivote actual (i, i)
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l; // Actualiza el índice del pivote más grande
                }
            }
            
            // Intercambia la fila actual 'i' con la fila 'p' (la fila del nuevo pivote)
            for(int m = i; m <= n-1; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }
            
            // Intercambia también los elementos correspondientes en el vector de resultados 'b'
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        // --- Eliminación hacia adelante ---
        // Pone a cero los elementos debajo de la diagonal (columna 'i')
        for(int j = i+1; j <= n-1; j++) {
            // Calcula el factor de eliminación
            factor = a[j][i] / a[i][i];
            
            // Actualiza los elementos de la fila 'j' (a[j] = a[j] - factor * a[i])
            for(int k = i; k <= n-1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            
            // Actualiza el elemento del vector 'b' correspondiente
            b[j] -= factor * b[i];
        }
    }

    // Imprimimos la matriz después de la eliminación Gaussiana (ya en forma triangular superior)
    printf("La matriz despues de la eliminacion Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("El vector b despues de la eliminacion Gaussiana es:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // --- Verificación del Determinante ---
    printf("------------------DETERMINANTE------------------\n");
    product = 1.0; // Inicializa el producto para calcular el determinante
    
    // El determinante de una matriz triangular es el producto de los elementos de la diagonal.
    for(int i = 0; i < n; i++) {
        product *= a[i][i];
    }

    printf("El determinante de la matriz es: %lf\n\n", product);

    // Si el determinante es cercano a cero, el sistema no tiene solución única.
    if(fabs(product) < 1e-10) {
        printf("Error: el determinante es cero, no hay solucion unica.\n");
        exit(1); // Sale del programa con código de error
    }

    // --- FASE 2: SUSTITUCIÓN HACIA ATRÁS (Resolución del sistema triangular) ---
    // Calcula el último valor de la incógnita (X[n-1])
    X[n-1] = b[n-1] / a[n-1][n-1];
    
    // Recorre las filas de abajo hacia arriba, desde la penúltima hasta la primera
    for(int i = n-2; i >= 0; i--) {
        // SUMATORIA: comienza con b[i] y luego resta los terminos conocidos.
        sum = b[i]; // Inicializa la suma con el valor del vector b

        // Resta los términos de las incógnitas ya calculadas (de derecha a izquierda)
        for(int j = i+1; j < n; j++) {
            sum -= a[i][j] * X[j];
        }
        
        // Calcula la incógnita actual X[i]
        X[i] = sum / a[i][i];
    }
}

#endif