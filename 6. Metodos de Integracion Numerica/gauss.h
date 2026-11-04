#ifndef ELIMINACION_GAUSSIANA_H
#define ELIMINACION_GAUSSIANA_H

// Incluye las librerías estándar necesarias
#include <stdio.h> // Para entrada/salida (printf)
#include <stdlib.h> // Para funciones generales (exit)
#include <math.h> // Para funciones matemáticas (fabs)

// Define el tamaño máximo de la matriz (100x100 + 1 para manejar índices desde 1 si fuera necesario, aunque el código usa 0-indexing)
#ifndef MAX_SIZE
#define MAX_SIZE 100
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
    double sign = 1.0; // Guarda el signo del determinante cuando se intercambian filas
    const double EPS_PIVOTE = 1e-12; // Tolerancia para considerar un pivote nulo
    const double EPS_DETERMINANTE = 1e-10; // Tolerancia para considerar el determinante nulo

    // Validar el tamaño del sistema para evitar índices fuera de rango
    if (n <= 0) {
        printf("Error: el tamaño del sistema debe ser mayor que 0.\n");
        exit(1);
    }

    // --- FASE 1: ELIMINACIÓN GAUSSIANA (Conversión a Matriz Triangular Superior) ---
    // Recorre las filas de la matriz para la eliminación
    for(int i = 0; i <= n-2; i++) {

        // --- Pivoteo Parcial ---
        // Busca SIEMPRE la fila 'p' con el mayor valor absoluto en la columna 'i'
        // (desde la fila 'i' hacia abajo), no solo cuando el pivote es chico.
        p = i;
        for(int l = i+1; l <= n-1; l++) {
            if(fabs(a[l][i]) > fabs(a[p][i])) {
                p = l;
            }
        }

        // Intercambia la fila actual 'i' con la fila 'p' (la fila del nuevo pivote)
        if(p != i) {
            for(int m = i; m <= n-1; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }

            // Intercambia también los elementos correspondientes en el vector de resultados 'b'
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;

            // Cada intercambio de filas cambia el signo del determinante
            sign = -sign;
        }

        // Si el pivote es nulo (o casi nulo), el sistema no tiene solución única
        if(fabs(a[i][i]) < EPS_PIVOTE) {
            printf("Error: pivote nulo en la columna %d, no hay solución única.\n", i);
            exit(1);
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

    // El determinante de una matriz triangular es el producto de los elementos de la diagonal.
    // Se multiplica por 'sign' para corregir el signo si hubo intercambios de filas.
    product = sign;
    for(int i = 0; i < n; i++) {
        product *= a[i][i];
    }

    printf("El determinante de la matriz es: %lf\n\n", product);

    // Si el determinante es cercano a cero, el sistema no tiene solución única.
    if(fabs(product) < EPS_DETERMINANTE) {
        printf("Error: el determinante es cero, no hay solucion unica.\n");
        exit(1); // Sale del programa con código de error
    }

    // --- FASE 2: SUSTITUCIÓN HACIA ATRÁS (Resolución del sistema triangular) ---
    // Calcula el último valor de la incógnita (X[n-1])
    X[n-1] = b[n-1] / a[n-1][n-1];

    // Recorre las filas de abajo hacia arriba, desde la penúltima hasta la primera
    for(int i = n-2; i >= 0; i--) {
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
