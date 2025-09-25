#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100  // Ajústalo si lo necesitas

void eliminacion_gaussiana(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int p;
    double factor, producto, suma, auxiliar;

    // Recorrer las filas de la matriz (Eliminación Gaussiana)
    for(int i = 0; i <= n-2; i++) {
        // Pivoteo parcial
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l;
                }
            }
            // Intercambiar filas de la matriz A
            for(int m = i; m <= n-1; m++) {
                auxiliar = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = auxiliar;
            }
            // Intercambiar elementos del vector b
            auxiliar = b[p];
            b[p] = b[i];
            b[i] = auxiliar;
        }

        // Hacer cero los elementos debajo de la diagonal
        for(int j = i+1; j <= n-1; j++) {
            factor = a[j][i] / a[i][i];
            for(int k = i; k <= n-1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Imprimimos la matriz después de la eliminación Gaussiana, es más conveniente
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

    // Verificación del determinante
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
    X[n-1] = b[n-1] / a[n-1][n-1];
    for(int i = n-2; i >= 0; i--) {
        suma = b[i];
        for(int j = i+1; j < n; j++) {
            suma -= a[i][j] * X[j];
        }
        X[i] = suma / a[i][i];
    }
}

#endif