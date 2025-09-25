#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100

void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int p;
    double factor, product, sum, aux;

    // Walk the rows of the matrix (Gaussian elimination)
    for(int i = 0; i <= n-2; i++) {
        // Partial pivoting
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l;
                }
            }
            for(int m = i; m <= n-1; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        // Make zero the elements below the diagonal
        for(int j = i+1; j <= n-1; j++) {
            factor = a[j][i] / a[i][i];
            for(int k = i; k <= n-1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // We print the matrix after Gaussian elimination, It's more convenient
    printf("The matrix after Gaussian elimination is:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("The vector b after Gaussian elimination is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // Determinant check
    product = 1.0;
    for(int i = 0; i < n; i++) {
        product *= a[i][i];
    }

    printf("------------------DETERMINANT------------------\n");
    printf("The determinant of the matrix is: %lf\n\n", product);

    if(fabs(product) < 1e-10) {
        printf("Error: determinant is zero, no unique solution.\n");
        exit(1);
    }

    // Back substitution
    X[n-1] = b[n-1] / a[n-1][n-1];
    for(int i = n-2; i >= 0; i--) {
        sum = b[i];
        for(int j = i+1; j < n; j++) {
            sum -= a[i][j] * X[j];
        }
        X[i] = sum / a[i][i];
    }
}

#endif
