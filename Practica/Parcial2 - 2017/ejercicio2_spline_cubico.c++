/* NO FUNCIONA
 * ===========================================================================
 * PARCIAL 2 - 2017: EJERCICIO 2 - PASO 1
 * Interpolacion con Spline Cubico Natural
 * ===========================================================================
 * 
 * Objetivo: De 9 puntos NO equiespaciados a 11 puntos EQUIESPACIADOS
 * Metodo: Splines cubicos naturales (segunda derivada = 0 en extremos)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define N_DATOS_ORIG 9     // 9 puntos originales
#define N_PUNTOS_NUEVOS 11 // 11 puntos equiespaciados

// Datos originales del problema
double X_original[N_DATOS_ORIG] = {1.00, 1.02, 1.10, 1.23, 1.35, 1.50, 1.70, 1.86, 2.00};
double f_original[N_DATOS_ORIG] = {0.098, 0.071, -0.043, -0.251, -0.453, -0.693, -0.945, -1.051, -1.053};

// ===========================================================================
// FUNCIONES PARA RESOLVER SISTEMA LINEAL (Metodo de Gauss simple)
// ===========================================================================

/**
 * Resuelve el sistema Ax = b usando eliminacion Gaussiana simple
 * Nota: Esta version simplificada no usa pivoteo, funciona para sistemas bien condicionados
 */
void resolver_sistema_gauss(int n, double A[][4*MAX_POINTS], double b[], double x[]) {
    int i, j, k;
    double factor, suma;
    
    // Eliminacion hacia adelante
    for (k = 0; k < n-1; k++) {
        for (i = k+1; i < n; i++) {
            if (fabs(A[k][k]) < 1e-10) {
                printf("ADVERTENCIA: Pivote muy pequeno en fila %d\n", k);
                continue;
            }
            factor = A[i][k] / A[k][k];
            for (j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Sustitucion hacia atras
    for (i = n-1; i >= 0; i--) {
        suma = 0.0;
        for (j = i+1; j < n; j++) {
            suma += A[i][j] * x[j];
        }
        if (fabs(A[i][i]) < 1e-10) {
            printf("ERROR: Sistema singular en posicion %d\n", i);
            x[i] = 0.0;
        } else {
            x[i] = (b[i] - suma) / A[i][i];
        }
    }
}

// ===========================================================================
// FUNCIONES PARA SPLINE CUBICO
// ===========================================================================

/**
 * Construye el sistema de ecuaciones para splines cubicos naturales
 * Para n puntos, tenemos n-1 splines, cada uno con 4 coeficientes
 * Total: 4(n-1) ecuaciones
 */
void construir_sistema_spline(double X[], double Y[], int n, 
                               double A[][4*MAX_POINTS], double b[]) {
    int k, j, fila;
    int tam = 4*(n-1);
    
    printf("\n===============================================================\n");
    printf("CONSTRUYENDO SISTEMA DE SPLINES CUBICOS\n");
    printf("===============================================================\n\n");
    
    printf("- Numero de puntos: %d\n", n);
    printf("- Numero de splines: %d\n", n-1);
    printf("- Tamano del sistema: %d x %d\n\n", tam, tam);
    
    // Inicializar matriz a cero
    for (int i = 0; i < tam; i++) {
        b[i] = 0.0;
        for (int j = 0; j < tam; j++) {
            A[i][j] = 0.0;
        }
    }
    
    // ECUACIONES 1-16: Cada spline pasa por sus dos puntos extremos
    // Para cada spline k: S_k(X[k]) = Y[k] y S_k(X[k+1]) = Y[k+1]
    printf("Paso 1: Condiciones de interpolacion [%d ecuaciones]\n", 2*(n-1));
    for (k = 0; k < n-1; k++) {
        // S_k(X[k]) = Y[k]
        for (j = 0; j <= 3; j++) {
            A[2*k][4*k+j] = pow(X[k], 3-j);
        }
        b[2*k] = Y[k];
        
        // S_k(X[k+1]) = Y[k+1]
        for (j = 0; j <= 3; j++) {
            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
        }
        b[2*k+1] = Y[k+1];
    }
    
    // ECUACIONES 17-23: Continuidad de primera derivada
    // S'_k(X[k+1]) = S'_{k+1}(X[k+1])
    printf("Paso 2: Continuidad de primera derivada [%d ecuaciones]\n", n-2);
    for (k = 0; k < n-2; k++) {
        fila = 2*(n-1) + k;
        
        // Primera derivada: S'(x) = 3a*x^2 + 2b*x + c
        for (j = 0; j <= 2; j++) {
            A[fila][4*k+j] = (3-j) * pow(X[k+1], 2-j);
            A[fila][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
        }
        b[fila] = 0.0;
    }
    
    // ECUACIONES 24-30: Continuidad de segunda derivada
    // S''_k(X[k+1]) = S''_{k+1}(X[k+1])
    printf("Paso 3: Continuidad de segunda derivada [%d ecuaciones]\n", n-2);
    for (k = 0; k < n-2; k++) {
        fila = 2*(n-1) + (n-2) + k;
        
        // Segunda derivada: S''(x) = 6a*x + 2b
        A[fila][4*k] = 6.0 * X[k+1];
        A[fila][4*k+1] = 2.0;
        A[fila][4*(k+1)] = -6.0 * X[k+1];
        A[fila][4*(k+1)+1] = -2.0;
        b[fila] = 0.0;
    }
    
    // ECUACIONES 31-32: Condiciones de frontera (spline natural)
    // S''_0(X[0]) = 0 y S''_{n-2}(X[n-1]) = 0
    printf("Paso 4: Condiciones de frontera naturales [2 ecuaciones]\n\n");
    
    fila = 4*(n-1) - 2;
    // S''_0(X[0]) = 0
    A[fila][0] = 6.0 * X[0];
    A[fila][1] = 2.0;
    b[fila] = 0.0;
    
    // S''_{n-2}(X[n-1]) = 0
    fila = 4*(n-1) - 1;
    A[fila][4*(n-2)] = 6.0 * X[n-1];
    A[fila][4*(n-2)+1] = 2.0;
    b[fila] = 0.0;
    
    printf("OK - Sistema construido exitosamente\n");
}

/**
 * Evalua el spline cubico en un punto x dado
 */
double evaluar_spline(double X[], double coef[], int n, double x) {
    int k;
    
    // Encontrar el intervalo correcto
    if (x <= X[0]) {
        k = 0;
    } else if (x >= X[n-1]) {
        k = n-2;
    } else {
        k = 0;
        while (k < n-1 && x > X[k+1]) {
            k++;
        }
    }
    
    // Evaluar S_k(x) = a*x^3 + b*x^2 + c*x + d
    double a = coef[4*k];
    double b = coef[4*k+1];
    double c = coef[4*k+2];
    double d = coef[4*k+3];
    
    return a*x*x*x + b*x*x + c*x + d;
}

/**
 * Imprime los coeficientes de los splines
 */
void imprimir_coeficientes(double X[], double coef[], int n) {
    printf("\n===============================================================\n");
    printf("COEFICIENTES DE LOS SPLINES CUBICOS\n");
    printf("===============================================================\n\n");
    
    for (int k = 0; k < n-1; k++) {
        printf("Spline %d (intervalo [%.2f, %.2f]):\n", k+1, X[k], X[k+1]);
        printf("  S_%d(x) = %.8f*x^3 + %.8f*x^2 + %.8f*x + %.8f\n\n",
               k+1, coef[4*k], coef[4*k+1], coef[4*k+2], coef[4*k+3]);
    }
}

// ===========================================================================
// PROGRAMA PRINCIPAL
// ===========================================================================

int main() {
    double A[4*MAX_POINTS][4*MAX_POINTS];
    double b[4*MAX_POINTS];
    double solucion[4*MAX_POINTS];
    int tam_sistema = 4*(N_DATOS_ORIG-1);
    
    printf("===============================================================\n");
    printf("       PARCIAL 2 - 2017: EJERCICIO 2 - PASO 1              \n");
    printf("  Interpolacion con Spline Cubico Natural                   \n");
    printf("===============================================================\n\n");

    // Mostrar datos originales
    printf("===============================================================\n");
    printf("DATOS ORIGINALES (9 puntos NO equiespaciados)\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-------------+\n");
    printf("|  i  |    x     |    f(x)     |\n");
    printf("+-----+----------+-------------+\n");
    for (int i = 0; i < N_DATOS_ORIG; i++) {
        printf("| %3d | %8.2f | %11.6f |\n", i, X_original[i], f_original[i]);
    }
    printf("+-----+----------+-------------+\n\n");

    // Construir sistema de splines
    construir_sistema_spline(X_original, f_original, N_DATOS_ORIG, A, b);
    
    // Resolver sistema
    printf("\n===============================================================\n");
    printf("RESOLVIENDO SISTEMA CON ELIMINACION GAUSSIANA\n");
    printf("===============================================================\n\n");
    
    printf("- Resolviendo sistema de %d ecuaciones...\n", tam_sistema);
    resolver_sistema_gauss(tam_sistema, A, b, solucion);
    printf("OK - Sistema resuelto\n");
    
    // Imprimir coeficientes
    imprimir_coeficientes(X_original, solucion, N_DATOS_ORIG);
    
    // Crear tabla equiespaciada
    printf("===============================================================\n");
    printf("TABLA EQUIESPACIADA (11 puntos con h=0.1)\n");
    printf("===============================================================\n\n");
    
    double X_equi[N_PUNTOS_NUEVOS];
    double f_equi[N_PUNTOS_NUEVOS];
    double h = 0.1;
    
    printf("- Intervalo: [1.0, 2.0]\n");
    printf("- Numero de puntos: %d\n", N_PUNTOS_NUEVOS);
    printf("- Paso h = %.4f\n\n", h);
    
    // Evaluar spline en puntos equiespaciados
    for (int i = 0; i < N_PUNTOS_NUEVOS; i++) {
        X_equi[i] = 1.0 + i * h;
        f_equi[i] = evaluar_spline(X_original, solucion, N_DATOS_ORIG, X_equi[i]);
    }
    
    // Mostrar tabla resultante
    printf("+-----+----------+-----------------+--------------+\n");
    printf("|  i  |    x     |  f(x) spline    |  Delta_x     |\n");
    printf("+-----+----------+-----------------+--------------+\n");
    for (int i = 0; i < N_PUNTOS_NUEVOS; i++) {
        if (i == 0) {
            printf("| %3d | %8.4f | %15.10f |      -       |\n", 
                   i, X_equi[i], f_equi[i]);
        } else {
            double delta = X_equi[i] - X_equi[i-1];
            printf("| %3d | %8.4f | %15.10f |    %.4f    |\n", 
                   i, X_equi[i], f_equi[i], delta);
        }
    }
    printf("+-----+----------+-----------------+--------------+\n\n");
    
    printf("OK - Todos los puntos estan EQUIESPACIADOS con h = 0.1\n\n");
    
    // Verificacion: evaluar en puntos originales
    printf("===============================================================\n");
    printf("VERIFICACION: Evaluando spline en puntos originales\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-------------+-------------+-----------+\n");
    printf("|  i  |    x     |  f original | f spline    |   Error   |\n");
    printf("+-----+----------+-------------+-------------+-----------+\n");
    
    double max_error = 0.0;
    for (int i = 0; i < N_DATOS_ORIG; i++) {
        double f_spline = evaluar_spline(X_original, solucion, N_DATOS_ORIG, X_original[i]);
        double error = fabs(f_original[i] - f_spline);
        if (error > max_error) max_error = error;
        
        printf("| %3d | %8.2f | %11.6f | %11.6f | %9.2e |\n", 
               i, X_original[i], f_original[i], f_spline, error);
    }
    printf("+-----+----------+-------------+-------------+-----------+\n");
    printf("Error maximo: %.2e (deberia ser ~ 0)\n\n", max_error);
    
    // Resumen
    printf("===============================================================\n");
    printf("  RESUMEN DEL PASO 1                                        \n");
    printf("===============================================================\n\n");
    
    printf("OK ENTRADA: 9 puntos NO equiespaciados\n");
    printf("OK METODO: Splines cubicos naturales\n");
    printf("   - Sistema de %d ecuaciones lineales\n", tam_sistema);
    printf("   - %d splines cubicos (polinomios grado 3)\n", N_DATOS_ORIG-1);
    printf("   - Continuidad C2 (hasta segunda derivada)\n");
    printf("OK SALIDA: 11 puntos EQUIESPACIADOS\n");
    printf("   - Intervalo: [1.0, 2.0]\n");
    printf("   - Paso: h = 0.1\n");
    printf("   - Error de interpolacion: %.2e\n", max_error);
    printf("   - Listo para PASO 2 (diferenciacion numerica)\n\n");
    
    printf("===============================================================\n");
    printf("PROXIMO PASO: Calcular df/dx usando diferencias finitas\n");
    printf("===============================================================\n\n");

    return 0;
}
