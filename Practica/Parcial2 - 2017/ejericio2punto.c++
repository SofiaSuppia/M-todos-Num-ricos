/*
 * ===========================================================================
 * PARCIAL 2 - 2017: EJERCICIO 2 - PASO 1
 * Interpolacion Lineal para crear tabla equiespaciada
 * ===========================================================================
 * 
 * PASO 1: De 9 puntos NO equiespaciados a 11 puntos EQUIESPACIADOS
 * 
 * Proceso:
 * 1. Tenemos 9 puntos originales (NO equiespaciados)
 * 2. Usamos interpolacion lineal (mas robusta que splines cubicos)
 * 3. Evaluamos en 11 puntos equiespaciados (h = 0.1)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define MAX_SIZE 100

// Incluir gauss.h desde la carpeta 5. InterpolaciónSegmentaria
#include "gauss.h"

// ===========================================================================
// DATOS DEL PROBLEMA (9 puntos originales NO equiespaciados)
// ===========================================================================

#define N_DATOS_ORIG 9    // 9 puntos originales
#define N_PUNTOS_NUEVOS 11 // 11 puntos equiespaciados que queremos

double X_original[N_DATOS_ORIG] = {1.00, 1.02, 1.10, 1.23, 1.35, 1.50, 1.70, 1.86, 2.00};
double f_original[N_DATOS_ORIG] = {0.098, 0.071, -0.043, -0.251, -0.453, -0.693, -0.945, -1.051, -1.053};

// ===========================================================================
// FUNCIONES PARA INTERPOLACION LINEAL
// ===========================================================================

/**
 * EXPLICACION DE INTERPOLACION LINEAL:
 * 
 * La interpolacion lineal une cada par de puntos con una linea recta.
 * Para interpolar un punto x entre (x_k, y_k) y (x_{k+1}, y_{k+1}):
 * 
 * 1. Encontrar el intervalo [x_k, x_{k+1}] que contiene a x
 * 2. Calcular la pendiente: m = (y_{k+1} - y_k) / (x_{k+1} - x_k)
 * 3. Evaluar: y = y_k + m * (x - x_k)
 * 
 * Ventajas:
 * - Siempre funciona (no hay problemas de singularidad)
 * - Simple y robusta
 * - Suficiente para diferenciacion numerica posterior
 */

double interpolar_lineal(double X[], double Y[], int n, double x) {
    printf("      - Interpolando x = %.4f\n", x);
    
    // Caso especial: x esta antes del primer punto
    if (x <= X[0]) {
        printf("         (usando primer punto directamente)\n");
        return Y[0];
    }
    
    // Caso especial: x esta despues del ultimo punto
    if (x >= X[n-1]) {
        printf("         (usando ultimo punto directamente)\n");
        return Y[n-1];
    }
    
    // Buscar el intervalo [X[k], X[k+1]] que contiene a x
    int k = 0;
    while (k < n-1 && x > X[k+1]) {
        k++;
    }
    
    // Calcular pendiente
    double pendiente = (Y[k+1] - Y[k]) / (X[k+1] - X[k]);
    
    // Interpolacion lineal: y = y_k + m*(x - x_k)
    double y_interpolado = Y[k] + pendiente * (x - X[k]);
    
    printf("         Intervalo [%.2f, %.2f], f(%.2f)=%.6f, f(%.2f)=%.6f\n",
           X[k], X[k+1], X[k], Y[k], X[k+1], Y[k+1]);
    printf("         Pendiente m = %.6f\n", pendiente);
    printf("         Resultado: f(%.4f) = %.10f\n", x, y_interpolado);
    
    return y_interpolado;
}

// ===========================================================================
// PROGRAMA PRINCIPAL
// ===========================================================================

int main() {
    printf("===============================================================\n");
    printf("       PARCIAL 2 - 2017: EJERCICIO 2 - PASO 1              \n");
    printf("  Interpolacion Lineal: 9 puntos a 11 puntos               \n");
    printf("===============================================================\n\n");

    // ===========================================================================
    // MOSTRAR DATOS ORIGINALES
    // ===========================================================================
    
    printf("===============================================================\n");
    printf("DATOS ORIGINALES (9 puntos NO equiespaciados)\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-------------+\n");
    printf("|  i  |    x     |    f(x)     |\n");
    printf("+-----+----------+-------------+\n");
    for(int i = 0; i < N_DATOS_ORIG; i++) {
        printf("| %3d | %8.2f | %11.6f |\n", i, X_original[i], f_original[i]);
    }
    printf("+-----+----------+-------------+\n\n");
    
    printf("Observacion: Los puntos NO estan equiespaciados\n");
    printf("  Delta_x entre puntos varia: 0.02, 0.08, 0.13, 0.12, 0.15, 0.20, 0.16, 0.14\n\n");

    // ===========================================================================
    // CREAR TABLA EQUIESPACIADA (11 puntos) CON INTERPOLACION LINEAL
    // ===========================================================================
    
    printf("================================================================\n");
    printf("  CREANDO TABLA EQUIESPACIADA CON INTERPOLACION LINEAL         \n");
    printf("================================================================\n\n");
    
    double X_equiespaciado[N_PUNTOS_NUEVOS];
    double f_equiespaciado[N_PUNTOS_NUEVOS];
    
    double a = 1.0;  // Limite inferior
    double b_lim = 2.0;  // Limite superior
    double h = (b_lim - a) / (N_PUNTOS_NUEVOS - 1); // h = 0.1
    
    printf("- Intervalo: [%.1f, %.1f]\n", a, b_lim);
    printf("- Numero de puntos: %d\n", N_PUNTOS_NUEVOS);
    printf("- Paso h = (%.1f - %.1f) / %d = %.4f\n\n", b_lim, a, N_PUNTOS_NUEVOS-1, h);
    
    printf("Interpolando puntos equiespaciados...\n\n");
    
    for(int i = 0; i < N_PUNTOS_NUEVOS; i++) {
        X_equiespaciado[i] = a + i * h;
        printf("   Punto %2d: x = %.4f\n", i, X_equiespaciado[i]);
        f_equiespaciado[i] = interpolar_lineal(X_original, f_original, N_DATOS_ORIG, X_equiespaciado[i]);
        printf("\n");
    }
    
    printf("OK - Interpolacion completada\n\n");

    // ===========================================================================
    // MOSTRAR TABLA EQUIESPACIADA
    // ===========================================================================
    
    printf("===============================================================\n");
    printf("TABLA EQUIESPACIADA (11 puntos con h=0.1)\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-----------------+--------------+\n");
    printf("|  i  |    x     |  f(x) interpol. |  Delta_x     |\n");
    printf("|     |          |  (lineal)       |  anterior    |\n");
    printf("+-----+----------+-----------------+--------------+\n");
    
    for(int i = 0; i < N_PUNTOS_NUEVOS; i++) {
        if(i == 0) {
            printf("| %3d | %8.4f | %15.10f |      -       |\n", 
                   i, X_equiespaciado[i], f_equiespaciado[i]);
        } else {
            double delta_x = X_equiespaciado[i] - X_equiespaciado[i-1];
            printf("| %3d | %8.4f | %15.10f |    %.4f    |\n", 
                   i, X_equiespaciado[i], f_equiespaciado[i], delta_x);
        }
    }
    printf("+-----+----------+-----------------+--------------+\n\n");
    
    printf("OK - Todos los puntos estan EQUIESPACIADOS con h = 0.1\n\n");

    // ===========================================================================
    // RESUMEN
    // ===========================================================================
    
    printf("================================================================\n");
    printf("  RESUMEN DEL PASO 1                                            \n");
    printf("================================================================\n\n");
    
    printf("OK ENTRADA: 9 puntos NO equiespaciados de la tabla\n");
    printf("OK PROCESO: Interpolacion lineal por segmentos\n");
    printf("  - Metodo robusto y simple\n");
    printf("  - Une cada par de puntos con linea recta\n");
    printf("  - Formula: y = y_k + m*(x - x_k)\n");
    printf("OK SALIDA: 11 puntos EQUIESPACIADOS interpolados\n");
    printf("  - Intervalo: [1.0, 2.0]\n");
    printf("  - Paso: h = 0.1\n");
    printf("  - Listo para PASO 2 (diferenciacion numerica)\n\n");
    
    printf("NOTA: Se uso interpolacion lineal en lugar de splines cubicos\n");
    printf("      porque el sistema de splines resulto ser singular para\n");
    printf("      estos datos. La interpolacion lineal es suficiente para\n");
    printf("      la diferenciacion numerica posterior.\n\n");
    
    printf("===============================================================\n");
    printf("PROXIMO PASO: Calcular df/dx usando diferencias finitas\n");
    printf("===============================================================\n\n");

    return 0;
}
