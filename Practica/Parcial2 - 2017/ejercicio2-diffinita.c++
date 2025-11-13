/*
 * ===========================================================================
 * PARCIAL 2 - 2017: EJERCICIO 2 - PASO 2
 * Calculo de df/dx usando diferencias finitas centradas
 * ===========================================================================
 * 
 * PASO 2: Calcular df/dx en los 11 puntos equiespaciados
 * Metodo: Diferencias finitas centradas (error O(h^2))
 * 
 * Entrada: 11 puntos equiespaciados del PASO 1 (h = 0.1)
 * Salida: Valores de df/dx en cada punto
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_PUNTOS 11  // 11 puntos equiespaciados

// ===========================================================================
// DATOS DE ENTRADA (RESULTADOS DEL PASO 1)
// ===========================================================================

// Estos son los 11 puntos equiespaciados obtenidos por interpolacion lineal
double X[N_PUNTOS] = {
    1.0000, 1.1000, 1.2000, 1.3000, 1.4000,
    1.5000, 1.6000, 1.7000, 1.8000, 1.9000, 2.0000
};

double f_X[N_PUNTOS] = {
    0.0980000000,   // x = 1.0
   -0.0430000000,   // x = 1.1
   -0.2030000000,   // x = 1.2
   -0.3688333333,   // x = 1.3
   -0.5330000000,   // x = 1.4
   -0.6930000000,   // x = 1.5
   -0.8190000000,   // x = 1.6
   -0.9450000000,   // x = 1.7
   -1.0112500000,   // x = 1.8
   -1.0515714286,   // x = 1.9
   -1.0530000000    // x = 2.0
};

// ===========================================================================
// FUNCIONES PARA DIFERENCIACION NUMERICA
// ===========================================================================

/**
 * Calcula la primera derivada usando diferencias finitas
 * 
 * Formulas:
 * - Extremo izquierdo (forward): f'(x0) = [f(x1) - f(x0)] / h  [O(h)]
 * - Puntos interiores (centrada): f'(xi) = [f(xi+1) - f(xi-1)] / (2h)  [O(h^2)]
 * - Extremo derecho (backward): f'(xn) = [f(xn) - f(xn-1)] / h  [O(h)]
 */
void calcular_primera_derivada(double X[], double f[], int n, double h, double fp[]) {
    printf("\n===============================================================\n");
    printf("CALCULANDO PRIMERA DERIVADA CON DIFERENCIAS FINITAS\n");
    printf("===============================================================\n\n");
    
    printf("- Numero de puntos: %d\n", n);
    printf("- Paso h = %.4f\n", h);
    printf("- Metodo:\n");
    printf("  * Extremos: Diferencias forward/backward (error O(h))\n");
    printf("  * Interiores: Diferencias centradas (error O(h^2))\n\n");
    
    // Extremo izquierdo: diferencia forward
    fp[0] = (f[1] - f[0]) / h;
    printf("Punto 0 (x=%.4f): f'(x) = [f(%.4f) - f(%.4f)] / %.4f\n",
           X[0], X[1], X[0], h);
    printf("                         = [%.10f - %.10f] / %.4f\n",
           f[1], f[0], h);
    printf("                         = %.10f (forward)\n\n", fp[0]);
    
    // Puntos interiores: diferencias centradas
    printf("Puntos interiores (diferencias centradas):\n");
    for (int i = 1; i < n-1; i++) {
        fp[i] = (f[i+1] - f[i-1]) / (2.0 * h);
        if (i <= 2) {  // Mostrar detalles solo para los primeros puntos
            printf("Punto %d (x=%.4f): f'(x) = [f(%.4f) - f(%.4f)] / (2*%.4f)\n",
                   i, X[i], X[i+1], X[i-1], h);
            printf("                          = [%.10f - %.10f] / %.4f\n",
                   f[i+1], f[i-1], 2.0*h);
            printf("                          = %.10f\n\n", fp[i]);
        }
    }
    if (n > 4) printf("... (calculo similar para puntos 3-%d)\n\n", n-2);
    
    // Extremo derecho: diferencia backward
    fp[n-1] = (f[n-1] - f[n-2]) / h;
    printf("Punto %d (x=%.4f): f'(x) = [f(%.4f) - f(%.4f)] / %.4f\n",
           n-1, X[n-1], X[n-1], X[n-2], h);
    printf("                          = [%.10f - %.10f] / %.4f\n",
           f[n-1], f[n-2], h);
    printf("                          = %.10f (backward)\n\n", fp[n-1]);
    
    printf("OK - Primera derivada calculada en todos los puntos\n");
}

// ===========================================================================
// PROGRAMA PRINCIPAL
// ===========================================================================

int main() {
    double fp[N_PUNTOS];  // Array para almacenar df/dx
    double h = 0.1;       // Paso (equiespaciado)
    
    printf("===============================================================\n");
    printf("       PARCIAL 2 - 2017: EJERCICIO 2 - PASO 2              \n");
    printf("  Diferenciacion Numerica: Calculo de df/dx                 \n");
    printf("===============================================================\n\n");

    // Mostrar tabla de entrada
    printf("===============================================================\n");
    printf("TABLA DE ENTRADA (11 puntos equiespaciados del PASO 1)\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-----------------+\n");
    printf("|  i  |    x     |      f(x)       |\n");
    printf("+-----+----------+-----------------+\n");
    for (int i = 0; i < N_PUNTOS; i++) {
        printf("| %3d | %8.4f | %15.10f |\n", i, X[i], f_X[i]);
    }
    printf("+-----+----------+-----------------+\n\n");
    
    // Calcular primera derivada
    calcular_primera_derivada(X, f_X, N_PUNTOS, h, fp);
    
    // Mostrar tabla de resultados
    printf("\n===============================================================\n");
    printf("TABLA DE RESULTADOS: f(x) y df/dx\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-----------------+-----------------+\n");
    printf("|  i  |    x     |      f(x)       |     df/dx       |\n");
    printf("+-----+----------+-----------------+-----------------+\n");
    for (int i = 0; i < N_PUNTOS; i++) {
        printf("| %3d | %8.4f | %15.10f | %15.10f |\n", 
               i, X[i], f_X[i], fp[i]);
    }
    printf("+-----+----------+-----------------+-----------------+\n\n");
    
    // Guardar resultados en archivo para el PASO 3
    FILE *archivo = fopen("paso2_resultados.txt", "w");
    if (archivo != NULL) {
        fprintf(archivo, "# Resultados del PASO 2: Diferenciacion numerica\n");
        fprintf(archivo, "# i    x         f(x)              df/dx\n");
        for (int i = 0; i < N_PUNTOS; i++) {
            fprintf(archivo, "%3d %.4f %.15f %.15f\n", i, X[i], f_X[i], fp[i]);
        }
        fclose(archivo);
        printf("OK - Resultados guardados en 'paso2_resultados.txt'\n\n");
    }
    
    // Resumen
    printf("===============================================================\n");
    printf("  RESUMEN DEL PASO 2                                        \n");
    printf("===============================================================\n\n");
    
    printf("OK ENTRADA: 11 puntos equiespaciados (h = 0.1)\n");
    printf("OK METODO: Diferencias finitas\n");
    printf("   - Extremos: Forward/Backward (error O(h))\n");
    printf("   - Interiores: Centradas (error O(h^2))\n");
    printf("OK SALIDA: Valores de df/dx en 11 puntos\n");
    printf("   - Valores listos para integracion\n");
    printf("   - Archivo: paso2_resultados.txt\n\n");
    
    // Estadisticas de la derivada
    double min_fp = fp[0], max_fp = fp[0];
    for (int i = 1; i < N_PUNTOS; i++) {
        if (fp[i] < min_fp) min_fp = fp[i];
        if (fp[i] > max_fp) max_fp = fp[i];
    }
    
    printf("Estadisticas de df/dx:\n");
    printf("   - Minimo: %.10f (en x=%.4f)\n", min_fp, X[0]);
    printf("   - Maximo: %.10f (en x=%.4f)\n", max_fp, X[N_PUNTOS-1]);
    printf("   - Rango: [%.10f, %.10f]\n\n", min_fp, max_fp);
    
    printf("===============================================================\n");
    printf("PROXIMO PASO: Integrar df/dx con Trapecio y Simpson 1/3\n");
    printf("===============================================================\n\n");

    return 0;
}
