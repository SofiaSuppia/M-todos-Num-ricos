#include <stdio.h>
#include <math.h>

#define NP 11

/*
 * ============================================================================
 * PROBLEMA 2 - SEGUNDO PARCIAL 2024
 * ============================================================================
 *
 * Ecuacion diferencial ordinaria:
 *
 *      dy/dx = x * exp(x^2) / y,      y(0) = 1.
 *
 * Solucion exacta:
 *
 *      y(x) = exp(x^2 / 2).
 *
 * Metodo propuesto por el examen (Adams-Bashforth de 2 pasos):
 *
 *      y[n+1] = y[n] + h/2 * (3*f(x[n], y[n]) - f(x[n-1], y[n-1]))
 *
 * para n >= 1. Para n = 0 se usa el dato inicial y para calcular y[1] se usa
 * un metodo de un solo paso. En este programa se usa Euler:
 *
 *      y[1] = y[0] + h*f(x[0], y[0]).
 *
 * El programa resuelve:
 *
 *      a) Implementa el metodo propuesto.
 *      b) Usa la solucion exacta para comparar.
 *      c) Calcula la solucion aproximada con h = 0.1 y muestra y(0.5), y(1).
 *      d) Repite con RK4.
 *      e) Guarda los datos para graficar las tres soluciones.
 *
 * La estructura de Euler/RK4 se adapta de:
 *
 *      8. EDO/euler_neum_midpoint_RK4.c++
 * ============================================================================
 */

/*
 * Funcion f(x,y) de la EDO.
 */
double f(double x, double y)
{
    return x * exp(x * x) / y;
}

/*
 * Solucion exacta de la EDO.
 *
 * Separando variables:
 *
 *      y dy = x exp(x^2) dx
 *
 *      y^2/2 = exp(x^2)/2 + C
 *
 * y(0)=1 implica C=0, por lo tanto:
 *
 *      y = exp(x^2/2).
 */
double y_exacta(double x)
{
    return exp(x * x / 2.0);
}

/*
 * Un paso del metodo de Euler.
 */
double paso_euler(double x, double y, double h)
{
    return y + h * f(x, y);
}

/*
 * Un paso del metodo de Adams-Bashforth de 2 pasos.
 *
 * Necesita los dos valores anteriores: y[n] e y[n-1].
 */
double paso_ab2(double x_n, double y_n, double x_nm1, double y_nm1, double h)
{
    return y_n + (h / 2.0) * (3.0 * f(x_n, y_n) - f(x_nm1, y_nm1));
}

/*
 * Un paso completo de Runge-Kutta de cuarto orden (RK4).
 */
double paso_rk4(double x, double y, double h)
{
    double k1 = f(x, y);
    double k2 = f(x + h / 2.0, y + h * k1 / 2.0);
    double k3 = f(x + h / 2.0, y + h * k2 / 2.0);
    double k4 = f(x + h, y + h * k3);

    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

int main(void)
{
    const double h = 0.1;
    const double x0 = 0.0;

    double x[NP];
    double y_ab2[NP];
    double y_rk4[NP];
    double y_ex[NP];
    double error_ab2[NP];
    double error_rk4[NP];

    FILE *archivo;
    int i;

    /* Grilla de puntos x = 0.0, 0.1, ..., 1.0. */
    for (i = 0; i < NP; i++) {
        x[i] = x0 + i * h;
    }

    /*
     * ------------------------------------------------------------------
     * Condiciones iniciales.
     * ------------------------------------------------------------------
     */
    y_ab2[0] = 1.0;
    y_rk4[0] = 1.0;

    /*
     * Para el metodo AB2 se necesita y[1]. Se calcula con Euler, como indica
     * el enunciado.
     */
    y_ab2[1] = paso_euler(x[0], y_ab2[0], h);

    /*
     * ------------------------------------------------------------------
     * Metodo propuesto: Adams-Bashforth de 2 pasos.
     * ------------------------------------------------------------------
     */
    for (i = 1; i < NP - 1; i++) {
        y_ab2[i + 1] = paso_ab2(x[i], y_ab2[i], x[i - 1], y_ab2[i - 1], h);
    }

    /*
     * ------------------------------------------------------------------
     * Metodo RK4.
     * ------------------------------------------------------------------
     */
    for (i = 0; i < NP - 1; i++) {
        y_rk4[i + 1] = paso_rk4(x[i], y_rk4[i], h);
    }

    /*
     * ------------------------------------------------------------------
     * Solucion exacta y errores absolutos.
     * ------------------------------------------------------------------
     */
    for (i = 0; i < NP; i++) {
        y_ex[i] = y_exacta(x[i]);
        error_ab2[i] = fabs(y_ab2[i] - y_ex[i]);
        error_rk4[i] = fabs(y_rk4[i] - y_ex[i]);
    }

    /*
     * ------------------------------------------------------------------
     * Salida por pantalla.
     * ------------------------------------------------------------------
     */
    printf("========================================================================\n");
    printf("PROBLEMA 2 - SEGUNDO PARCIAL 2024\n");
    printf("========================================================================\n");
    printf("EDO: dy/dx = x*exp(x^2)/y,   y(0) = 1\n");
    printf("Solucion exacta: y(x) = exp(x^2/2)\n");
    printf("Paso h = %.1f\n\n", h);

    printf("Metodo propuesto (Adams-Bashforth 2 pasos):\n");
    printf("  y[n+1] = y[n] + h/2 * (3*f(x[n],y[n]) - f(x[n-1],y[n-1]))\n");
    printf("  y[1] se calcula con Euler: y[1] = y[0] + h*f(x[0],y[0])\n\n");

    printf("Tabla de resultados:\n");
    printf("    i        x          y exacta       y AB2          error AB2      y RK4          error RK4\n");

    for (i = 0; i < NP; i++) {
        printf("  %3d   %8.4f   %13.10f   %13.10f   %12.6e   %13.10f   %12.6e\n",
               i, x[i], y_ex[i], y_ab2[i], error_ab2[i], y_rk4[i], error_rk4[i]);
    }

    /*
     * ------------------------------------------------------------------
     * Inciso c): valores pedidos en x = 0.5 y x = 1.0.
     * ------------------------------------------------------------------
     */
    printf("\n------------------ INCISO c) VALORES PEDIDOS ---------------------------\n");
    printf("En x = 0.5:\n");
    printf("  y exacta     = %.10f\n", y_ex[5]);
    printf("  y AB2        = %.10f\n", y_ab2[5]);
    printf("  error AB2    = %.10e\n", error_ab2[5]);

    printf("\nEn x = 1.0:\n");
    printf("  y exacta     = %.10f\n", y_ex[10]);
    printf("  y AB2        = %.10f\n", y_ab2[10]);
    printf("  error AB2    = %.10e\n", error_ab2[10]);

    /*
     * ------------------------------------------------------------------
     * Inciso d): comparacion con RK4.
     * ------------------------------------------------------------------
     */
    printf("\n------------------ INCISO d) COMPARACION CON RK4 -----------------------\n");
    printf("En x = 0.5:\n");
    printf("  y RK4        = %.10f\n", y_rk4[5]);
    printf("  error RK4    = %.10e\n", error_rk4[5]);

    printf("\nEn x = 1.0:\n");
    printf("  y RK4        = %.10f\n", y_rk4[10]);
    printf("  error RK4    = %.10e\n", error_rk4[10]);

    /*
     * ------------------------------------------------------------------
     * Inciso e): archivo CSV para el grafico de las tres soluciones.
     * ------------------------------------------------------------------
     */
    archivo = fopen("resultados_problema2.csv", "w");

    if (archivo == NULL) {
        printf("\nError: no se pudo crear resultados_problema2.csv\n");
        return 1;
    }

    fprintf(archivo, "x,y_exacta,y_ab2,y_rk4,error_ab2,error_rk4\n");

    for (i = 0; i < NP; i++) {
        fprintf(archivo, "%.10f,%.10f,%.10f,%.10f,%.10e,%.10e\n",
                x[i], y_ex[i], y_ab2[i], y_rk4[i], error_ab2[i], error_rk4[i]);
    }

    fclose(archivo);

    printf("\nDatos guardados en resultados_problema2.csv\n");
    printf("Use grafico_problema2.py para generar el grafico de las tres soluciones.\n");

    return 0;
}
