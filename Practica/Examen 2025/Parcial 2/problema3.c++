#include <stdio.h>
#include <math.h>

/*
 * ============================================================================
 * PROBLEMA 3 - SEGUNDO PARCIAL 2025
 * ============================================================================
 *
 * EDO:
 *
 *      dy/dx = y - x^2 + 1
 *      y(0)  = 0.5
 *      x en [0, 2]
 *
 * Solucion exacta dada:
 *
 *      y(x) = (x + 1)^2 - 0.5 * e^x
 *
 * Incisos:
 *
 *      a) Euler y punto medio a mano, primer paso con h = 0.2.
 *      b) RK4 3/8 con h = 0.1; informar y(1.5) e y(2.0).
 *      c) Generar datos para graficar solucion exacta y aproximada.
 *
 * El esquema RK4 3/8 usado es:
 *
 *      k1 = f(x_i, y_i)
 *      k2 = f(x_i + h/3, y_i + h*k1/3)
 *      k3 = f(x_i + 2h/3, y_i - h*k1/3 + h*k2)
 *      k4 = f(x_i + h, y_i + h*k1 - h*k2 + h*k3)
 *
 *      y_{i+1} = y_i + (h/8)*(k1 + 3k2 + 3k3 + k4)
 * ============================================================================
 */

/* Funcion de la EDO. */
double f(double x, double y)
{
    return y - x * x + 1.0;
}

/* Solucion exacta. */
double y_exacta(double x)
{
    return (x + 1.0) * (x + 1.0) - 0.5 * exp(x);
}

int main(void)
{
    const double y0 = 0.5;
    const double h = 0.1;
    const int pasos = 20;

    double x;
    double y;
    double k1;
    double k2;
    double k3;
    double k4;
    double y_nuevo;
    double error;
    double y_15 = 0.0;
    double y_20 = 0.0;
    int i;

    FILE *archivo;

    printf("========================================================================\n");
    printf("PROBLEMA 3 - SEGUNDO PARCIAL 2025\n");
    printf("========================================================================\n");
    printf("dy/dx = y - x^2 + 1,  y(0) = 0.5,  x en [0,2]\n");
    printf("y(x) = (x+1)^2 - 0.5*e^x\n\n");

    printf("------------------ INCISO a) PRIMER PASO CON h = 0.2 ------------------\n");

    /* Euler */
    double k_euler = f(0.0, y0);
    double y_euler = y0 + 0.2 * k_euler;
    double y_ex_02 = y_exacta(0.2);

    printf("Euler:\n");
    printf("  k1 = f(0, 0.5) = %.10f\n", k_euler);
    printf("  y(0.2) aprox. = %.10f\n", y_euler);
    printf("  y(0.2) exacta = %.10f\n", y_ex_02);
    printf("  error absoluto = %.10f\n\n", fabs(y_euler - y_ex_02));

    /* Punto medio */
    double k1_medio = f(0.0, y0);
    double x_medio = 0.0 + 0.2 / 2.0;
    double y_medio_pred = y0 + (0.2 / 2.0) * k1_medio;
    double k2_medio = f(x_medio, y_medio_pred);
    double y_medio = y0 + 0.2 * k2_medio;

    printf("Punto medio:\n");
    printf("  k1 = f(0, 0.5) = %.10f\n", k1_medio);
    printf("  x_pred = 0.1, y_pred = %.10f\n", y_medio_pred);
    printf("  k2 = f(0.1, %.10f) = %.10f\n", y_medio_pred, k2_medio);
    printf("  y(0.2) aprox. = %.10f\n", y_medio);
    printf("  y(0.2) exacta = %.10f\n", y_ex_02);
    printf("  error absoluto = %.10f\n\n", fabs(y_medio - y_ex_02));

    printf("------------------ INCISO b) RK4 3/8 CON h = 0.1 ---------------------\n");
    printf("  i      x_i        y_i            k1            k2            k3            k4         y_(i+1)       exacta        error\n");

    archivo = fopen("resultados_rk4.csv", "w");
    if (archivo == NULL) {
        printf("ERROR: no se pudo crear resultados_rk4.csv\n");
        return 1;
    }
    fprintf(archivo, "x,y_numerica,y_exacta,error\n");

    x = 0.0;
    y = y0;
    fprintf(archivo, "%.10f,%.10f,%.10f,%.10f\n",
            x, y, y_exacta(x), fabs(y - y_exacta(x)));

    for (i = 1; i <= pasos; i++) {
        k1 = f(x, y);
        k2 = f(x + h / 3.0, y + h * k1 / 3.0);
        k3 = f(x + 2.0 * h / 3.0, y + h * (-k1 / 3.0 + k2));
        k4 = f(x + h, y + h * (k1 - k2 + k3));

        y_nuevo = y + (h / 8.0) * (k1 + 3.0 * k2 + 3.0 * k3 + k4);
        error = fabs(y_nuevo - y_exacta(x + h));

        if (i == 15) {
            y_15 = y_nuevo;
        }
        if (i == 20) {
            y_20 = y_nuevo;
        }

        printf("%3d  %8.4f  %12.10f  %12.10f  %12.10f  %12.10f  %12.10f  %12.10f  %12.10f  %12.10f\n",
               i, x, y, k1, k2, k3, k4, y_nuevo, y_exacta(x + h), error);

        fprintf(archivo, "%.10f,%.10f,%.10f,%.10f\n",
                x + h, y_nuevo, y_exacta(x + h), error);

        x += h;
        y = y_nuevo;
    }

    fclose(archivo);

    printf("\nResultados pedidos:\n");
    printf("  y(1.5) aproximada = %.10f\n", y_15);
    printf("  y(1.5) exacta     = %.10f\n", y_exacta(1.5));
    printf("  error absoluto    = %.10f\n", fabs(y_15 - y_exacta(1.5)));
    printf("  y(2.0) aproximada = %.10f\n", y_20);
    printf("  y(2.0) exacta     = %.10f\n", y_exacta(2.0));
    printf("  error absoluto    = %.10f\n", fabs(y_20 - y_exacta(2.0)));

    printf("\nSe genero el archivo resultados_rk4.csv para el grafico.\n");

    return 0;
}
