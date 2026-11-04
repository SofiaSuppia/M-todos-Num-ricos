#include <stdio.h>
#include <math.h>

/*
 * ============================================================================
 * PROBLEMA 2 - SEGUNDO PARCIAL 2025
 * ============================================================================
 *
 * Sea:
 *
 *      f(x) = sqrt(1 + x^2)
 *
 * Se pide calcular f'(1) y completar una tabla de derivadas aproximadas en
 * x = 0, 0.1, 0.2, ..., 1.0.
 *
 * Derivada exacta:
 *
 *      f'(x) = x / sqrt(1 + x^2)
 *      f'(1) = 1 / sqrt(2) = 0.7071067812
 *
 * Incisos:
 *
 *      b) Diferencia finita centrada de orden h^2 con h = 0.5.
 *      c) La misma formula con h = 0.01.
 *      d) Tabla con h = 0.1:
 *             primer punto: diferencia hacia adelante,
 *             ultimo punto: diferencia hacia atras,
 *             puntos internos: diferencia centrada.
 *      e) Trapecio compuesto usando la tabla del inciso d.
 *
 * En el inciso e se usan los valores de la tabla redondeados a 4 decimales,
 * tal como se escriben en papel.
 * ============================================================================
 */

/* Funcion f(x) = sqrt(1 + x^2). */
double f(double x)
{
    return sqrt(1.0 + x * x);
}

/* Derivada exacta f'(x) = x / sqrt(1 + x^2). */
double f_prima_exacta(double x)
{
    return x / sqrt(1.0 + x * x);
}

/* Redondea a 4 cifras decimales, para la tabla pedida en el inciso d. */
double redondear4(double x)
{
    return round(x * 10000.0) / 10000.0;
}

int main(void)
{
    const double exacta_en_1 = 1.0 / sqrt(2.0);
    const double h_centrada_1 = 0.5;
    const double h_centrada_2 = 0.01;
    const double h_tabla = 0.1;

    double centrada_1;
    double centrada_2;
    double error_1;
    double error_2;

    double x[11];
    double derivada[11];
    double derivada_tabla[11];
    double suma_interior;
    double I_trapecio;
    double I_exacta;
    double error_trapecio;
    int i;

    centrada_1 = (f(1.0 + h_centrada_1) - f(1.0 - h_centrada_1))
               / (2.0 * h_centrada_1);
    centrada_2 = (f(1.0 + h_centrada_2) - f(1.0 - h_centrada_2))
               / (2.0 * h_centrada_2);

    error_1 = fabs(centrada_1 - exacta_en_1);
    error_2 = fabs(centrada_2 - exacta_en_1);

    printf("========================================================================\n");
    printf("PROBLEMA 2 - SEGUNDO PARCIAL 2025\n");
    printf("========================================================================\n");
    printf("f(x)  = sqrt(1 + x^2)\n");
    printf("f'(x) = x / sqrt(1 + x^2)\n");
    printf("f'(1) = 1/sqrt(2) = %.10f\n\n", exacta_en_1);

    printf("------------------ INCISO b) CENTRADA CON h = 0.5 ---------------------\n");
    printf("f(1.5) = %.10f\n", f(1.5));
    printf("f(0.5) = %.10f\n", f(0.5));
    printf("f'(1) aproximada = %.10f\n", centrada_1);
    printf("f'(1) aproximada (6 decimales) = %.6f\n", centrada_1);
    printf("Error absoluto exacto = %.10f\n\n", error_1);

    printf("------------------ INCISO c) CENTRADA CON h = 0.01 --------------------\n");
    printf("f(1.01) = %.10f\n", f(1.01));
    printf("f(0.99) = %.10f\n", f(0.99));
    printf("f'(1) aproximada = %.10f\n", centrada_2);
    printf("f'(1) aproximada (6 decimales) = %.6f\n", centrada_2);
    printf("Error absoluto exacto = %.10f\n", error_2);
    printf("Cociente de errores E(0.5)/E(0.01) = %.6f\n", error_1 / error_2);
    printf("Cociente esperado (0.5/0.01)^2 = %.6f\n\n", 2500.0);

    printf("------------------ INCISO d) TABLA CON h = 0.1 ------------------------\n");
    for (i = 0; i <= 10; i++) {
        x[i] = i * h_tabla;
    }

    for (i = 0; i <= 10; i++) {
        if (i == 0) {
            derivada[i] = (f(x[1]) - f(x[0])) / h_tabla;
        } else if (i == 10) {
            derivada[i] = (f(x[10]) - f(x[9])) / h_tabla;
        } else {
            derivada[i] = (f(x[i + 1]) - f(x[i - 1])) / (2.0 * h_tabla);
        }

        derivada_tabla[i] = redondear4(derivada[i]);
    }

    printf("    i        x[i]      f'(x) aprox.      f'(x) exacta\n");
    for (i = 0; i <= 10; i++) {
        printf("  %3d   %10.4f   %14.4f   %14.10f\n",
               i, x[i], derivada_tabla[i], f_prima_exacta(x[i]));
    }

    printf("\n------------------ INCISO e) TRAPECIO COMPUESTO -----------------------\n");
    suma_interior = 0.0;
    for (i = 1; i <= 9; i++) {
        suma_interior += derivada_tabla[i];
    }

    I_trapecio = (h_tabla / 2.0)
               * (derivada_tabla[0] + 2.0 * suma_interior + derivada_tabla[10]);

    I_exacta = f(1.0) - f(0.0);
    error_trapecio = fabs(I_trapecio - I_exacta);

    printf("Suma de puntos internos = %.10f\n", suma_interior);
    printf("I_Trapecio = %.10f\n", I_trapecio);
    printf("I_exacta   = f(1) - f(0) = sqrt(2) - 1 = %.10f\n", I_exacta);
    printf("Error absoluto exacto = %.10f\n\n", error_trapecio);

    printf("Resumen:\n");
    printf("  f'(1) exacta             = %.10f\n", exacta_en_1);
    printf("  Centrada h = 0.5         = %.10f  (error %.10f)\n",
           centrada_1, error_1);
    printf("  Centrada h = 0.01        = %.10f  (error %.10f)\n",
           centrada_2, error_2);
    printf("  Trapecio con tabla h=0.1 = %.10f  (error %.10f)\n",
           I_trapecio, error_trapecio);

    return 0;
}
