#include <stdio.h>
#include <math.h>

#define N 7

/*
 * ============================================================================
 * PROBLEMA 1 - SEGUNDO PARCIAL 2025
 * ============================================================================
 *
 * Integral exacta:
 *
 *      I = integral de 0 a 1 de e^x * sin(x) dx
 *
 * Valor exacto dado por el examen:
 *
 *      I_exacta = 0.9093306736
 *
 * El programa resuelve:
 *
 *      a) Cuadratura de Gauss-Legendre de 2 puntos.
 *      b) Simpson 1/3 simple.
 *      c) Tabla de 5 puntos equiespaciados usando spline cubico natural.
 *      d) Simpson compuesto usando la tabla del inciso c.
 *
 * Para el inciso c se usan los 7 pares de datos dados en el examen. El spline
 * cubico natural se construye con las segundas derivadas M[i] = S''(X[i]).
 *
 * En el inciso d se redondea la tabla del spline a 4 cifras decimales, porque
 * el examen pide integrar "a partir de la tabla del inciso c".
 * ============================================================================
 */

/* Funcion del problema 1. */
double f(double x)
{
    return exp(x) * sin(x);
}

/*
 * Resuelve un sistema lineal A*x = b por eliminacion gaussiana con pivoteo
 * parcial. Se usa para calcular las segundas derivadas del spline natural.
 */
void gauss(int n, double A[N][N], double b[N], double x[N])
{
    int i, j, k, pivote;

    for (k = 0; k < n - 1; k++) {
        pivote = k;
        for (i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[pivote][k])) {
                pivote = i;
            }
        }

        if (pivote != k) {
            for (j = k; j < n; j++) {
                double temp = A[k][j];
                A[k][j] = A[pivote][j];
                A[pivote][j] = temp;
            }
            double temp = b[k];
            b[k] = b[pivote];
            b[pivote] = temp;
        }

        for (i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (i = n - 1; i >= 0; i--) {
        double suma = b[i];
        for (j = i + 1; j < n; j++) {
            suma -= A[i][j] * x[j];
        }
        x[i] = suma / A[i][i];
    }
}

/*
 * Calcula las segundas derivadas M[i] del spline cubico natural que pasa por
 * los puntos (X[i], Y[i]).
 *
 * Para i = 1, ..., N-2 se plantea:
 *
 *   h[i-1]*M[i-1] + 2*(h[i-1]+h[i])*M[i] + h[i]*M[i+1]
 *
 *      = 6 * [ (Y[i+1]-Y[i])/h[i] - (Y[i]-Y[i-1])/h[i-1] ]
 *
 * con M[0] = M[N-1] = 0 por ser spline natural.
 */
void calcular_segundas_derivadas(double X[N], double Y[N], double M[N])
{
    double h[N - 1];
    double A[N][N];
    double b[N];
    double solucion[N];
    int i, j, m;

    for (i = 0; i < N - 1; i++) {
        h[i] = X[i + 1] - X[i];
    }

    for (i = 0; i < N; i++) {
        M[i] = 0.0;
        for (j = 0; j < N; j++) {
            A[i][j] = 0.0;
        }
        b[i] = 0.0;
    }

    m = N - 2;
    for (i = 1; i < N - 1; i++) {
        int fila = i - 1;

        if (i - 1 >= 1) {
            A[fila][fila - 1] = h[i - 1];
        }

        A[fila][fila] = 2.0 * (h[i - 1] + h[i]);

        if (i + 1 <= N - 2) {
            A[fila][fila + 1] = h[i];
        }

        b[fila] = 6.0 * ((Y[i + 1] - Y[i]) / h[i]
                       - (Y[i] - Y[i - 1]) / h[i - 1]);
    }

    gauss(m, A, b, solucion);

    for (i = 0; i < m; i++) {
        M[i + 1] = solucion[i];
    }
}

/*
 * Evalua el spline cubico natural en un punto x dentro del intervalo.
 *
 * Para x en [X[i], X[i+1]]:
 *
 *   S(x) = (X[i+1]-x)/h * Y[i] + (x-X[i])/h * Y[i+1]
 *
 *          + [ ((X[i+1]-x)^3/h - h*(X[i+1]-x)) * M[i]
 *            + ((x-X[i])^3/h - h*(x-X[i])) * M[i+1] ] / 6
 */
double evaluar_spline(double X[N], double Y[N], double M[N], double x)
{
    int i;
    double h, a, b;

    for (i = 0; i < N - 1; i++) {
        if (x >= X[i] && x <= X[i + 1]) {
            break;
        }
    }

    if (i == N - 1) {
        i = N - 2;
    }

    h = X[i + 1] - X[i];
    a = (X[i + 1] - x) / h;
    b = (x - X[i]) / h;

    return a * Y[i] + b * Y[i + 1]
         + ((a * a * a - a) * M[i] + (b * b * b - b) * M[i + 1]) * h * h / 6.0;
}

/* Redondea a 4 cifras decimales, para construir la tabla pedida. */
double redondear4(double x)
{
    return round(x * 10000.0) / 10000.0;
}

int main(void)
{
    const double I_exacta = 0.9093306736;

    double X[N] = {0.00, 0.20, 0.43, 0.55, 0.85, 0.90, 1.00};
    double Y[N] = {0.0000, 0.2426, 0.6408, 0.9059, 1.7577, 1.9267, 2.2873};
    double M[N];
    double x_nuevos[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    double y_nuevos[5];
    double y_tabla[5];
    int i;

    double raiz3 = 1.0 / sqrt(3.0);
    double x1 = 0.5 * (1.0 - raiz3);
    double x2 = 0.5 * (1.0 + raiz3);
    double I_gauss = 0.5 * (f(x1) + f(x2));
    double I_simpson_simple = (1.0 / 6.0) * (f(0.0) + 4.0 * f(0.5) + f(1.0));

    double h_simpson = 0.25;
    double I_simpson_compuesto;

    printf("========================================================================\n");
    printf("PROBLEMA 1 - SEGUNDO PARCIAL 2025\n");
    printf("========================================================================\n");
    printf("I = integral de 0 a 1 de e^x * sin(x) dx\n");
    printf("Valor exacto dado: %.10f\n\n", I_exacta);

    printf("------------------ INCISO a) GAUSS-LEGENDRE 2 PUNTOS ------------------\n");
    printf("x1 = %.10f\n", x1);
    printf("x2 = %.10f\n", x2);
    printf("f(x1) = %.10f\n", f(x1));
    printf("f(x2) = %.10f\n", f(x2));
    printf("I_Gauss = %.10f\n", I_gauss);
    printf("I_Gauss (6 decimales) = %.6f\n", I_gauss);
    printf("Error absoluto exacto = %.10f\n\n", fabs(I_gauss - I_exacta));

    printf("------------------ INCISO b) SIMPSON 1/3 SIMPLE -----------------------\n");
    printf("f(0.0) = %.10f\n", f(0.0));
    printf("f(0.5) = %.10f\n", f(0.5));
    printf("f(1.0) = %.10f\n", f(1.0));
    printf("I_Simpson_simple = %.10f\n", I_simpson_simple);
    printf("I_Simpson_simple (6 decimales) = %.6f\n", I_simpson_simple);
    printf("Error absoluto exacto = %.10f\n\n", fabs(I_simpson_simple - I_exacta));

    printf("------------------ INCISO c) SPLINE CUBICO NATURAL --------------------\n");
    calcular_segundas_derivadas(X, Y, M);

    printf("Segundas derivadas M[i] = S''(X[i]):\n");
    for (i = 0; i < N; i++) {
        printf("  M[%d] = %.10f\n", i, M[i]);
    }

    printf("\nTabla original:\n");
    printf("    i        X[i]         Y[i]\n");
    for (i = 0; i < N; i++) {
        printf("  %3d   %10.4f   %10.4f\n", i, X[i], Y[i]);
    }

    printf("\nTabla nueva de 5 datos equiespaciados (4 decimales):\n");
    printf("    i        X[i]         Y[i]\n");
    for (i = 0; i < 5; i++) {
        y_nuevos[i] = evaluar_spline(X, Y, M, x_nuevos[i]);
        y_tabla[i] = redondear4(y_nuevos[i]);
        printf("  %3d   %10.4f   %10.4f\n", i, x_nuevos[i], y_tabla[i]);
    }

    printf("\n------------------ INCISO d) SIMPSON COMPUESTO ------------------------\n");
    printf("Se usa la tabla redondeada del inciso c, con h = 0.25:\n");
    I_simpson_compuesto = (h_simpson / 3.0) *
        (y_tabla[0] + 4.0 * y_tabla[1] + 2.0 * y_tabla[2]
         + 4.0 * y_tabla[3] + y_tabla[4]);

    printf("I_Simpson_compuesto = %.10f\n", I_simpson_compuesto);
    printf("Error absoluto exacto = %.10f\n\n",
           fabs(I_simpson_compuesto - I_exacta));

    printf("Comparacion con el valor exacto:\n");
    printf("  I_exacta              = %.10f\n", I_exacta);
    printf("  Gauss-Legendre        = %.10f\n", I_gauss);
    printf("  Simpson 1/3 simple    = %.10f\n", I_simpson_simple);
    printf("  Simpson compuesto     = %.10f\n", I_simpson_compuesto);

    return 0;
}
