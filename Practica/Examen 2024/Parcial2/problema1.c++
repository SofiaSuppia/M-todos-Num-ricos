#include <stdio.h>
#include <math.h>

#define N 5

/*
 * ============================================================================
 * PROBLEMA 1 - SEGUNDO PARCIAL 2024
 * ============================================================================
 *
 * La tabla dada corresponde a la funcion
 *
 *      f(x) = ln(1 + x^2) - sin(x)
 *
 * con los cinco puntos:
 *
 *      x = 1.00, 1.20, 1.45, 1.78, 2.00
 *      y = -0.1483, -0.0400, 0.1395, 0.4493, 0.7001
 *
 * El programa resuelve:
 *
 *      a) Spline cubico natural y polinomio del tramo [1.20, 1.45].
 *      b) Tabla de 5 puntos equiespaciados en [1, 2] con h = 0.25.
 *      c) Tabla de derivadas por diferencias finitas de menor orden.
 *      d) Integral de f'(x) en [1, 2] por trapecio compuesto.
 *      e) La misma integral por Gauss-Legendre de 2 puntos.
 *      f) La misma integral por Gauss-Legendre de 6 puntos.
 *
 * La logica del spline y de la eliminacion gaussiana se adapta de los codigos
 * del repositorio:
 *
 *      - 5. InterpolacionSegmentaria/splines.c++
 *      - 5. InterpolacionSegmentaria/gauss.h
 *
 * La logica de Gauss-Legendre se adapta de:
 *
 *      - 7. Diferenciacion/Integracion por Cuadrados/gauss_legendre.c++
 *
 * Los datos de la tabla original estan redondeados a 4 decimales. Por eso, en
 * los incisos b) y c) se conserva esa precision para que el desarrollo a mano
 * coincida con la salida del programa.
 * ============================================================================
 */

/* Funcion exacta del enunciado. */
double f(double x)
{
    return log(1.0 + x * x) - sin(x);
}

/* Derivada exacta de f(x), usada en los incisos d), e) y f). */
double df(double x)
{
    return 2.0 * x / (1.0 + x * x) - cos(x);
}

/* Redondea a 4 cifras decimales, igual que los datos de la tabla. */
double redondear4(double x)
{
    return round(x * 10000.0) / 10000.0;
}

/*
 * Resuelve un sistema lineal A*x = b mediante eliminacion gaussiana con
 * pivoteo parcial.
 *
 * Se usa para calcular las segundas derivadas M[i] del spline cubico natural.
 * La implementacion sigue la idea de gauss.h, pero sin imprimir la matriz,
 * porque en este programa se necesita una salida limpia para la captura.
 */
void gauss(int n, double A[N][N], double b[N], double x[N])
{
    int i, j, k, pivote;

    for (k = 0; k < n - 1; k++) {
        /* Busca el pivote de mayor valor absoluto en la columna k. */
        pivote = k;
        for (i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > fabs(A[pivote][k])) {
                pivote = i;
            }
        }

        /* Intercambia la fila k con la fila del pivote. */
        if (pivote != k) {
            for (j = k; j < n; j++) {
                double aux = A[k][j];
                A[k][j] = A[pivote][j];
                A[pivote][j] = aux;
            }

            double aux = b[k];
            b[k] = b[pivote];
            b[pivote] = aux;
        }

        /* Eliminacion hacia adelante. */
        for (i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];

            for (j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }

            b[i] -= factor * b[k];
        }
    }

    /* Sustitucion hacia atras. */
    for (i = n - 1; i >= 0; i--) {
        double suma = b[i];

        for (j = i + 1; j < n; j++) {
            suma -= A[i][j] * x[j];
        }

        x[i] = suma / A[i][i];
    }
}

/*
 * Calcula las segundas derivadas M[i] = S''(X[i]) del spline cubico natural.
 *
 * Para los nodos interiores i = 1, ..., N-2 se plantea:
 *
 *   h[i-1]*M[i-1] + 2*(h[i-1]+h[i])*M[i] + h[i]*M[i+1]
 *
 *      = 6 * [ (Y[i+1]-Y[i])/h[i] - (Y[i]-Y[i-1])/h[i-1] ]
 *
 * y en los extremos se impone la condicion natural:
 *
 *      M[0] = 0,   M[N-1] = 0.
 */
void calcular_segundas_derivadas(const double X[N], const double Y[N], double M[N])
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
 * Evalua el spline cubico natural en un punto x.
 *
 * Si x pertenece al intervalo [X[i], X[i+1]]:
 *
 *   S(x) = (X[i+1]-x)/h * Y[i] + (x-X[i])/h * Y[i+1]
 *
 *          + [ ((X[i+1]-x)^3/h - h*(X[i+1]-x)) * M[i]
 *            + ((x-X[i])^3/h - h*(x-X[i])) * M[i+1] ] / 6
 */
double evaluar_spline(const double X[N], const double Y[N], const double M[N], double x)
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

/*
 * Obtiene los coeficientes del polinomio del tramo i.
 *
 * coef_local[0..3] corresponden a:
 *
 *      S(t) = c3*t^3 + c2*t^2 + c1*t + c0,   t = x - X[i]
 *
 * coef_global[0..3] corresponden a:
 *
 *      S(x) = c3*x^3 + c2*x^2 + c1*x + c0
 */
void coeficientes_tramo(const double X[N], const double Y[N], const double M[N],
                        int i, double coef_local[4], double coef_global[4])
{
    double h = X[i + 1] - X[i];
    double x1 = X[i];
    double y1 = Y[i];
    double y2 = Y[i + 1];
    double m1 = M[i];
    double m2 = M[i + 1];

    coef_local[0] = (m2 - m1) / (6.0 * h);
    coef_local[1] = m1 / 2.0;
    coef_local[2] = (y2 - y1) / h - h * (2.0 * m1 + m2) / 6.0;
    coef_local[3] = y1;

    coef_global[0] = coef_local[0];
    coef_global[1] = coef_local[1] - 3.0 * coef_local[0] * x1;
    coef_global[2] = coef_local[2] - 2.0 * coef_local[1] * x1
                   + 3.0 * coef_local[0] * x1 * x1;
    coef_global[3] = coef_local[3] - coef_local[2] * x1
                   + coef_local[1] * x1 * x1
                   - coef_local[0] * x1 * x1 * x1;
}

int main(void)
{
    const double X[N] = {1.00, 1.20, 1.45, 1.78, 2.00};
    const double Y[N] = {-0.1483, -0.0400, 0.1395, 0.4493, 0.7001};

    const double x_eq[N] = {1.00, 1.25, 1.50, 1.75, 2.00};
    const double h_eq = 0.25;

    double M[N];
    double y_spline[N];
    double y_tabla[N];
    double y_exacta[N];
    double error[N];
    double derivada[N];
    double derivada4[N];

    double coef_local[4];
    double coef_global[4];

    double I_exacta;
    double I_trapecio;
    double I_gauss2;
    double I_gauss6;

    int i;

    printf("========================================================================\n");
    printf("PROBLEMA 1 - SEGUNDO PARCIAL 2024\n");
    printf("========================================================================\n");
    printf("f(x) = ln(1 + x^2) - sin(x)\n\n");

    printf("Datos originales:\n");
    printf("    i        X[i]         Y[i]\n");
    for (i = 0; i < N; i++) {
        printf("  %3d   %10.4f   %10.4f\n", i, X[i], Y[i]);
    }

    /*
     * ------------------------------------------------------------------
     * INCISO a): spline cubico natural y polinomio del tramo [1.20,1.45]
     * ------------------------------------------------------------------
     */
    printf("\n------------------ INCISO a) SPLINE CUBICO NATURAL ---------------------\n");
    calcular_segundas_derivadas(X, Y, M);

    printf("Segundas derivadas M[i] = S''(X[i]):\n");
    for (i = 0; i < N; i++) {
        printf("  M[%d] = %12.8f\n", i, M[i]);
    }

    /*
     * El tramo [1.20, 1.45] es el intervalo i = 1.
     */
    coeficientes_tramo(X, Y, M, 1, coef_local, coef_global);

    printf("\nPolinomio del tramo [1.20, 1.45]:\n");
    printf("Con t = x - 1.20:\n");
    printf("  S(t) = %.10f*t^3 %+.10f*t^2 %+.10f*t %+.10f\n",
           coef_local[0], coef_local[1], coef_local[2], coef_local[3]);
    printf("En potencias de x:\n");
    printf("  S(x) = %.10f*x^3 %+.10f*x^2 %+.10f*x %+.10f\n",
           coef_global[0], coef_global[1], coef_global[2], coef_global[3]);

    printf("\nVerificacion en los extremos del tramo:\n");
    {
        double t1 = 1.20 - X[1];
        double t2 = 1.45 - X[1];
        double s1 = ((coef_local[0] * t1 + coef_local[1]) * t1 + coef_local[2]) * t1
                  + coef_local[3];
        double s2 = ((coef_local[0] * t2 + coef_local[1]) * t2 + coef_local[2]) * t2
                  + coef_local[3];

        printf("  S(1.20) = %.10f   (dato Y[1] = %.10f)\n", s1, Y[1]);
        printf("  S(1.45) = %.10f   (dato Y[2] = %.10f)\n", s2, Y[2]);
    }

    /*
     * ------------------------------------------------------------------
     * INCISO b): tabla de 5 puntos equiespaciados en [1,2]
     * ------------------------------------------------------------------
     */
    printf("\n------------------ INCISO b) TABLA EQUIESPACIADA -----------------------\n");
    printf("    i        x          f(x) exacta    S(x) spline     error exacto\n");

    for (i = 0; i < N; i++) {
        y_exacta[i] = f(x_eq[i]);
        y_spline[i] = evaluar_spline(X, Y, M, x_eq[i]);
        y_tabla[i] = redondear4(y_spline[i]);
        error[i] = fabs(y_tabla[i] - y_exacta[i]);

        printf("  %3d   %8.4f   %14.8f   %14.8f   %14.8e\n",
               i, x_eq[i], y_exacta[i], y_tabla[i], error[i]);
    }

    /*
     * ------------------------------------------------------------------
     * INCISO c): derivadas por diferencias finitas de menor orden
     * ------------------------------------------------------------------
     *
     * Se usan los valores redondeados del spline, porque son los que se
     * escriben en la tabla del inciso b) y los que se usarian a mano.
     *
     * Primer punto: diferencia hacia adelante.
     * Puntos interiores: diferencia centrada.
     * Ultimo punto: diferencia hacia atras.
     */
    printf("\n------------------ INCISO c) TABLA DE DERIVADAS -----------------------\n");

    derivada[0] = (y_tabla[1] - y_tabla[0]) / h_eq;
    derivada[N - 1] = (y_tabla[N - 1] - y_tabla[N - 2]) / h_eq;

    for (i = 1; i < N - 1; i++) {
        derivada[i] = (y_tabla[i + 1] - y_tabla[i - 1]) / (2.0 * h_eq);
    }

    printf("    i        x          S(x) tabla     f'(x) aproximada   f'(x) exacta\n");

    for (i = 0; i < N; i++) {
        derivada4[i] = redondear4(derivada[i]);
        printf("  %3d   %8.4f   %12.4f   %16.4f   %14.8f\n",
               i, x_eq[i], y_tabla[i], derivada4[i], df(x_eq[i]));
    }

    /*
     * ------------------------------------------------------------------
     * INCISO d): integral por trapecio compuesto usando la tabla c)
     * ------------------------------------------------------------------
     *
     * Integral pedida:
     *
     *      I = integral de 1 a 2 de f'(x) dx = f(2) - f(1).
     */
    printf("\n------------------ INCISO d) TRAPECIO COMPUESTO ------------------------\n");

    I_exacta = f(2.0) - f(1.0);

    I_trapecio = h_eq * (0.5 * derivada4[0] + derivada4[1] + derivada4[2]
                       + derivada4[3] + 0.5 * derivada4[4]);

    printf("Integral exacta: I = f(2) - f(1) = %.10f\n", I_exacta);
    printf("Integral por trapecio compuesto  = %.10f\n", I_trapecio);
    printf("Error absoluto exacto            = %.10e\n", fabs(I_trapecio - I_exacta));

    /*
     * ------------------------------------------------------------------
     * INCISO e): Gauss-Legendre de 2 puntos
     * ------------------------------------------------------------------
     *
     * Para integrar f'(x) en [1,2] se transforma el intervalo a [-1,1]:
     *
     *      x = 1.5 + 0.5*t,   dx = 0.5*dt
     *
     * Los nodos y pesos para 2 puntos son:
     *
     *      t1 = -1/sqrt(3),  t2 = 1/sqrt(3),  w1 = w2 = 1.
     */
    printf("\n------------------ INCISO e) GAUSS-LEGENDRE 2 PUNTOS ------------------\n");
    {
        double raiz3 = 1.0 / sqrt(3.0);
        double t1 = -raiz3;
        double t2 = raiz3;
        double x1 = 1.5 + 0.5 * t1;
        double x2 = 1.5 + 0.5 * t2;
        double g1 = df(x1);
        double g2 = df(x2);

        I_gauss2 = 0.5 * (g1 + g2);

        printf("Nodos de Gauss:\n");
        printf("  t1 = %.10f   ->   x1 = %.10f\n", t1, x1);
        printf("  t2 = %.10f   ->   x2 = %.10f\n", t2, x2);
        printf("Valores de f'(x):\n");
        printf("  f'(x1) = %.10f\n", g1);
        printf("  f'(x2) = %.10f\n", g2);
        printf("Integral por Gauss de 2 puntos = %.10f\n", I_gauss2);
        printf("Error absoluto exacto          = %.10e\n", fabs(I_gauss2 - I_exacta));
        printf("Error porcentual exacto        = %.10f %%\n",
               fabs(I_gauss2 - I_exacta) / fabs(I_exacta) * 100.0);
    }

    /*
     * ------------------------------------------------------------------
     * INCISO f): Gauss-Legendre de 6 puntos
     * ------------------------------------------------------------------
     *
     * Se usan los mismos nodos y pesos que en gauss_legendre.c++.
     * La transformacion al intervalo [1,2] es la misma que en el inciso e).
     */
    printf("\n------------------ INCISO f) GAUSS-LEGENDRE 6 PUNTOS ------------------\n");
    {
        const double t[6] = {-0.932469514, -0.661209386, -0.238619186,
                              0.238619186,  0.661209386,  0.932469514};
        const double w[6] = {0.1713245, 0.3607616, 0.4679139,
                             0.4679139, 0.3607616, 0.1713245};
        double suma = 0.0;

        printf("    i        t_i            x_i            f'(x_i)          peso w_i\n");

        for (i = 0; i < 6; i++) {
            double xi = 1.5 + 0.5 * t[i];
            double gi = df(xi);

            suma += w[i] * gi;

            printf("  %3d   %12.9f   %12.9f   %14.10f   %12.7f\n",
                   i, t[i], xi, gi, w[i]);
        }

        I_gauss6 = 0.5 * suma;

        printf("Integral por Gauss de 6 puntos = %.10f\n", I_gauss6);
        printf("Error absoluto exacto          = %.10e\n", fabs(I_gauss6 - I_exacta));
        printf("Error porcentual exacto        = %.10e %%\n",
               fabs(I_gauss6 - I_exacta) / fabs(I_exacta) * 100.0);
    }

    /*
     * Resumen final para comparar rapidamente los metodos.
     */
    printf("\n------------------ RESUMEN FINAL ---------------------------------------\n");
    printf("I exacta                 = %.10f\n", I_exacta);
    printf("Trapecio compuesto       = %.10f   error = %.10e\n",
           I_trapecio, fabs(I_trapecio - I_exacta));
    printf("Gauss-Legendre 2 puntos  = %.10f   error = %.10e\n",
           I_gauss2, fabs(I_gauss2 - I_exacta));
    printf("Gauss-Legendre 6 puntos  = %.10f   error = %.10e\n",
           I_gauss6, fabs(I_gauss6 - I_exacta));

    return 0;
}
