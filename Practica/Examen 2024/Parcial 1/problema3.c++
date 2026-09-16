#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

#define N 3

// ============================================================
// Muestra el sistema como matriz aumentada [ A | b ]
// ============================================================
void mostrarSistema(double A[][N], double b[]) {
    for (int i = 0; i < N; i++) {
        printf("|");
        for (int j = 0; j < N; j++) {
            printf("%10.6lf ", A[i][j]);
        }
        printf("| %10.6lf\n", b[i]);
    }
    printf("\n");
}

// ============================================================
// Norma euclidea de la diferencia de dos vectores (u - v).
// Es la misma definicion de error que usa el codigo
// JacobiYGaussSeidel.c++ del tema 2.
// ============================================================
double normaEuclidea(double u[], double v[]) {
    double suma = 0.0;
    for (int i = 0; i < N; i++) {
        suma += pow(u[i] - v[i], 2);
    }
    return sqrt(suma);
}

int main() {

    // ============================================================
    // SISTEMA DEL PARCIAL 1 - 2024 (Problema 3)
    //   2*x2 + 5*x3  = -14        -> la ecuacion 1 no tiene x1
    //   x1 - 4*x2 + 2*x3 = -19
    //   5*x1 + x2 - 2*x3 = 16
    // ============================================================
    double A[N][N] = {
        { 0,  2,  5},
        { 1, -4,  2},
        { 5,  1, -2}
    };
    double b[N] = { -14, -19, 16 };

    // Solucion exacta (x1, x2, x3) conocida = (1, 3, -4)
    double solExacta[N] = { 1, 3, -4 };

    printf("=========================================================\n");
    printf(" PROBLEMA 3 - SISTEMA DE ECUACIONES LINEALES\n");
    printf("   2*x2 + 5*x3      = -14\n");
    printf("   x1   - 4*x2 + 2*x3 = -19\n");
    printf("   5*x1 +   x2  - 2*x3 = 16\n");
    printf("=========================================================\n\n");
    printf("Solucion exacta conocida: x = (1, 3, -4)\n\n");

    // ============================================================
    // PARTE a) ELIMINACION GAUSSIANA (con pivoteo parcial)
    // ============================================================
    printf("=========== PARTE a) ELIMINACION GAUSSIANA ===========\n\n");

    double Aa[N][N], ba[N];
    for (int i = 0; i < N; i++) {
        ba[i] = b[i];
        for (int j = 0; j < N; j++) Aa[i][j] = A[i][j];
    }

    printf("Matriz original:\n");
    mostrarSistema(Aa, ba);

    // --- Verificacion de si hace falta pivoteo ---
    printf("El elemento a[1][1] = %lf.\n", Aa[0][0]);
    printf("Como a[1][1] = 0, la eliminacion simple dividiria por cero\n");
    printf("-> SI es necesario emplear pivoteo parcial.\n\n");

    // --- Pivoteo parcial: buscar el mayor valor absoluto en la columna 0 ---
    int filaPiv = 0;
    for (int i = 1; i < N; i++) {
        if (fabs(Aa[i][0]) > fabs(Aa[filaPiv][0])) filaPiv = i;
    }
    if (filaPiv != 0) {
        printf("PIVOTEO: se intercambian la fila 1 y la fila %d "
               "(pivote = %.6lf).\n\n", filaPiv + 1, Aa[filaPiv][0]);
        for (int j = 0; j < N; j++) {
            double aux = Aa[0][j]; Aa[0][j] = Aa[filaPiv][j]; Aa[filaPiv][j] = aux;
        }
        double aux = ba[0]; ba[0] = ba[filaPiv]; ba[filaPiv] = aux;
        mostrarSistema(Aa, ba);
    }

    // --- Eliminacion hacia adelante ---
    for (int j = 1; j < N; j++) {
        double m = Aa[j][0] / Aa[0][0];
        printf("m[%d][1] = %.6lf / %.6lf = %.6lf\n", j + 1, Aa[j][0], Aa[0][0], m);
        for (int k = 0; k < N; k++) Aa[j][k] -= m * Aa[0][k];
        ba[j] -= m * ba[0];
    }
    printf("\nLuego de la eliminacion en la columna 1:\n");
    mostrarSistema(Aa, ba);

    double m32 = Aa[2][1] / Aa[1][1];
    printf("m[3][2] = %.6lf / %.6lf = %.6lf\n", Aa[2][1], Aa[1][1], m32);
    for (int k = 1; k < N; k++) Aa[2][k] -= m32 * Aa[1][k];
    ba[2] -= m32 * ba[1];

    printf("\nLuego de la eliminacion en la columna 2:\n");
    mostrarSistema(Aa, ba);

    // --- Retrosustitucion ---
    double x[N];
    x[2] = ba[2] / Aa[2][2];
    x[1] = (ba[1] - Aa[1][2] * x[2]) / Aa[1][1];
    x[0] = (ba[0] - Aa[0][1] * x[1] - Aa[0][2] * x[2]) / Aa[0][0];

    printf("RETROSUSTITUCION:\n");
    printf("  x3 = %.6lf / %.6lf = %.6lf\n", ba[2], Aa[2][2], x[2]);
    printf("  x2 = (%.6lf - %.6lf*%.6lf) / %.6lf = %.6lf\n",
           ba[1], Aa[1][2], x[2], Aa[1][1], x[1]);
    printf("  x1 = (%.6lf - %.6lf*%.6lf - %.6lf*%.6lf) / %.6lf = %.6lf\n\n",
           ba[0], Aa[0][1], x[1], Aa[0][2], x[2], Aa[0][0], x[0]);

    printf("SOLUCION EXACTA: x1 = %.6lf   x2 = %.6lf   x3 = %.6lf\n",
           x[0], x[1], x[2]);
    printf("Coincide con x = (1, 3, -4)\n\n");

    // ============================================================
    // PARTE b) GAUSS-SEIDEL (posibilidad + primera iteracion a mano)
    // ============================================================
    printf("=========== PARTE b) GAUSS-SEIDEL ===========\n\n");

    printf("El metodo NO puede usarse directamente con el orden del\n");
    printf("enunciado: la ecuacion 1 no tiene x1 -> a[1][1] = 0, y ademas\n");
    printf("esa fila no es diagonalmente dominante.\n\n");

    printf("Con intercambio de filas (E3, E2, E1) el sistema queda:\n");
    double Ag[N][N] = {
        { 5,  1, -2},
        { 1, -4,  2},
        { 0,  2,  5}
    };
    double bg[N] = { 16, -19, -14 };
    mostrarSistema(Ag, bg);

    printf("Verificacion de dominancia diagonal estricta (condicion\n");
    printf("suficiente de convergencia de Gauss-Seidel):\n");
    for (int i = 0; i < N; i++) {
        double suma = 0.0;
        for (int j = 0; j < N; j++) if (j != i) suma += fabs(Ag[i][j]);
        printf("  fila %d: |%.1lf|  vs  %.1lf  ->  %sdominante\n",
               i + 1, Ag[i][i], suma,
               (fabs(Ag[i][i]) > suma) ? "SI es " : "NO es ");
    }
    printf("\nComo el sistema reordenado es estrictamente diagonalmente\n");
    printf("dominante, Gauss-Seidel converge. SI es posible resolverlo\n");
    printf("mediante el intercambio de filas.\n\n");

    printf("Despejes (x0 = (0, 0, 0)):\n");
    printf("  x1 = (16 - x2 + 2*x3)/5\n");
    printf("  x2 = (19 + x1 + 2*x3)/4\n");
    printf("  x3 = (-14 - 2*x2)/5\n\n");

    double x0 = 0.0, x1g = 0.0, x2g = 0.0, x3g = 0.0;
    x1g = (16 - x2g + 2 * x3g) / 5;
    printf("  ITERACION 1:\n");
    printf("    x1 = (16 - %.1lf + 2*%.1lf)/5 = %.6lf\n", x2g, x3g, x1g);
    x2g = (19 + x1g + 2 * x3g) / 4;
    printf("    x2 = (19 + %.6lf + 2*%.1lf)/4 = %.6lf\n", x1g, x3g, x2g);
    x3g = (-14 - 2 * x2g) / 5;
    printf("    x3 = (-14 - 2*%.6lf)/5 = %.6lf\n\n", x2g, x3g);

    double x1it[N] = { x1g, x2g, x3g };
    printf("  x^(1) = (%.6lf, %.6lf, %.6lf)\n", x1g, x2g, x3g);
    printf("  Error absoluto exacto ||x^(1) - x_exacta|| = %.6lf\n\n",
           normaEuclidea(x1it, solExacta));

    // ============================================================
    // PARTE c) GAUSS-SEIDEL POR CODIGO hasta |error| <= 10^-4
    // ============================================================
    printf("=========== PARTE c) GAUSS-SEIDEL POR CODIGO ===========\n\n");

    printf("Ingrese tolerancia (error absoluto): ");
    double tol;
    scanf("%lf", &tol);

    double xv[N] = { 0, 0, 0 };   // x0 = (0, 0, 0)
    double xn[N];
    double errorEst;
    int iter = 0;

    printf("\n  iter   error estimado    error exacto\n");
    printf("  ----------------------------------------\n");
    do {
        iter++;
        for (int i = 0; i < N; i++) {
            double suma = 0.0;
            for (int j = 0; j < i; j++) suma += Ag[i][j] * xn[j];
            for (int j = i + 1; j < N; j++) suma += Ag[i][j] * xv[j];
            xn[i] = (bg[i] - suma) / Ag[i][i];
        }
        errorEst = normaEuclidea(xn, xv);
        printf("  %3d    %.10lf    %.10lf\n", iter, errorEst,
               normaEuclidea(xn, solExacta));
        for (int i = 0; i < N; i++) xv[i] = xn[i];
    } while (errorEst > tol && iter < 10000);

    printf("  ----------------------------------------\n\n");
    printf("SOLUCION OBTENIDA:\n");
    printf("  x1 = %.10lf\n  x2 = %.10lf\n  x3 = %.10lf\n", xn[0], xn[1], xn[2]);
    printf("  Error absoluto estimado = %.10lf\n", errorEst);
    printf("  Error absoluto exacto   = %.10lf\n", normaEuclidea(xn, solExacta));
    printf("  Numero de iteraciones   = %d\n", iter);

    return 0;
}