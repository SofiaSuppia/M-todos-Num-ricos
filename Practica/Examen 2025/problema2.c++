// ============================================================================
//  ENUNCIADO  (Problema n.° 2 — Examen 2025)
//
//  Dado el siguiente sistema de ecuaciones:
//
//      4·x1 +  x2                     −  x6 =  0
//      2·x1 + 10·x2 −  x3        − 2·x5    =  9
//             3·x2 +  8·x3        −  x5 − 3·x6 =  7
//             − x2         + 4·x4       +  x6 = 20
//       x1   − 2·x2                 + 5·x5     = 22
//                    3·x3  − 2·x4        + 6·x6 = 37
//
//  a)  Obtenga la solucion del sistema empleando eliminacion gaussiana.
//      ¿Es necesario emplear pivoteo? Justifique.
//
//  b)  ¿Es posible resolver el sistema utilizando el metodo de Gauss-Seidel?
//      (Justifique). En caso afirmativo, calcule, mediante el codigo
//      correspondiente, la solucion con un error absoluto de 10^-4 tomando
//      como valor inicial x0j = 1 (j = 1,...,6). Escriba la solucion obtenida
//      (con 10 cifras decimales), error absoluto estimado y el numero de
//      iteraciones empleadas.
// ============================================================================
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

#define N 6

// ============================================================
// Muestra el sistema como matriz aumentada [ A | b ]
// ============================================================
void mostrarSistema(double A[][N], double b[]) {
    for (int i = 0; i < N; i++) {
        printf("|");
        for (int j = 0; j < N; j++) {
            printf("%8.4lf ", A[i][j]);
        }
        printf("| %9.4lf\n", b[i]);
    }
    printf("\n");
}

// ============================================================
// Error absoluto entre dos vectores (norma euclidea) — igual que
// computeError() en JacobiYGaussSeidel.c++ del tema 2:
//   error = sqrt( sum( (xn[i] - xv[i])^2 ) )
// ============================================================
double normaEuclidea(double xn[], double xv[]) {
    double error = 0.0;
    for (int i = 0; i < N; i++) {
        error += (xn[i] - xv[i]) * (xn[i] - xv[i]);
    }
    return sqrt(error);
}

// ============================================================
// Verifica el resultado calculando el residuo ||A*x - b||
// ============================================================
void verificarResiduo(double A[][N], double b[], double x[]) {
    double residuo = 0.0;
    for (int i = 0; i < N; i++) {
        double s = 0.0;
        for (int j = 0; j < N; j++) {
            s += A[i][j] * x[j];
        }
        residuo += (s - b[i]) * (s - b[i]);
    }
    residuo = sqrt(residuo);
    printf("  Residuo ||A*x - b|| = %.12lf\n", residuo);
}

int main() {

    // ============================================================
    // SISTEMA DEL EXAMEN 2025 (Problema 2), orden 6x6
    // ============================================================
    double A[N][N] = {
        { 4,  1,  0,  0,  0, -1},
        { 2, 10, -1,  0, -2,  0},
        { 0,  3,  8,  0, -1, -3},
        { 0, -1,  0,  4,  0,  1},
        { 1, -2,  0,  0,  5,  0},
        { 0,  0,  3, -2,  0,  6}
    };
    double b[N] = { 0, 9, 7, 20, 22, 37 };

    printf("===========================================================\n");
    printf(" PROBLEMA 2 - SISTEMA DE ECUACIONES LINEALES (6x6)\n");
    printf("===========================================================\n");
    printf("  4·x1 +  x2                           -  x6 =  0\n");
    printf("  2·x1 + 10·x2 -  x3             - 2·x5    =  9\n");
    printf("         3·x2 +  8·x3             -  x5 - 3·x6 =  7\n");
    printf("         - x2          + 4·x4        +   x6 = 20\n");
    printf("   x1   - 2·x2                   +  5·x5    = 22\n");
    printf("                3·x3  - 2·x4          +  6·x6 = 37\n\n");

    // ============================================================
    // PARTE a) ELIMINACION GAUSSIANA
    // ============================================================
    printf("=========== PARTE a) ELIMINACION GAUSSIANA ===========\n\n");

    double Aa[N][N], ba[N];
    for (int i = 0; i < N; i++) {
        ba[i] = b[i];
        for (int j = 0; j < N; j++) Aa[i][j] = A[i][j];
    }

    printf("Sistema original:\n");
    mostrarSistema(Aa, ba);

    // --- Justificacion de por que NO hace falta pivoteo ---
    printf("¿Es necesario emplear pivoteo?\n");
    printf("  Los elementos de la diagonal son: 4, 10, 8, 4, 5, 6\n");
    printf("  => todos son DISTINTOS de cero, no se divide por cero.\n");
    printf("  Ademas, la matriz es estrictamente diagonalmente dominante,\n");
    printf("  por lo que en cada columna el mayor valor absoluto ya esta en\n");
    printf("  la diagonal (ni el pivoteo parcial haria un intercambio).\n");
    printf("  => NO es necesario emplear pivoteo.\n\n");

    // --- Triangulacion columnas 1..N-1 (sin pivoteo) ---
    for (int i = 0; i < N - 1; i++) {
        printf("=== PASO: columna %d (pivote a[%d][%d] = %.6lf) ===\n",
               i + 1, i + 1, i + 1, Aa[i][i]);

        for (int j = i + 1; j < N; j++) {
            if (Aa[j][i] == 0.0) continue;
            double m = Aa[j][i] / Aa[i][i];
            printf("  factor m[%d][%d] = %.6lf / %.6lf = %.6lf\n",
                   j + 1, i + 1, Aa[j][i], Aa[i][i], m);
            for (int k = i; k < N; k++) {
                Aa[j][k] -= m * Aa[i][k];
            }
            ba[j] -= m * ba[i];
        }

        printf("  Sistema luego de eliminar la columna %d:\n", i + 1);
        mostrarSistema(Aa, ba);
    }

    printf("MATRIZ TRIANGULAR SUPERIOR FINAL:\n");
    mostrarSistema(Aa, ba);

    // --- Determinante (producto de la diagonal; sin intercambio no
    //     hace falta corregir signo) ---
    double det = 1.0;
    for (int i = 0; i < N; i++) det *= Aa[i][i];
    printf("DETERMINANTE = %.6lf (producto de la diagonal)\n\n", det);

    // --- Retrosustitucion paso a paso ---
    double x[N];
    printf("RETROSUSTITUCION:\n");
    x[N - 1] = ba[N - 1] / Aa[N - 1][N - 1];
    printf("  x6 = %.10lf\n", x[N - 1]);

    for (int i = N - 2; i >= 0; i--) {
        double suma = ba[i];
        for (int j = i + 1; j < N; j++) {
            suma -= Aa[i][j] * x[j];
        }
        x[i] = suma / Aa[i][i];
        printf("  x%d = %.10lf\n", i + 1, x[i]);
    }

    printf("\nSOLUCION POR ELIMINACION GAUSSIANA:\n");
    for (int i = 0; i < N; i++) {
        printf("  x%d = %.10lf\n", i + 1, x[i]);
    }
    verificarResiduo(A, b, x);
    printf("\n");

    // ============================================================
    // PARTE b) GAUSS-SEIDEL
    // ============================================================
    printf("=========== PARTE b) GAUSS-SEIDEL ===========\n\n");

    // --- Justificacion de convergencia: dominancia por fila ---
    printf("¿Es posible resolver el sistema con Gauss-Seidel?\n");
    printf("  Se verifica la dominancia estricta por fila (|a[i][i]| > suma):\n\n");
    for (int i = 0; i < N; i++) {
        double suma = 0.0;
        for (int j = 0; j < N; j++) {
            if (j != i) suma += fabs(A[i][j]);
        }
        printf("  Fila %d: |%g| > %.2lf  =>  %s\n",
               i + 1, A[i][i], suma,
               (fabs(A[i][i]) > suma) ? "domina" : "NO domina");
    }
    printf("\n  La matriz es ESTRICTAMENTE DIAGONALMENTE DOMINANTE en todas");
    printf("\n  las filas => Gauss-Seidel CONVERGE para cualquier valor inicial.");
    printf("\n  => SI es posible.\n\n");

    // --- Gauss-Seidel con x0j = 1, tolerancia de error absoluto 1e-4 ---
    double xn[N], xv[N];
    for (int j = 0; j < N; j++) xv[j] = 1.0;   // valor inicial x0j = 1

    double tolerance = 1e-4;
    double error, old_error = 1.0e9;
    int iter = 0;

    printf("  iter |      x1            x2            x3            x4"
           "            x5            x6    |       error\n");
    printf("  -----+-----------------------------------------------------------------------"
           "-----------+---------------\n");
    printf("    0  |");
    for (int j = 0; j < N; j++) printf(" %13.8lf", xv[j]);
    printf(" |        ----\n");

    do {
        iter++;

        // Gauss-Seidel: al calcular xn[i] se usan los valores NUEVOS ya
        // calculados (j < i) y los VIEJOS de la iteracion anterior (j > i)
        for (int i = 0; i < N; i++) {
            double suma = 0.0;
            for (int j = 0; j < i; j++)     suma += A[i][j] * xn[j];
            for (int j = i + 1; j < N; j++) suma += A[i][j] * xv[j];
            xn[i] = (b[i] - suma) / A[i][i];
        }

        error = normaEuclidea(xn, xv);

        // se copian los nuevos para la proxima iteracion
        for (int j = 0; j < N; j++) xv[j] = xn[j];

        printf("  %3d  |", iter);
        for (int j = 0; j < N; j++) printf(" %13.8lf", xn[j]);
        printf(" | %17.12lf\n", error);

    } while (error > tolerance && iter < 10000);

    printf("\nRESULTADOS DE GAUSS-SEIDEL:\n");
    printf("  Solucion (10 cifras decimales):\n");
    for (int j = 0; j < N; j++) {
        printf("    x%d = %.10lf\n", j + 1, xn[j]);
    }
    printf("  Error absoluto estimado: %.10le\n", error);
    printf("  Numero de iteraciones:   %d\n", iter);
    verificarResiduo(A, b, xn);
    printf("\n");

    return 0;
}