// ============================================================================
//  ENUNCIADO  (Problema n.° 3 — Examen 2025)
//
//  Los datos de la tabla pueden ser ajustados por una funcion de la forma
//
//        f(x) = a·sen(x) + b·cos(x)
//
//  El metodo de minimos cuadrados se basa en la minimizacion del error
//  cuadratico medio o, equivalentemente, en la minimizacion de la funcion
//
//        ε = Σ(k=1..N) ( yk - f(xk) )²
//
//  con N el numero de puntos a ajustar. Esto permite obtener, mediante
//  derivacion respecto a los parametros a y b, los coeficientes a y b que
//  mejor ajustan los datos.
//
//  a) Obtener, para la funcion dada y los valores tabulados, el sistema de
//     ecuaciones lineales que deben satisfacer los coeficientes a y b para
//     que ε sea minimo. Escribir la matriz de coeficientes y el vector
//     independiente correspondientes.
//
//  b) Resolver el sistema planteado, empleando eliminacion gaussiana, y
//     escribir el vector solucion, asi como los valores de St, Sr y el
//     coeficiente de correlacion.
//
//  Datos:
//     x      y
//    1.0   -0.236
//    1.2    0.209
//    1.5    0.853
//    2.0    1.746
//    2.6    2.231
//    2.9    2.163
//    3.0    2.110
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50
#define MAX_TAMANO 100
#include "gauss.h"   // eliminacion_gaussiana(n, A, b, solucion)

// ============================================================
// Lee los pares (xi, yi) desde un archivo con formato:
//   primera linea: numero de puntos n
//   luego n lineas: xi yi
// (igual que Guia5/regresionMinimosCuadrados.cpp)
// ============================================================
int leer_puntos_datos(const char* nombre_archivo, double X[], double Y[], int* n) {
    FILE *p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return 0;
    }
    if (fscanf(p_archivo, "%d", n) != 1) {
        printf("Error: No se pudo leer el numero de puntos de datos\n");
        fclose(p_archivo);
        return 0;
    }
    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("Error: Numero de puntos invalido (%d). Maximo: %d\n", *n, MAX_PUNTOS);
        fclose(p_archivo);
        return 0;
    }
    for (int i = 0; i < *n; i++) {
        if (fscanf(p_archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el punto %d\n", i + 1);
            fclose(p_archivo);
            return 0;
        }
    }
    fclose(p_archivo);
    printf("Archivo '%s' abierto correctamente. Se leyeron %d puntos.\n\n", nombre_archivo, *n);
    return 1;
}

// ============================================================
// Muestra los puntos de datos en una tabla
// ============================================================
void imprimir_puntos_datos(double X[], double Y[], int n) {
    printf("Puntos de datos:\n");
    printf("================\n");
    printf("     i  |        Xi      |       Yi\n");
    printf("--------+---------------+--------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d   | %12.6f   | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

int main() {
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    double A[MAX_TAMANO + 1][MAX_TAMANO + 1];
    double b[MAX_TAMANO + 1], solucion[MAX_TAMANO + 1];
    int n;

    printf("===========================================================\n");
    printf(" PROBLEMA 3 - AJUSTE POR MINIMOS CUADRADOS (SENO/COSENO)\n");
    printf(" f(x) = a·sen(x) + b·cos(x)\n");
    printf("===========================================================\n\n");

    // --- Lectura de datos ---
    if (!leer_puntos_datos("data.txt", X, Y, &n)) {
        printf("Error al leer datos. Saliendo.\n");
        return 1;
    }
    imprimir_puntos_datos(X, Y, n);

    // ============================================================
    // PARTE a) SISTEMA DE ECUACIONES NORMALES (2x2)
    //
    // Se minimiza  ε = Σ ( yk - a·sen(xk) - b·cos(xk) )²
    //
    //   ∂ε/∂a = 0  ->  a·Σsen²  + b·Σsen·cos = Σ y·sen
    //   ∂ε/∂b = 0  ->  a·Σsen·cos + b·Σcos²  = Σ y·cos
    //
    // Matriz de coeficientes y vector independiente:
    //   A = [ Σsen²    Σsen·cos ]       b = [ Σ y·sen ]
    //       [ Σsen·cos Σcos²    ]           [ Σ y·cos ]
    // ============================================================
    double s_sen2 = 0.0, s_cos2 = 0.0, s_sencos = 0.0;
    double s_ysen = 0.0, s_ycos = 0.0;

    printf("Tabla de calculo de sumatorias (sen, cos y terminos cruzados):\n");
    printf("    i |     Xi     |     Yi     |   sen(xi)  |   cos(xi)  "
           "|  sen²(xi)  |  sen·cos  |   cos²(xi)  |   y·sen    |   y·cos\n");
    printf("------+-----------+-----------+------------+------------+"
           "-----------+-----------+------------+------------+------------\n");

    for (int i = 0; i < n; i++) {
        double s = sin(X[i]);
        double c = cos(X[i]);
        s_sen2   += s * s;
        s_cos2   += c * c;
        s_sencos += s * c;
        s_ysen   += Y[i] * s;
        s_ycos   += Y[i] * c;

        printf("%4d  | %9.4lf | %9.4lf | %10.6lf | %10.6lf | %11.7lf "
               "| %9.6lf | %11.7lf | %10.6lf | %10.6lf\n",
               i + 1, X[i], Y[i], s, c, s * s, s * c, c * c, Y[i] * s, Y[i] * c);
    }

    printf("\nSumatorias:\n");
    printf("  Σsen²  = %.10lf\n", s_sen2);
    printf("  Σcos²  = %.10lf\n", s_cos2);
    printf("  Σsen·cos = %.10lf\n", s_sencos);
    printf("  Σy·sen = %.10lf\n", s_ysen);
    printf("  Σy·cos = %.10lf\n", s_ycos);

    printf("\nSISTEMA DE ECUACIONES NORMALES (parte a):\n");
    printf("  | %.6lf   %.6lf | | a |   | %.6lf |\n", s_sen2, s_sencos, s_ysen);
    printf("  | %.6lf   %.6lf | | b | = | %.6lf |\n", s_sencos, s_cos2, s_ycos);

    A[0][0] = s_sen2;   A[0][1] = s_sencos;
    A[1][0] = s_sencos; A[1][1] = s_cos2;
    b[0] = s_ysen;      b[1] = s_ycos;

    // ============================================================
    // PARTE b) RESOLUCION CON ELIMINACION GAUSSIANA (gauss.h)
    // ============================================================
    printf("\n=========== PARTE b) ELIMINACION GAUSSIANA ===========\n\n");
    eliminacion_gaussiana(2, A, b, solucion);

    double a = solucion[0];   // coeficiente de sen(x)
    double bcoef = solucion[1]; // coeficiente de cos(x)

    printf("------------------SOLUCION------------------\n");
    printf("El vector solucion es:\n");
    printf("  a = %.10lf\n", a);
    printf("  b = %.10lf\n", bcoef);
    printf("\n  f(x) = %.10lf·sen(x) + %.10lf·cos(x)\n\n", a, bcoef);

    // ============================================================
    // St, Sr y coeficiente de correlacion r
    // ============================================================
    double suma_y = 0.0, promedio_y, St = 0.0, Sr = 0.0;

    for (int i = 0; i < n; i++) suma_y += Y[i];
    promedio_y = suma_y / n;

    for (int i = 0; i < n; i++) {
        St += pow(Y[i] - promedio_y, 2);
    }

    printf("Predicciones, residuos y verificacion de Sr:\n");
    printf("    i |     Xi     |     Yi     |      f(xi)     |   ei = yi-f(xi) |   ei²\n");
    printf("------+-----------+-----------+----------------+------------------+----------\n");
    for (int i = 0; i < n; i++) {
        double fxi = a * sin(X[i]) + bcoef * cos(X[i]);
        Sr += pow(Y[i] - fxi, 2);
        printf("%4d  | %9.4lf | %9.4lf | %14.10lf | %16.10lf | %.12lf\n",
               i + 1, X[i], Y[i], fxi, Y[i] - fxi, (Y[i] - fxi) * (Y[i] - fxi));
    }

    double r = sqrt((St - Sr) / St);

    printf("\n------------------RESULTADOS------------------\n");
    printf("  Promedio de y:  yp = %.10lf\n", promedio_y);
    printf("  St (suma total de cuadrados) : %.10lf\n", St);
    printf("  Sr (suma de residuos al cuadrado): %.10lf\n", Sr);
    printf("  Coeficiente de correlacion r = sqrt((St - Sr)/St) = %.10lf\n", r);

    double distancia = fabs(r - 1.0);
    if (distancia < 0.05)      printf("  El ajuste es EXCELENTE (r muy cercano a 1).\n");
    else if (distancia < 0.15) printf("  El ajuste es muy bueno.\n");
    else if (distancia < 0.3)  printf("  El ajuste es aceptable.\n");
    else                       printf("  El ajuste es pobre o inadecuado.\n");

    return 0;
}