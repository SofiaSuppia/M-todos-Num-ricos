#include <cstdio>
#include <cmath>
#include <iostream>

#define MAX_SIZE 50

using namespace std;

// ============================================================
// Funcion para leer la matriz aumentada desde un archivo.
// Formato esperado del archivo (una fila por ecuacion):
//   a11 a12 ... a1n b1
//   a21 a22 ... a2n b2
//   ...
//   an1 an2 ... ann bn
//
// @param filename  nombre del archivo a leer
// @param a         matriz de coeficientes (n x n)
// @param b         vector de terminos independientes
// @param n         puntero donde se guarda el tamanio del sistema
// @return          true si se leyo correctamente, false si hubo error
// ============================================================
bool read_array_file(const char* filename, double a[][MAX_SIZE], double b[], int* n) {
    FILE *fp;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se puede abrir el archivo '%s'\n", filename);
        return false;
    }

    printf("Archivo '%s' abierto correctamente\n\n", filename);

    // Se cuentan las filas para saber el tamanio del sistema
    int rows = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        rows++;
    }

    *n = rows;
    printf("Tamanio del sistema: %d x %d\n", *n, *n);

    if (*n > MAX_SIZE) {
        printf("Error: sistema demasiado grande (%d). Maximo permitido: %d\n", *n, MAX_SIZE);
        fclose(fp);
        return false;
    }

    // Se cierra y reabre el archivo para reiniciar el puntero de lectura
    fclose(fp);
    fp = fopen(filename, "r");

    // Se lee la matriz aumentada: n coeficientes + 1 termino independiente por fila
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("Error al leer elemento a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        if (fscanf(fp, "%lf", &b[i]) != 1) {
            printf("Error al leer termino independiente b[%d]\n", i);
            fclose(fp);
            return false;
        }
    }

    fclose(fp);
    printf("Matriz leida correctamente\n\n");
    return true;
}

// ============================================================
// Verifica si la matriz es diagonalmente dominante.
// Para cada fila i: |a[i][i]| > suma(|a[i][j]|, j!=i)
// Tambien verifica que no haya ceros en la diagonal.
// ============================================================
bool checkDiagonalDominance(double a[][MAX_SIZE], int n) {
    for (int i = 0; i < n; i++) {
        // Si hay cero en la diagonal, no se puede dividir
        if (fabs(a[i][i]) == 0.0) {
            printf("Error: cero en la diagonal en a[%d][%d]\n", i, i);
            return false;
        }

        // Suma de los valores absolutos de los elementos fuera de la diagonal
        double suma = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                suma += fabs(a[i][j]);
            }
        }

        // Si no es diagonalmente dominante, solo avisa (no cancela)
        if (fabs(a[i][i]) < suma) {
            printf("Aviso: la matriz NO es diagonalmente dominante en la fila %d\n", i);
        }
    }
    return true;
}

// ============================================================
// Inicializa un vector con todos sus elementos en cero
// ============================================================
void initializeVector(double v[], int n) {
    for (int i = 0; i < n; i++) {
        v[i] = 0.0;
    }
}

// ============================================================
// Muestra un vector completo en una sola linea (para la tabla de
// iteraciones).
// ============================================================
void mostrarVector(double xn[], int n) {
    for (int j = 0; j < n; j++) {
        printf("%12.8lf ", xn[j]);
    }
}

// ============================================================
// Calcula el error entre dos vectores (norma euclidea)
// error = sqrt( sum( (xn[i] - xv[i])^2 ) )
// ============================================================
double computeError(double xn[], double xv[], int n) {
    double error = 0.0;
    for (int i = 0; i < n; i++) {
        error += pow(xn[i] - xv[i], 2);
    }
    return sqrt(error);
}

// ============================================================
// Cabecera de la tabla de iteraciones: iter | x[0] x[1] ... | error
// ============================================================
void imprimirCabecera(int n) {
    printf("\n  iter | ");
    for (int j = 0; j < n; j++) {
        printf("    x[%d]      ", j);
    }
    printf("|      error\n");
    printf("  -----------------------------------------------------\n");
}

// ============================================================
// Verifica el resultado calculando el residuo  ||A*x - b|| .
// Si el residuo es muy chico, la solucion es correcta.
// ============================================================
void verificarResiduo(double a[][MAX_SIZE], double b[], double xn[], int n) {
    double residuo = 0.0;
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            s += a[i][j] * xn[j];
        }
        residuo += pow(s - b[i], 2);
    }
    residuo = sqrt(residuo);
    printf("Residuo ||A*x - b|| = %.12lf\n", residuo);
}

// ============================================================
// Metodo de Jacobi
//
// En Jacobi se usan SOLO los valores de la iteracion anterior (xv)
// para calcular todos los valores nuevos (xn). Despues se copian.
//
// Formula: xn[i] = (b[i] - sum(a[i][j]*xv[j], j!=i)) / a[i][i]
//
// Diferencia con Gauss-Seidel: Jacobi usa siempre los valores
// VIEJOS, Gauss-Seidel usa los valores NUEVOS que ya se calcularon.
// ============================================================
void jacobiMethod(double a[][MAX_SIZE], double b[], double xv[], double xn[], int n) {
    double tolerance;
    int iterations = 0;
    double old_error, new_error;

    initializeVector(xv, n);

    printf("\nIngrese tolerancia: ");
    scanf("%lf", &tolerance);

    old_error = 1000.0;

    imprimirCabecera(n);

    do {
        iterations++;

        // Se calculan TODOS los xn usando SOLO los valores viejos xv
        for (int i = 0; i < n; i++) {
            double suma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    suma += a[i][j] * xv[j];
                }
            }
            xn[i] = (b[i] - suma) / a[i][i];
        }

        // Se calcula el error entre las dos iteraciones
        new_error = computeError(xn, xv, n);

        // Verificacion de convergencia: si el error crece, no converge
        if (new_error > old_error) {
            printf("\nEl metodo NO converge (el error crece).\n");
            printf("Iteraciones: %d, Error: %.12lf\n", iterations, new_error);
            return;
        }

        old_error = new_error;

        // RECIEN ACA se copian los valores nuevos a xv para la proxima iteracion
        for (int i = 0; i < n; i++) {
            xv[i] = xn[i];
        }

        printf("  %3d  | ", iterations);
        mostrarVector(xn, n);
        printf("| %.12lf\n", new_error);

    } while (new_error > tolerance && iterations < 10000);

    // Se muestra la solucion
    printf("\n--- SOLUCION (Jacobi) ---\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.12lf\n", i, xn[i]);
    }
    printf("Error: %.12lf\n", new_error);
    printf("Iteraciones: %d\n", iterations);
    verificarResiduo(a, b, xn, n);
}

// ============================================================
// Metodo de Gauss-Seidel
//
// A diferencia de Jacobi, cuando se calcula xn[i] se usan:
//   - Los valores NUEVOS (xn) para j < i (ya calculados en esta iteracion)
//   - Los valores VIEJOS (xv) para j > i (aun no calculados)
//
// Esto hace que converja mas rapido que Jacobi.
// ============================================================
void gaussSeidelMethod(double a[][MAX_SIZE], double b[], double xv[], double xn[], int n) {
    double tolerance;
    int iterations = 0;
    double old_error, new_error;

    initializeVector(xv, n);

    printf("\nIngrese tolerancia: ");
    scanf("%lf", &tolerance);

    old_error = 1000.0;

    imprimirCabecera(n);

    do {
        iterations++;

        for (int i = 0; i < n; i++) {
            double suma = 0.0;

            // Valores NUEVOS (ya calculados en esta iteracion): j < i
            for (int j = 0; j < i; j++) {
                suma += a[i][j] * xn[j];
            }

            // Valores VIEJOS (de la iteracion anterior): j > i
            for (int j = i + 1; j < n; j++) {
                suma += a[i][j] * xv[j];
            }

            xn[i] = (b[i] - suma) / a[i][i];
        }

        // Se calcula el error
        new_error = computeError(xn, xv, n);

        // Verificacion de convergencia
        if (new_error > old_error) {
            printf("\nEl metodo NO converge (el error crece).\n");
            printf("Iteraciones: %d, Error: %.12lf\n", iterations, new_error);
            return;
        }

        old_error = new_error;

        // Se copian los valores nuevos para la proxima iteracion
        for (int i = 0; i < n; i++) {
            xv[i] = xn[i];
        }

        printf("  %3d  | ", iterations);
        mostrarVector(xn, n);
        printf("| %.12lf\n", new_error);

    } while (new_error > tolerance && iterations < 10000);

    // Se muestra la solucion
    printf("\n--- SOLUCION (Gauss-Seidel) ---\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.12lf\n", i, xn[i]);
    }
    printf("Error: %.12lf\n", new_error);
    printf("Iteraciones: %d\n", iterations);
    verificarResiduo(a, b, xn, n);
}

// ============================================================
// Programa principal
// ============================================================
int main(int argc, char const *argv[]) {

    // --- PASO 1: LECTURA DE DATOS ---
    int n;
    double a[MAX_SIZE][MAX_SIZE];
    double b[MAX_SIZE];

    if (!read_array_file("matriz.dat", a, b, &n)) {
        return 1;
    }

    // Se imprime el sistema original
    printf("Sistema original de ecuaciones:\n");
    printf("================================\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.6lf ", a[i][j]);
        }
        printf("| %10.6lf\n", b[i]);
    }
    printf("\n");

    // --- PASO 2: VERIFICACION DE DIAGONAL DOMINANTE ---
    if (!checkDiagonalDominance(a, n)) {
        printf("No se puede continuar. La matriz tiene ceros en la diagonal.\n");
        return 1;
    }
    printf("Verificacion completada.\n\n");

    // --- PASO 3: ELEGIR METODO ---
    int opcion;
    double xv[MAX_SIZE], xn[MAX_SIZE];

    do {
        printf("\n--- MENU ---");
        printf("\n1- Jacobi");
        printf("\n2- Gauss-Seidel");
        printf("\n3- Salir");
        printf("\n------------");
        printf("\nOPCION: ");
        scanf("%d", &opcion);

        switch (opcion) {
            case 1:
                jacobiMethod(a, b, xv, xn, n);
                break;
            case 2:
                gaussSeidelMethod(a, b, xv, xn, n);
                break;
            case 3:
                printf("Saliendo...\n");
                break;
            default:
                printf("Opcion invalida.\n");
        }
    } while (opcion != 3);

    return 0;
}

//4 -1 1 7
//4 -8 1 -21
//-2 1 5 15