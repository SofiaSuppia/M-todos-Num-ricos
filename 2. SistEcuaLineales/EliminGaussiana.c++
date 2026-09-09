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

    // Se abre el archivo
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

    // El sistema debe ser cuadrado: n filas = n columnas
    *n = rows;
    printf("Tamanio del sistema: %d x %d\n", *n, *n);

    // Se verifica que no supere el tamanio maximo permitido
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
        // Se leen los coeficientes de la fila i
        for (int j = 0; j < *n; j++) {
            if (fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("Error al leer elemento a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        // Se lee el termino independiente b[i]
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
// Programa principal: Eliminacion Gaussiana con Pivoteo
//
// Resuelve un sistema de ecuaciones lineales Ax = b
// usando el metodo de eliminacion gaussiana con pivoteo parcial.
//
// Pasos del algoritmo:
//   1. Leer datos (matriz A y vector b desde archivo)
//   2. Pivoteo (intercambiar filas si el pivote es muy chico)
//   3. Triangulacion (hacer ceros debajo de la diagonal)
//   4. Calcular determinante (producto de la diagonal)
//   5. Retrosustitucion (resolver de abajo hacia arriba)
//   6. Imprimir resultados
// ============================================================
int main(int argc, char const *argv[]) {

    // --- PASO 1: LECTURA DE DATOS ---
    // Se lee la matriz aumentada (A|b) desde el archivo
    int n;  // tamanio del sistema (n ecuaciones, n incognitas)
    double a[MAX_SIZE][MAX_SIZE];  // matriz de coeficientes
    double b[MAX_SIZE];            // vector de terminos independientes

    if (!read_array_file("matriz.dat", a, b, &n)) {
        return 1;
    }

    // Se imprime el sistema original para verificar que se leyo correctamente
    printf("Sistema original de ecuaciones:\n");
    printf("================================\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.6lf ", a[i][j]);
        }
        printf("| %10.6lf\n", b[i]);
    }
    printf("\n");

    // --- PASO 2: PIVOTEO ---
    // Recorre cada columna (pivote). Si el elemento de la diagonal a[i][i]
    // es muy chico (casi cero), se busca en la columna i la fila con el
    // elemento de mayor valor absoluto y se intercambian filas.
    // Esto evita division por cero y mejora la precision numerica.
    int p;
    double aux;

    for (int i = 0; i < n - 1; i++) {
        p = i;

        // Si el pivote es casi cero, buscar mejor fila
        if (fabs(a[i][i]) < 1e-5) {
            // Se busca la fila con el mayor valor absoluto en la columna i
            for (int l = i + 1; l < n; l++) {
                if (fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l;  // p guarda la fila del mejor pivote
                }
            }

            // Se intercambian las filas i y p en la matriz
            for (int m = i; m < n; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }

            // Se intercambian las filas i y p en el vector b
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        // --- PASO 3: TRIANGULACION (Eliminacion Gaussiana) ---
        // Para cada fila j debajo del pivote i, se calcula:
        //   factor = a[j][i] / a[i][i]
        // Luego se resta 'factor' veces la fila pivote de la fila j,
        // tanto en la matriz A como en el vector b.
        // Al terminar, la matriz queda en forma triangular superior.
        for (int j = i + 1; j < n; j++) {
            double factor = a[j][i] / a[i][i];

            for (int k = i; k < n; k++) {
                a[j][k] = a[j][k] - factor * a[i][k];
            }
            b[j] = b[j] - factor * b[i];
        }
    }

    // Se imprime la matriz y el vector b despues de la triangulacion
    printf("Matriz despues de la eliminacion gaussiana:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("Vector b despues de la eliminacion gaussiana:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // --- PASO 4: CALCULO DEL DETERMINANTE ---
    // El determinante de una matriz triangular es el producto
    // de los elementos de la diagonal principal.
    // Si det = 0, la matriz es singular y no tiene solucion unica.
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det = det * a[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", det);

    if (det == 0) {
        printf("El determinante es cero, el sistema no tiene solucion unica.\n");
        return 1;
    }

    // --- PASO 5: RETROSUSTITUCION ---
    // Se resuelve el sistema triangular superior Ax = b
    // desde la ultima ecuacion hasta la primera:
    //   x[n-1] = b[n-1] / a[n-1][n-1]
    //   x[i] = (b[i] - sum(a[i][j]*x[j], j=i+1..n-1)) / a[i][i]
    double x[MAX_SIZE];
    double suma;

    x[n - 1] = b[n - 1] / a[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--) {
        suma = b[i];
        for (int j = i + 1; j < n; j++) {
            suma = suma - a[i][j] * x[j];
        }
        x[i] = suma / a[i][i];
    }

    // --- PASO 6: IMPRESION DE RESULTADOS ---
    printf("------------------SOLUCION------------------\n");
    printf("La solucion del sistema es:\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %lf\n", i, x[i]);
    }

    return 0;
}
