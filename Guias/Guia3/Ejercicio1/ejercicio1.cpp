#include <cstdio>
#include <stdlib.h>
#include <cmath>

// Ahora puedes manejar matrices de hasta 50x50
#define MAX_SIZE 50

/**
 * Lee una matriz aumentada de un archivo de texto llamado data.dat
 * El formato esperado es:
 * a11 a12 ... a1n b1
 * a21 a22 ... a2n b2
 * ...
 * an1 an2 ... ann bn
 * @param nombre_archivo Nombre del archivo a leer
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param n Puntero al tamaño del sistema (número de ecuaciones)
 * @return true si el archivo se leyó correctamente, false en caso contrario
 * */
bool leer_archivo_arreglo(const char* nombre_archivo, double a[][MAX_SIZE+1], double b[], int* n);


/* La matriz a resolver es:
4x1 - x2 + 2x3 + 3x4 = 20
0x1 -2x2 + 7x3 - 4x4 = -7
0x1 + 0x2 + 6x3 + 5x4 = 4
0x1 + 0x2 + 0x3 + 3x4 = 6

En otros términos:
A = | 4 -1  2  3 |
    | 0 -2  7 -4 |
    | 0  0  6  5 |
    | 0  0  0  3 |
b = | 20 |
    | -7 |
    |  4 |
    |  6 |

Aplicaremos eliminación de Gauss con pivoteo parcial para resolver el sistema de ecuaciones.

La solución es:
x1 = 3.0
x2 = -4.0
x3 = -1.0
x4 = 2.0

Determinante de A = -144.0

*/

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma, aux;

    // Definiendo arreglos usando la variable global MAX_SIZE
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Leer matriz desde archivo usando la función
    if(!leer_archivo_arreglo("C:\\Users\\Juan Carlos\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia3\\Ejercicio1\\data.dat", a, b, &n)) {
        return 1;
    }

    printf("Sistema de ecuaciones original:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);       // 10 espacios totales, 6 decimales
    }
    printf("\n");
    // Recorrer las filas de la matriz (eliminación de Gauss)
    for(int i = 1; i <= n-1; i++) {

        // ¿Qué sucede si a[i][i] es cercano a cero?
        // Usaremos pivoteo parcial para evitar división por cero o inestabilidad numérica
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l; // Encontramos la fila con el elemento más grande en la columna i
                }
            }
            for(int m = i; m <= n; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux; // Intercambiar las filas p e i
            }
            aux = b[p];
            b[p] = b[i];
            b[i] = aux; // Intercambiar el término independiente
        }

        // Hace cero los elementos debajo de la diagonal en la columna actual
        for(int j = i+1; j <= n; j++) {
            factor = a[j][i] / a[i][i];

            // Recorre las columnas en la fila j
            for(int k = i; k <= n; k++) {
                a[j][k] = a[j][k] - factor * a[i][k];
            }
            b[j] = b[j] - factor * b[i];
        }
    }

    // Imprimimos la matriz después de la eliminación de Gauss, es más conveniente
    printf("La matriz después de la eliminación de Gauss es:\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("El vector b después de la eliminación de Gauss es:\n");
    for(int i = 1; i <= n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");


    // Verificamos el determinante de la matriz
    producto = 1.0;
    for(int i = 1; i <= n; i++) {
        producto = producto * a[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", producto);

    if(producto == 0) {
        printf("El determinante de la matriz es cero, el sistema no tiene una solución única.\n");
        exit(0); // Salir con código de error
    }

    // Realizamos la sustitución hacia atrás para encontrar la solución
    X[n] = b[n] / a[n][n];

    for(int i = n-1; i >= 1; i--) {
        suma = b[i];
        for(int j = i+1; j <= n; j++) {
            suma = suma - a[i][j] * X[j];
        }
        suma = suma / a[i][i];
        X[i] = suma;
    }
    printf("------------------SOLUCION------------------\n");
    printf("La solución del sistema es:\n");
    for(int i = 1; i <= n; i++) {
        printf("X[%d] = %lf\n", i, X[i]);
    }


    return 0;
}

bool leer_archivo_arreglo(const char* nombre_archivo, double a[][MAX_SIZE+1], double b[], int* n) {
    FILE *fp;
    char c;

    // Abrir archivo de datos
    fp = fopen(nombre_archivo, "r");
    if (fp == NULL) {
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }

    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar filas en el archivo
    int filas = 0;
    while(!feof(fp)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), fp) != NULL) {
            filas++;
        }
    }

    // El sistema debe ser cuadrado, por lo tanto n = filas
    *n = filas;
    printf("📊 Tamaño del sistema: %d x %d\n", *n, *n);

    // Cerrar y reabrir el archivo para reiniciar el puntero
    fclose(fp);
    fp = fopen(nombre_archivo, "r");

    // Verificar el tamaño máximo usando la variable global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_SIZE);
        fclose(fp);
        return false;
    }

    // Leer la matriz aumentada del archivo
    // Formato esperado: cada fila contiene n coeficientes + 1 término independiente
    // Ejemplo para 3x3: a11 a12 a13 b1
    //                   a21 a22 a23 b2
    //                   a31 a32 a33 b3

    int i, j;
    for(i = 1; i <= *n; i++) {
        // Leyendo los coeficientes de la matriz
        for(j = 1; j <= *n; j++) {
            if(fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("❌ Error al leer el elemento a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        // Leer el término independiente
        if(fscanf(fp, "%lf", &b[i]) != 1) {
            printf("❌ Error al leer el término independiente b[%d]\n", i);
            fclose(fp);
            return false;
        }
    }

    fclose(fp);
    printf("✅ Arreglo leído correctamente del archivo\n\n");

    return true;
}