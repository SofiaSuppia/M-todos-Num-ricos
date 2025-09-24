#include <cstdio>
#include <stdlib.h>
#include <cmath>

// Ahora puedes manejar matrices de hasta 50x50
#define MAX_SIZE 50

/**
 * Lee una matriz aumentada desde un archivo de texto llamado data.dat
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

/**
 * Calcula la norma de Frobenius de una matriz.
 * @param a Matriz de coeficientes
 * @param n Tamaño de la matriz
 * @return La norma de Frobenius
 */
double norma_frobenius(double a[][MAX_SIZE+1], int n);

/* La matriz a resolver es:

En esta matriz usaremos c1 = 20514

20514c1 4424c2 978c3 224c4 = 20514
4424c1 978c2 224c3 54c4 = 4424
978c1 224c2 54c3 14c4 = 978
224c1 54c2 14c3 4c4 = 224


En otros términos:
A = | 20514 4424 978 224 |
    | 4424 978 224 54    |
    | 978 224 54 14      |
    | 224 54 14 4        |

b = | 20514 |
    | 4424  |
    | 978   |
    | 224   |
    

Aplicaremos eliminación de Gauss con pivoteo parcial para resolver el sistema de ecuaciones.

La solución es:
    c1 = 1.0
    c2 = 0.0
    c3 = 0.0
    c4 = 0.0

Determinante de A = 144.0

||A|| = 21518.53

----------------------------------------------------
Ahora usaremos la misma matriz pero con c1 = 20515

20515c1 4424c2 978c3 224c4 = 20541
4424c1 978c2 224c3 54c4 = 4424
978c1 224c2 54c3 14c4 = 978
224c1 54c2 14c3 4c4 = 224


En otros términos:
A = | 20515 4424 978 224 |
    | 4424 978 224 54    |
    | 978 224 54 14      |
    | 224 54 14 4        |

b = | 20541 |
    | 4424  |
    | 978   |
    | 224   |

La solución es:
    c1 = 0.6428
    c2 = 3.75
    c3 = -12.39
    c4 = 12.75

Determinante de A = 224.0

||A|| = 21519.48


Nótese que la solución de entrada difiere significativamente en ambos casos, donde solo se considera una pequeña perturbación adicional en el nuevo sistema de ecuaciones.
Estos problemas se encuentran frecuentemente en la solución de sistemas de ecuaciones, donde la matriz de coeficientes se llama mal condicionada.
Esto ocurre cuando la matriz de coeficientes es "casi singular", es decir, cuando su determinante es pequeño en comparación con la norma euclidiana de la matriz de coeficientes.
En lenguaje matemático:
Det(A) <<<< ||A|| <-- Esto significa que la matriz está mal condicionada.

*/

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma, aux;
    
    // Definiendo arreglos usando el tamaño máximo global
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Leer el arreglo desde el archivo usando la función
    if(!leer_archivo_arreglo("C:\\Users\\Juan Carlos\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia3\\Ejercicio6\\data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Sistema de ecuaciones original:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);      // 10 espacios totales, 6 decimales
    }
    printf("\n");

    // Calculamos la norma de Frobenius de la matriz
    double normaA = norma_frobenius(a, n);
    printf("------------------NORMA------------------\n");
    printf("La norma de Frobenius de la matriz es: %lf\n\n", normaA);


    // Recorrer las filas de la matriz (eliminación de Gauss)
    for(int i = 1; i <= n-1; i++) {

        // ¿Qué sucede si a[i][i] es cercano a cero?
        // Usaremos pivoteo parcial para evitar la división por cero o inestabilidad numérica
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
                a[i][m] = aux; // Intercambiar filas p e i
            }
            aux = b[p];
            b[p] = b[i];
            b[i] = aux; // Intercambiar el término independiente
        }

        // Hace cero los elementos debajo de la diagonal en la columna actual
        for(int j = i+1; j <= n; j++) {
            factor = a[j][i] / a[i][i]; // Sin el signo negativo
            
            // Recorre las columnas en la fila j
            for(int k = i; k <= n; k++) {
                a[j][k] = a[j][k] - factor * a[i][k]; // Corregido: a[j][k] - factor * a[i][k]
            }
            b[j] = b[j] - factor * b[i]; // Corregido: b[j] - factor * b[i]
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
        // Limpiar los valores muy pequeños que deberían ser cero
        if(fabs(X[i]) < 1e-10) {
            X[i] = 0.0;
        }
        printf("X[%d] = %.6lf\n", i, X[i]);
    }
    

    return 0;
}

bool leer_archivo_arreglo(const char* nombre_archivo, double a[][MAX_SIZE+1], double b[], int* n) {
    FILE *fp;
    char c;
    
    // Abrir el archivo de datos
    fp = fopen(nombre_archivo, "r");
    if (fp == NULL) {
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }
    
    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar las filas en el archivo
    /* int filas = 0;
    while((c = fgetc(fp)) != EOF) {
        if(c == '\n') {
            filas++;
        }
    } */

    // Esta alternativa funciona mejor que la anterior
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

    // Cerrar y volver a abrir el archivo para reiniciar el puntero
    fclose(fp);
    fp = fopen(nombre_archivo, "r");
    
    // Verificar el tamaño máximo usando el tamaño global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_SIZE);
        fclose(fp);
        return false;
    }

    // Leer la matriz aumentada desde el archivo
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
    printf("✅ Arreglo leído correctamente desde el archivo\n\n");
    
    return true;
}

double norma_frobenius(double a[][MAX_SIZE+1], int n) {
    double suma = 0.0;
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            suma += a[i][j] * a[i][j];
        }
    }
    return sqrt(suma);
}