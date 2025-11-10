#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50 // Máximo número de puntos de datos
#define MAX_TAMANO 100 // Máximo tamaño para la matriz A (4*(n-1) x 4*(n-1))
#include "gauss.h" // Incluye la cabecera para la eliminación gaussiana

/**
 * Función para leer pares de datos Xi, Yi desde un archivo
 * Formato esperado: la primera línea contiene el número de puntos, luego pares xi yi
 * @param nombre_archivo Nombre del archivo a leer
 * @param X Arreglo para almacenar los valores de X
 * @param Y Arreglo para almacenar los valores de Y
 * @param n Puntero para almacenar el número de puntos de datos leídos
 * @return 1 si tiene éxito, 0 en caso contrario
 */
int leer_puntos_datos(const char* nombre_archivo, double X[], double Y[], int* n);

/**
 * Función para mostrar los puntos de datos
 * @param X Arreglo de valores X
 * @param Y Arreglo de valores Y
 * @param n Número de puntos de datos
 */
void imprimir_puntos_datos(double X[], double Y[], int n);

/**
 * Función para imprimir los coeficientes de los splines cúbicos
 * @param X Arreglo de valores X (puntos de datos)
 * @param solucion Arreglo que contiene los coeficientes de la solución (a, b, c, d)
 * @param n Número de puntos de datos
 */
void imprimir_splines_cubicos(double X[], double solucion[], int n);

/**
 * Función para evaluar el spline cúbico en un punto x dado
 * @param X Arreglo de valores X (puntos de datos)
 * @param solucion Arreglo de coeficientes del spline
 * @param n Número de puntos de datos
 * @param x El valor x a evaluar
 * @return El valor y interpolado
 */
double evaluar_spline(double X[], double solucion[], int n, double x);

/**
 * Función para calcular e imprimir los splines lineales
 * @param X Arreglo de valores X (puntos de datos)
 * @param Y Arreglo de valores Y (puntos de datos)
 * @param n Número de puntos de datos
 */
void spline_lineal(double X[], double Y[], int n);

int main(int argc, char const *argv[]) {
    // Arreglos para los puntos de datos a leer del archivo de texto
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para la matriz de coeficientes y el vector de términos independientes
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion[MAX_TAMANO+1];
    // n = número de puntos de datos
    int n, opcion, fila1, fila2;
    
    // Leer puntos de datos desde el archivo
    if (!leer_puntos_datos("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Suppia_Sofia_25-09-25\\Ejercicio4\\data.dat", X, Y, &n)) {
        printf("Error al leer datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    // Inicializar A y b a cero
    for (int i = 0; i < 4*(n-1); i++) {
        b[i] = 0.0;
        for (int j = 0; j < 4*(n-1); j++) {
            A[i][j] = 0.0;
        }
    }

    printf("Elija una opción para la interpolación con splines:\n");
    printf("1. Spline Lineal\n");
    printf("2. Spline Cúbico (Natural)\n");
    scanf("%d", &opcion);

    switch(opcion) {
        case 1:
            spline_lineal(X, Y, n);
            // Evaluamos puntos x_gorro para verificar si el spline funciona correctamente
            printf("¿Desea evaluar un punto x_gorro? (1 para sí, 0 para no): ");
            scanf("%d", &opcion);
            while(opcion) {
                double x_gorro;
                printf("Ingrese el valor de x_gorro: ");
                scanf("%lf", &x_gorro);
                // Encontrar el intervalo correcto para x_gorro
                if(x_gorro < X[0] || x_gorro > X[n-1]) {
                    printf("x_gorro está fuera de los límites [%lf, %lf]. Ingrese un valor dentro del rango.\n", X[0], X[n-1]);
                } else {
                    // Localizar el intervalo [X[k], X[k+1]]
                    int k = 0;
                    while(k < n-1 && x_gorro > X[k+1]) {
                        k++;
                    }
                    // Fórmula de interpolación lineal: y = Yk + m * (x - Xk)
                    double y_gorro = Y[k] + (Y[k+1] - Y[k]) * (x_gorro - X[k]) / (X[k+1] - X[k]);
                    printf("El valor interpolado en x_gorro = %.3f es y_gorro = %.3f\n", x_gorro, y_gorro);
                }
                printf("¿Desea evaluar otro punto x_gorro? (1 para sí, 0 para no): ");
                scanf("%d", &opcion);
            }
        break;
        case 2:
            // Para n puntos, hay n-1 splines cúbicos, cada uno con 4 coeficientes (a,b,c,d)
            // El sistema es de 4(n-1) x 4(n-1)

            // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)]

            // 1. Las primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos puntos finales
            for(int k = 0; k < n-1; k++) {
                // Spline k pasa por el punto (X[k], Y[k]) -> a*X[k]^3 + b*X[k]^2 + c*X[k] + d = Y[k]
                for(int j = 0; j <= 3; j++) {
                    A[2*k][4*k+j] = pow(X[k], 3-j); // a, b, c, d (potencias 3, 2, 1, 0)
                }
                b[2*k] = Y[k];

                // Spline k pasa por el punto (X[k+1], Y[k+1]) -> a*X[k+1]^3 + b*X[k+1]^2 + c*X[k+1] + d = Y[k+1]
                for(int j = 0; j <= 3; j++) {
                    A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                }
                b[2*k+1] = Y[k+1];
            }
        
            // 2. Las siguientes n-2 ecuaciones: Continuidad de la primera derivada en los puntos interiores X[k+1]
            for(int k = 0; k < n-2; k++) {
                int fila = 2*(n-1) + k;
                // Primera derivada del spline k en X[k+1]
                for(int j = 0; j <= 2; j++) {
                    // Derivadas de a*x^3, b*x^2, c*x: 3a*x^2, 2b*x, c
                    A[fila][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                }
                // Primera derivada del spline k+1 en X[k+1]
                for(int j = 0; j <= 2; j++) {
                    A[fila][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j); // Restamos el lado derecho (que se iguala a 0)
                }
                b[fila] = 0.0;
            }
        
            // 3. Las siguientes n-2 ecuaciones: Continuidad de la segunda derivada en los puntos interiores X[k+1]
            for(int k = 0; k < n-2; k++) {
                int fila = 2*(n-1) + (n-2) + k;
                // Segunda derivada del spline k en X[k+1]: 6a*x + 2b
                A[fila][4*k] = 6 * X[k+1];
                A[fila][4*k+1] = 2;
                // Segunda derivada del spline k+1 en X[k+1]: -(6a*x + 2b)
                A[fila][4*(k+1)] = -6 * X[k+1];
                A[fila][4*(k+1)+1] = -2;
                b[fila] = 0.0;
            }
        
            // 4. Las dos condiciones de frontera (2 ecuaciones adicionales): Spline Natural
            // La segunda derivada es cero en los puntos finales X[0] y X[n-1]
            fila1 = 4*(n-1) - 2; // Penúltima fila
            fila2 = 4*(n-1) - 1; // Última fila

            // Segunda derivada = 0 en X[0] (primer spline): 6a*X[0] + 2b = 0
            A[fila1][0] = 6 * X[0];
            A[fila1][1] = 2;
            b[fila1] = 0.0;

            // Segunda derivada = 0 en X[n-1] (último spline): 6a*X[n-1] + 2b = 0
            A[fila2][4*(n-2)] = 6 * X[n-1];
            A[fila2][4*(n-2)+1] = 2;
            b[fila2] = 0.0;
        
            // Usamos la función de gauss.h para resolver el sistema con Eliminación Gaussiana
            eliminacion_gaussiana(4*(n-1), A, b, solucion);
        
            // Imprimir los coeficientes del spline cúbico
            imprimir_splines_cubicos(X, solucion, n);
        
            // Evaluamos puntos x_gorro para verificar si el spline funciona correctamente
            printf("¿Desea evaluar un punto x_gorro? (1 para sí, 0 para no): ");
            scanf("%d", &opcion);
            while(opcion) {
                double x_gorro;
                printf("Ingrese el valor de x_gorro: ");
                scanf("%lf", &x_gorro);
                double y_gorro = evaluar_spline(X, solucion, n, x_gorro);
                printf("El valor interpolado en x_gorro = %.3f es y_gorro = %.3f\n", x_gorro, y_gorro);
                printf("¿Desea evaluar otro punto x_gorro? (1 para sí, 0 para no): ");
                scanf("%d", &opcion);
            }
            break;
        default:
            printf("Opción inválida. Intente de nuevo.\n");
            break;
    }

    return 0;
}


// Implementación de funciones auxiliares

int leer_puntos_datos(const char* nombre_archivo, double X[], double Y[], int* n) {
    FILE *p_archivo;
    
    p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return 0;
    }
    
    printf("Archivo '%s' abierto exitosamente\n", nombre_archivo);
    
    // Leer número de puntos de datos
    if (fscanf(p_archivo, "%d", n) != 1) {
        printf("Error: No se pudo leer el número de puntos de datos\n");
        fclose(p_archivo);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("Error: Número de puntos inválido (%d). Máximo permitido: %d\n", *n, MAX_PUNTOS);
        fclose(p_archivo);
        return 0;
    }
    
    // Leer puntos de datos
    for (int i = 0; i < *n; i++) {
        if (fscanf(p_archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el punto de dato %d\n", i + 1);
            fclose(p_archivo);
            return 0;
        }
    }
    
    fclose(p_archivo);
    printf("Se leyeron %d puntos de datos exitosamente\n\n", *n);
    return 1;
}


void imprimir_puntos_datos(double X[], double Y[], int n) {
    printf("Puntos de Datos:\n");
    printf("=============\n");
    printf("     i  |        Xi      |       Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

void imprimir_splines_cubicos(double X[], double solucion[], int n) {
    printf("------------------SOLUCIÓN------------------\n");
    printf("Coeficientes de los Splines Cúbicos:\n");
    for(int k = 0; k < n-1; k++) {
        // Los coeficientes para el spline k están en solucion[4k], solucion[4k+1], etc.
        printf("Spline %d (desde X[%d]=%.3f hasta X[%d]=%.3f):\n", k+1, k, X[k], k+1, X[k+1]);
        printf("  S%d(x) = %.6fx³ + %.6fx² + %.6fx + %.6f\n", 
                k+1, solucion[4*k], solucion[4*k+1], solucion[4*k+2], solucion[4*k+3]);
    }
    printf("\n");
}

double evaluar_spline(double X[], double solucion[], int n, double x) {
    int k;

    // Encontrar el intervalo correcto del spline
    if (x <= X[0]) {
        k = 0;  // Usar el primer spline si x es menor que el primer punto
    } else if (x >= X[n-1]) {
        k = n - 2; // Usar el último spline si x es mayor que el último punto
    } else {
        // Buscar el intervalo [X[k], X[k+1]]
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluar S_k(x) = a*x^3 + b*x^2 + c*x + d
    double y = solucion[4*k] * pow(x, 3) +
                solucion[4*k+1] * pow(x, 2) +
                solucion[4*k+2] * x +
                solucion[4*k+3];

    return y;
}

void spline_lineal(double X[], double Y[], int n) {
    printf("------------------SPLINES LINEALES------------------\n");
    for (int k = 0; k < n - 1; k++) {
        // Calcular la pendiente (m) para el intervalo k
        double mk = (Y[k+1] - Y[k]) / (X[k+1] - X[k]);
        printf("Spline %d (desde X[%d]=%.3f hasta X[%d]=%.3f):\n", k+1, k, X[k], k+1, X[k+1]);
        // Forma punto-pendiente: f(X) = Yk + mk * (X - Xk)
        printf("  f_%d(X) = %.6f + %.6f * (X - %.6f)\n\n", k+1, Y[k], mk, X[k]);
    }
}