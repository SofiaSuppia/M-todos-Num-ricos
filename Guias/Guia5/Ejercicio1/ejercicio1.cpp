#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50 // Número máximo de puntos de datos
#define MAX_TAMANO 100 // Tamaño máximo para la matriz de coeficientes

#include "gauss.h" // Incluye funciones para Eliminación Gaussiana

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
 * Función para imprimir el polinomio interpolante Pn(x)
 * @param a Arreglo de coeficientes (a0, a1, a2, ...)
 * @param n Grado del polinomio + 1 (n-1 es la potencia más alta)
 */
void imprimir_polinomio(double a[], int n);

/**
 * Función para calcular el error absoluto entre valores reales y estimados
 * @param fx El valor real de la función evaluada en X̂
 * @param Pn El valor estimado del polinomio
 * @return El error absoluto
 */
double calcular_error(double fx, double Pn);

/**
 * Función para evaluar la función real f(x)
 * @param x El punto en el que se debe evaluar la función (en nuestro caso, X̂)
 * @return El valor de la función en x
 */
double evaluar_funcion(double x);

/* Este ejercicio debe hacerse con cualquier polinomio de Lagrange de cualquier grado.
Pero la implementación se hizo con grado 2 debido a los cálculos extensos */

int main(int argc, char const *argv[]) {
    double X_gorro, suma, producto, error, fx;
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para el cálculo de coeficientes del polinomio
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion[MAX_TAMANO+1];
    // Solución del polinomio interpolante
    double a[MAX_TAMANO+1];
    int n, opcion;

    // Leer puntos de datos desde el archivo
    if (!leer_puntos_datos("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia5\\Ejercicio1\\data.txt", X, Y, &n)) {
        printf("Error al leer los datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    printf("Elija una opción:\n");
    printf("1. Interpolación de Lagrange\n");
    printf("2. Polinomio Interpolante (usando eliminación Gaussiana)\n");
    printf("0. Salir\n");
    scanf("%d", &opcion);
    
    switch(opcion) {
        case 1:
            /* Si X_gorro está dentro del rango de los puntos de datos, es una interpolación
            Si X_gorro está fuera del rango, es una extrapolación */
            // Interpolación de Lagrange
            printf("Introduzca su X̂ (X_gorro) a ser interpolado: ");
            scanf("%lf", &X_gorro);

            // Calcular Pn(X̂)
            suma = 0.0;
            for(int k = 0; k < n; k++) {        
                producto = 1.0;
                for(int i = 0; i < n; i++) {    
                    if(i != k) {
                        // Cnk(X̂) = (X̂ - Xi) / (Xk - Xi)
                        producto = producto * ((X_gorro - X[i]) / (X[k] - X[i])); 
                    }
                }
                // Imprimir coeficiente Cnk(X̂)
                // Podemos ver que suma Cnk = 1
                printf("C%d%d(%.3f) = %.6f\n", n-1, k, X_gorro, producto);
                
                // Pn(X̂) = Σ Yk * Cnk(X̂)
                suma = suma + (Y[k] * producto);
            }
            // error = |f(X̂) - Pn(X̂)|
            fx = 0.0; // Aquí puede definir la función real f(X̂) si se conoce
            error = 0.0; // Si fx es conocida, calcular error

            printf("\nEl valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);
            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", calcular_error(evaluar_funcion(X_gorro), suma));
            break;
        case 2:
            // Polinomio Interpolante
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    // Matriz de Vandermonde: A[i][j] = X[i]^j
                    A[i][j] = pow(X[i], j);
                }
                // Vector de términos independientes: b[i] = Y[i]
                b[i] = Y[i];
            }
        
            // Usamos la función de gauss.h para resolver el sistema con Eliminación Gaussiana
            // Nota: Se asume que 'gauss_elimination' es el nombre de la función en 'gauss.h'
            eliminacion_gaussiana(n, A, b, solucion); 
        
            // Copiamos la solución a a[i] para darle relevancia a nuestro contexto (coeficientes a0, a1, ...)
            for(int i = 0; i < n; i++) {
                a[i] = solucion[i];
            }
        
            printf("------------------SOLUCIÓN------------------\n");
            printf("La solución del sistema (coeficientes del polinomio) es:\n");
            for(int i = 0; i < n; i++) {
                printf("a[%d] = %lf\n", i, a[i]);
            }
            
            printf("\n------------------POLINOMIO INTERPOLANTE------------------\n");
            imprimir_polinomio(a, n);

            // Opcional: Evaluar el polinomio en X_gorro para interpolar
            printf("Introduzca su X̂ (X_gorro) a ser interpolado: ");
            scanf("%lf", &X_gorro);
            suma = 0.0;
            for(int i = 0; i < n; i++) {
                suma = suma + a[i] * pow(X_gorro, i);
            }

            printf("El valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);
            break;
        default: 
            printf("Saliendo del programa.\n");
            return 0;
    }

    return 0;
}

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
        printf("Error: Número de puntos inválido (%d)\n", *n);
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
    printf("    i   |       Xi      |       Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

void imprimir_polinomio(double a[], int n) {
    printf("Pn(x) = ");
    
    // Manejar el primer término (término constante a0)
    if (fabs(a[0]) > 1e-10) {  // Evitar imprimir valores muy pequeños como cero
        printf("%.6f", a[0]);
    } else {
        printf("0");
    }
    
    // Manejar el resto de los términos
    for (int i = 1; i < n; i++) {
        if (fabs(a[i]) > 1e-10) {  // Solo imprimir si el coeficiente es significativo
            // Imprimir signo
            if (a[i] > 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
            
            // Imprimir coeficiente (valor absoluto ya que el signo ya fue impreso)
            double coeficiente = fabs(a[i]);
            if (coeficiente != 1.0) {
                printf("%.6f", coeficiente);
            }
            
            // Imprimir parte variable
            if (i == 1) {
                printf("x");
            } else {
                printf("x^%d", i);
            }
        }
    }
    printf("\n\n");
}

double calcular_error(double fx, double Pn) {
    return fabs(fx - Pn);
}

double evaluar_funcion(double x) {
    // f(x) = x + (2 / x)
    return x + (2 / x);
}