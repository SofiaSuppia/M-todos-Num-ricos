#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 20
#define MAX_TAMANO 100 // Usado para la matriz A en la Eliminación Gaussiana
#include "gauss.h" // Se asume que este archivo incluye la función gauss_elimination

/**
 * Función para leer pares de datos Xi, Yi de un archivo
 * Formato esperado: la primera línea contiene el número de puntos, luego los pares xi yi
 * @param nombre_archivo Nombre del archivo a leer
 * @param X Arreglo para almacenar los valores de X
 * @param Y Arreglo para almacenar los valores de Y
 * @param n Puntero para almacenar el número de puntos de datos leídos
 * @return 1 si tiene éxito, 0 en caso contrario
 */
int leer_puntos_datos(const char* nombre_archivo, double X[], double Y[], int* n);

/**
 * Función para mostrar los puntos de datos
 * @param X Arreglo de valores de X
 * @param Y Arreglo de valores de Y
 * @param n Número de puntos de datos
 */
void imprimir_puntos_datos(double X[], double Y[], int n);

/**
 * Función para imprimir el polinomio interpolador Pn(x)
 * @param a Arreglo de coeficientes
 * @param n Grado del polinomio (n-1 es la potencia más alta)
 */
void imprimir_polinomio(double a[], int n);

/**
 * Función para calcular el error absoluto entre valores reales y estimados
 * @param f_x El valor real de la función evaluado en X_gorro
 * @param P_n El valor del polinomio estimado
 * @return El error absoluto
 */
double calcular_error(double f_x, double P_n);

/**
 * Función para evaluar la función real f(x)
 * @param x El punto en el que se debe evaluar la función (en nuestro caso, X_gorro)
 * @return El valor de la función en x
 */
double evaluar_funcion(double x);

/**
 * Función para calcular el polinomio de Lagrange en un punto dado
 * @param x El punto en el que se debe evaluar el polinomio
 * @param X Arreglo de valores de X
 * @param Y Arreglo de valores de Y
 * @param n Número de puntos de datos
 * @return El valor del polinomio de Lagrange en x
 */
double polinomio_lagrange(double x, double X[], double Y[], int n);

/**
 * Función para calcular los coeficientes del polinomio de Lagrange expandido (forma canónica)
 * @param X Arreglo de valores de X
 * @param Y Arreglo de valores de Y
 * @param n Número de puntos de datos
 * @param coefs Arreglo para almacenar los coeficientes calculados
 */
void coeficientes_lagrange(double X[], double Y[], int n, double coefs[]);

/**
 * Función para imprimir la forma simbólica del polinomio de Lagrange expandido
 * @param coefs Arreglo de coeficientes
 * @param n Número de puntos de datos
 */
void imprimir_lagrange_expandido(double coefs[], int n);

/*
    Polinomio de Lagrange (ejemplo de resultado):
    -15.33 + 12.6089*x + 28.0355*x² - 13.056*x³
*/

int main(int argc, char const *argv[]) {
    double X_gorro, suma, producto, error_actual, f_x, P_n;
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para el cálculo de coeficientes polinomiales (forma canónica)
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion_sistema[MAX_TAMANO+1];
    // Solución (coeficientes) del polinomio interpolador
    double a[MAX_TAMANO+1];
    int n, opcion;

    // Leer puntos de datos del archivo
    if (!leer_puntos_datos("exercise3.txt", X, Y, &n)) {
        printf("Fallo al leer los datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    printf("Elija una opción:\n");
    printf("1. Interpolación de Lagrange (evaluación)\n");
    printf("2. Polinomio Interpolador (por Eliminación Gaussiana)\n");
    printf("3. Polinomio de Lagrange (Forma Expandida e Interpolación)\n");
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
                // Se puede ver que la suma de los Cnk es 1
                printf("C%d%d(%.3f) = %.6f\n", n-1, k, X_gorro, producto);
                
                // Pn(X̂) = Σ Yk * Cnk(X̂)
                suma = suma + (Y[k] * producto);
            }

            // f_x es la función real. Si no se conoce, se usa f_x = 0.0
            f_x = evaluar_funcion(X_gorro);
            error_actual = calcular_error(f_x, suma);

            printf("\nEl valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);
            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error_actual);
            break;

        case 2:
            // Polinomio Interpolador (a0 + a1*x + ... + a_n-1 * x^n-1)
            // Llenar la matriz A (matriz de Vandermonde) y el vector b
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    A[i][j] = pow(X[i], j); // x^0, x^1, x^2, ...
                }
                b[i] = Y[i];
            }
        
            // Usamos la función de gauss.h para resolver el sistema con Eliminación Gaussiana
            // La solución se almacena en el vector 'solucion_sistema'
            // Nota: Se asume que 'gauss_elimination' usa indexación 0 si se usa MAX_SIZE sin +1
            //       o que la matriz A y b son ajustadas para 1-based indexing dentro de la función.
            //       Aquí se usa 0-based indexing para A y b.
            // Para mantener la consistencia con el código original, se mantienen los límites n.
            // (Si gauss_elimination usa 1-based indexing, esto podría causar problemas de índice)
            gauss_elimination(n, A, b, solucion_sistema);
        
            // Copiamos la solución a 'a[i]' para darle relevancia a nuestro contexto (coeficientes a_i)
            for(int i = 0; i < n; i++) {
                a[i] = solucion_sistema[i];
            }
        
            printf("------------------SOLUCIÓN DEL SISTEMA------------------\n");
            printf("La solución del sistema es:\n");
            for(int i = 0; i < n; i++) {
                printf("a[%d] = %lf\n", i, a[i]);
            }
            
            printf("\n------------------POLINOMIO INTERPOLADOR------------------\n");
            // Pn(x) = a0 + a1*x + a2*x^2 + ... + a(n-1)*x^(n-1)
            imprimir_polinomio(a, n);

            printf("Introduzca su X̂ (X_gorro) a ser interpolado: ");
            scanf("%lf", &X_gorro);
            
            // Evaluar Pn(X̂)
            suma = 0.0;
            for(int i = 0; i < n; i++) {
                suma = suma + a[i] * pow(X_gorro, i);
            }
            
            f_x = evaluar_funcion(X_gorro);
            error_actual = calcular_error(f_x, suma);

            printf("El valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);
            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error_actual);
            break;

        case 3:
            /* Si X_gorro está dentro del rango de los puntos de datos, es una interpolación
               Si X_gorro está fuera del rango, es una extrapolación */
            // Interpolación de Lagrange
            printf("Introduzca su X̂ (X_gorro) a ser interpolado: ");
            scanf("%lf", &X_gorro);

            P_n = polinomio_lagrange(X_gorro, X, Y, n);

            // Calcular coeficientes del polinomio expandido
            double coefs[MAX_PUNTOS];
            coeficientes_lagrange(X, Y, n, coefs);
            
            printf("------------------POLINOMIO DE LAGRANGE EXPANDIDO------------------\n");
            imprimir_lagrange_expandido(coefs, n);

            // error = |f(T) - Pn(T)| 
            // Nota: Para datos experimentales no hay función teórica exacta
            f_x = evaluar_funcion(X_gorro);
            error_actual = calcular_error(f_x, P_n);

            printf("\nEl valor interpolado en X̂ = %.2f es: %.6f \n", X_gorro, P_n);
            if (f_x != 0.0) {
                printf("\nerror = |f(x) - Pn(x)| = %lf\n", error_actual);
            } else {
                printf("\nNota: El error no se calcula - se utilizan datos experimentales sin una función teórica exacta\n");
            }
            break;

        default: 
            printf("Saliendo del programa.\n");
            return 0;
    }

    return 0;
}

// --------------------------------------------------------------------------------------------------
// --- IMPLEMENTACIÓN DE FUNCIONES ---
// --------------------------------------------------------------------------------------------------

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
    
    // Leer puntos de datos (pares xi yi)
    for (int i = 0; i < *n; i++) {
        if (fscanf(p_archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el punto de datos %d\n", i + 1);
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
    printf("    i  |       Xi      |       Yi      \n");
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
            // Imprimir el signo
            if (a[i] > 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
            
            // Imprimir el coeficiente (valor absoluto ya que el signo ya se imprimió)
            double coef = fabs(a[i]);
            if (coef != 1.0) {
                printf("%.6f", coef);
            }
            
            // Imprimir la parte variable
            if (i == 1) {
                printf("x");
            } else {
                printf("x^%d", i);
            }
        }
    }
    printf("\n\n");
}

double polinomio_lagrange(double x, double X[], double Y[], int n) {
    double suma, producto;

    // Calcular Pn(x)
    suma = 0.0;
    for(int k = 0; k < n; k++) {          
        producto = 1.0;
        for(int i = 0; i < n; i++) {    
            if(i != k) {
                // Lk(x) = Π (x - Xi) / (Xk - Xi)
                producto = producto * ((x - X[i]) / (X[k] - X[i])); 
            }
        }
        // Imprimir el coeficiente Lk(x)
        printf("L%d%d(%.3f) = %.6f\n", n-1, k, x, producto);
        
        // Pn(x) = Σ Yk * Lk(x)
        suma = suma + (Y[k] * producto);
    }

    return suma;
}

double calcular_error(double f_x, double P_n) {
    return fabs(f_x - P_n);
}

double evaluar_funcion(double x) {
    // Definición de la función real f(x) (ejemplo: x + 2/x)
    return x + (2/x);
}


void coeficientes_lagrange(double X[], double Y[], int n, double coefs[]) {
    // Inicializar coeficientes a cero
    for (int i = 0; i < n; i++)
        coefs[i] = 0.0;

    // Calcular el polinomio de base Lk(x) y sumarlo a la forma canónica
    for (int k = 0; k < n; k++) {
        double base_polinomio[MAX_PUNTOS] = {1.0}; // Coeficientes del producto Π(x-Xi)
        int grado = 0;

        double denominador = 1.0;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;

            // Multiplica (x - X[i]) al polinomio base
            for (int j = grado; j >= 0; j--)
                base_polinomio[j+1] += base_polinomio[j];
            grado++;

            for (int j = grado; j >= 0; j--)
                base_polinomio[j] = (j>0 ? base_polinomio[j-1] : 0) - X[i]*(j>=0?base_polinomio[j]:0);

            denominador *= (X[k] - X[i]);
        }

        // El factor de ponderación: Yk / Lk(Xk)
        double factor = Y[k] / denominador;
        
        // Sumar los coeficientes ponderados a la forma canónica final
        for (int j = 0; j <= grado; j++)
            coefs[j] += base_polinomio[j] * factor;
    }
}

void imprimir_lagrange_expandido(double coefs[], int n) {
    printf("P_%d(x) = ", n-1);
    int es_primero = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(coefs[i]) < 1e-12) continue; // Si es casi cero, se ignora

        if (!es_primero) {
            // Imprimir signo (+ o -)
            if (coefs[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coefs[i] < 0) {
            // Imprimir signo (-) para el primer término
            printf("-");
        }
        
        // Imprimir coeficiente (valor absoluto)
        printf("%.6f", fabs(coefs[i]));
        
        // Imprimir la parte de la variable
        if (i > 0) {
             if (i == 1) printf("*x");
             else printf("*x^%d", i);
        }
        es_primero = 0;
    }
    printf("\n");
}