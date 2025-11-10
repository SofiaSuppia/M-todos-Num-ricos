#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 20 
#define MAX_TAMANO 100 
#include "gauss.h" 

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
 * @param a Arreglo de coeficientes
 * @param n Grado del polinomio + 1 (n-1 es la potencia más alta)
 */
void imprimir_polinomio(double a[], int n);

/**
 * Función para calcular el error absoluto entre valores reales y estimados
 * @param f_real El valor real de la función evaluado en X̂
 * @param Pn El valor estimado del polinomio
 * @return El error absoluto
 */
double calcular_error(double f_real, double Pn);

/**
 * Función para evaluar la función real f(x)
 * @param x El punto en el que se debe evaluar la función (en nuestro caso, X̂)
 * @return El valor de la función en x
 */
double evaluar_funcion(double x);

/**
 * Función para calcular el polinomio de Lagrange en un punto dado
 * @param x El punto en el que se debe evaluar el polinomio
 * @param X Arreglo de valores X
 * @param Y Arreglo de valores Y
 * @param n Número de puntos de datos
 * @return El valor del polinomio de Lagrange en x
 */
double polinomio_lagrange(double x, double X[], double Y[], int n);

/**
 * Función para calcular los coeficientes del polinomio de Lagrange (forma expandida)
 * @param X Arreglo de valores X
 * @param Y Arreglo de valores Y
 * @param n Número de puntos de datos
 * @param coefs Arreglo para almacenar los coeficientes calculados
 */
void coeficientes_lagrange(double X[], double Y[], int n, double coefs[]);

/**
 * Función para imprimir la forma simbólica expandida del polinomio de Lagrange
 * @param coefs Arreglo de coeficientes
 * @param n Número de puntos de datos
 */
void imprimir_lagrange_expandido(double coefs[], int n);

int main(int argc, char const *argv[]) {
    double X_gorro, suma, producto, error_calc, f_real, Pn;
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para el cálculo de coeficientes del polinomio (Matriz de Vandermonde)
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion[MAX_TAMANO+1];
    // Solución del polinomio interpolante (coeficientes)
    double a[MAX_TAMANO+1];
    int n, opcion;

    // Leer puntos de datos desde el archivo
    if (!leer_puntos_datos("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Suppia_Sofia_25-09-25\\Ejercicio4\\data.dat", X, Y, &n)) {
        printf("Error al leer los datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    printf("Elija una opcion:\n");
    printf("1. Interpolacion de Lagrange (calculo directo)\n");
    printf("2. Polinomio Interpolante (usando Eliminacion Gaussiana)\n");
    printf("3. Polinomio de Lagrange Expandido (muestra coeficientes)\n");
    printf("0. Salir\n");
    scanf("%d", &opcion);
    
    switch(opcion) {
        case 1:
            /* Si X_gorro está dentro del rango de los puntos de datos, es una interpolación.
               Si X_gorro está fuera del rango, es una extrapolación. */
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
                // Podemos ver que la suma de Cnk = 1
                printf("C%d%d(%.3f) = %.6f\n", n-1, k, X_gorro, producto);
                
                // Pn(X̂) = Σ Yk * Cnk(X̂)
                suma = suma + (Y[k] * producto);
            }

            // error = |f(X̂) - Pn(X̂)|
            f_real = 0.0; // Aquí se puede definir la función real f(X̂) si se conoce
            error_calc = 0.0; // Si f_real es conocida, calcular error

            f_real = evaluar_funcion(X_gorro);
            error_calc = calcular_error(f_real, suma);

            printf("\nEl valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);
            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error_calc);
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
            eliminacion_gaussiana(n, A, b, solucion);
        
            // Copiamos la solución a a[i] para darle relevancia a nuestro contexto (coeficientes a0, a1, ...)
            for(int i = 0; i < n; i++) {
                a[i] = solucion[i];
            }
        
            printf("------------------SOLUCION------------------\n");
            printf("La solucion del sistema es:\n");
            for(int i = 0; i < n; i++) {
                printf("a[%d] = %lf\n", i, a[i]);
            }
            
            printf("\n------------------POLINOMIO INTERPOLANTE------------------\n");
            // Pn(x) = a0 + a1*x + a2*x^2 + ... + a(n-1)*x^(n-1)
            imprimir_polinomio(a, n);

            // f(x) está definida en la función evaluar_funcion
            f_real = 0.0;

            // error = |f(X̂) - Pn(X̂)|
            error_calc = 0.0; // Si f_real es conocida, calcular error

            printf("Introduzca su X̂ a ser interpolado: ");
            scanf("%lf", &X_gorro);
            suma = 0.0;
            for(int i = 0; i < n; i++) {
                suma = suma + a[i] * pow(X_gorro, i);
            }
            f_real = evaluar_funcion(X_gorro);
            error_calc = calcular_error(f_real, suma);

            printf("El valor interpolado en X̂ = %lf es: %lf\n", X_gorro, suma);

            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error_calc);
            break;
            
        case 3:
            /* Si X_gorro está dentro del rango de los puntos de datos, es una interpolación.
               Si X_gorro está fuera del rango, es una extrapolación. */
            // Interpolación de Lagrange
            printf("Introduzca su X̂ a ser interpolado: ");
            scanf("%lf", &X_gorro);

            Pn = polinomio_lagrange(X_gorro, X, Y, n);

            // Calcular coeficientes del polinomio expandido
            double coefs[MAX_PUNTOS];
            coeficientes_lagrange(X, Y, n, coefs);
            
            printf("------------------POLINOMIO DE LAGRANGE------------------\n");
            imprimir_lagrange_expandido(coefs, n);

            // error = |f(x) - Pn(x)| 
            // Nota: Para datos experimentales no hay función teórica exacta
            f_real = evaluar_funcion(X_gorro);
            error_calc = calcular_error(f_real, Pn);

            // El mensaje de salida se ha ajustado para ser genérico (no solo temperatura)
            printf("\nEl valor interpolado en x_gorro = %.2f es: %.6f \n", X_gorro, Pn);
            if (f_real != 0.0) {
                printf("\nerror = |f(x) - Pn(x)| = %lf\n", error_calc);
            } else {
                printf("\nNota: El error no se calcula - se usan datos experimentales sin funcion teorica exacta\n");
            }
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
        printf("Error: No se pudo leer el numero de puntos de datos\n");
        fclose(p_archivo);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("Error: Numero de puntos invalido (%d). Maximo permitido: %d\n", *n, MAX_PUNTOS);
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

double polinomio_lagrange(double x, double X[], double Y[], int n) {
    double suma, producto;

    // Calcular Pn(X̂)
    suma = 0.0;
    for(int k = 0; k < n; k++) {        
        producto = 1.0;
        for(int i = 0; i < n; i++) {    
            if(i != k) {
                // Cnk(X̂) = (X̂ - Xi) / (Xk - Xi)
                producto = producto * ((x - X[i]) / (X[k] - X[i])); 
            }
        }
        // Imprimir coeficiente Cnk(X̂)
        // Podemos ver que la suma de Cnk = 1
        printf("C%d%d(%.3f) = %.6f\n", n-1, k, x, producto);
        
        // Pn(X̂) = Σ Yk * Cnk(X̂)
        suma = suma + (Y[k] * producto);
    }

    return suma;
}



void coeficientes_lagrange(double X[], double Y[], int n, double coefs[]) {
    // Inicializar coeficientes a cero
    for (int i = 0; i < n; i++)
        coefs[i] = 0.0;

    for (int k = 0; k < n; k++) {
        double base[MAX_PUNTOS] = {1.0};  
        int grado = 0;

        double denominador = 1.0;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;

            // Multiplica (x - X[i]) al polinomio base
            for (int j = grado; j >= 0; j--)
                base[j+1] += base[j];
            grado++;

            for (int j = grado; j >= 0; j--)
                base[j] = (j>0 ? base[j-1] : 0) - X[i]*(j>=0?base[j]:0);

            denominador *= (X[k] - X[i]);
        }

        double factor = Y[k] / denominador;
        for (int j = 0; j <= grado; j++)
            coefs[j] += base[j] * factor;
    }
}

void imprimir_lagrange_expandido(double coefs[], int n) {
    printf("P_%d(x) = ", n-1);
    int es_el_primero = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(coefs[i]) < 1e-12) continue; // Ignorar coeficientes muy cercanos a cero
        
        if (!es_el_primero) {
            if (coefs[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coefs[i] < 0) {
            printf("-");
        }
        
        printf("%.6f", fabs(coefs[i]));
        if (i > 0) printf("*x^%d", i);
        es_el_primero = 0;
    }
    printf("\n");
}

double calcular_error(double f_real, double Pn) {
    return fabs(f_real - Pn);
}

// Nota: Para datos experimentales, no hay función teórica exacta
// Retornamos 0 para indicar que no se puede calcular error teórico
double evaluar_funcion(double x) {
    // No hay función teórica conocida para estos datos experimentales
    // Se podría definir una si fuera conocida
    return 0.0; // Indica que no hay función teórica
}