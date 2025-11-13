#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Maximo numero de puntos de datos a leer (para el Spline)
#define MAX_POINTS 20 
// Tamano maximo de la matriz para la Eliminacion Gaussiana
#define MAX_SIZE 100 

// Incluye la funcion de Eliminacion Gaussiana para resolver el sistema Ax=b
#include "gauss.h"

/**
 * @brief Funcion para definir la funcion f(x) o, si hay datos, la evalua usando el Spline.
 * @param x El punto en el que se evalua la funcion.
 * @return El valor de la funcion en x.
 */
double f(double x);

/**
 * @brief Funcion para leer pares de datos (Xi, Yi) de un archivo.
 * Formato esperado: la primera linea contiene el numero de puntos, seguido de pares xi yi.
 * @param filename Nombre del archivo a leer.
 * @param X Array para almacenar los valores de X.
 * @param Y Array para almacenar los valores de Y.
 * @param n Puntero para almacenar el numero de puntos de datos leidos.
 * @return 1 si tiene exito, 0 en caso contrario.
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * @brief Funcion para mostrar los puntos de datos.
 * @param X Array de valores de X.
 * @param Y Array de valores de Y.
 * @param n Numero de puntos de datos.
 */
void print_data_points(double X[], double Y[], int n);

/**
 * @brief Funcion para evaluar el Spline Cubico en un punto x dado.
 * @param X Array de valores X (puntos de datos).
 * @param solution Array de coeficientes del Spline.
 * @param n Numero de puntos de datos.
 * @param x El valor x a evaluar.
 * @return El valor y interpolado.
 */
double evaluate_spline(double X[], double solution[], int n, double x);

/* Comentarios sobre Cuadratura de Gauss-Legendre (Resumen de uso y precision):
    2 puntos: Exacto para polinomios de grado <= 3. Ideal para funciones muy suaves.
    3 puntos: Exacto para polinomios de grado <= 5. Optimo equilibrio entre precision y eficiencia.
    4 puntos: Exacto para polinomios de grado <= 7. Para mayor curvatura y complejidad.
    5 puntos: Exacto para polinomios de grado <= 9. Para funciones altamente oscilatorias.
    6 puntos: Exacto para polinomios de grado <= 11. Para funciones extremadamente complejas.
*/

// Variables globales para almacenar los datos del Spline para la integracion
double spline_X[MAX_POINTS];
double spline_solution[MAX_SIZE + 1];
int spline_n = 0;

int main(int argc, char const *argv[]) {
    // Factores de ponderacion (pesos) y argumentos (raices) para Gauss-Legendre
    double c0, c1, c2, c3, c4, c5;
    double x0, x1, x2, x3, x4, x5;
    // Numero de puntos para la Cuadratura de Gauss-Legendre
    int number_of_points;
    // Resultado de la integral
    double integral;

    // Arrays para el calculo de los coeficientes del Spline (sistema Ax=b)
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Puntos de datos leidos del archivo de texto
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Numero de puntos de datos
    int n;
    // Limites de integracion [a, b]
    double a_limit, b_limit;
    // Opcion: funcion explicita o tabla de datos
    int choice;
    char continuar;

    printf("========================================================\n");
    printf("  INTEGRACION NUMERICA - CUADRATURA GAUSS-LEGENDRE\n");
    printf("========================================================\n");
    printf("Funcion actual: f(x) = 2*x^3\n\n");

    do {
        // Resetear spline_n para cada iteracion
        spline_n = 0;

        printf("Tiene una funcion o una tabla de datos?\n");
        printf("1. Tengo una tabla de datos\n");
        printf("2. Tengo una funcion (f(x) por defecto sera 2*x^3)\n");
        printf("--------------------------------------------------------\n");
        printf("Ingrese su eleccion: ");
        scanf("%d", &choice);

    if(choice == 1) {
        // --- PROCESAMIENTO DEL SPLINE CUBICO ---
        
        // 1. Lectura de datos
        if (!read_data_points("../data.txt", X, Y, &n)) {
            printf("Error al leer los datos del archivo. Saliendo.\n");
            continue;
        }
        print_data_points(X, Y, n);

        // 2. Inicializar la matriz A y el vector b del sistema lineal (tamaño 4*(n-1))
        for (int i = 0; i < 4*(n-1); i++) {
            b[i] = 0.0;
            for (int j = 0; j < 4*(n-1); j++) {
                A[i][j] = 0.0;
            }
        }

        // --- CONSTRUCCION DE LAS ECUACIONES DEL SISTEMA Ax=b ---

        // 2(n-1) Ecuaciones: Cada Spline S_k debe pasar por sus dos puntos finales (X_k, Y_k) y (X_{k+1}, Y_{k+1})
        for(int k = 0; k < n-1; k++) {
            // S_k(X_k) = Y_k
            for(int j = 0; j <= 3; j++) {
                A[2*k][4*k+j] = pow(X[k], 3-j); // Coeficientes: X[k]^3, X[k]^2, X[k]^1, X[k]^0
            }
            b[2*k] = Y[k];

            // S_k(X_{k+1}) = Y_{k+1}
            for(int j = 0; j <= 3; j++) {
                A[2*k+1][4*k+j] = pow(X[k+1], 3-j); // Coeficientes: X[k+1]^3, X[k+1]^2, X[k+1]^1, X[k+1]^0
            }
            b[2*k+1] = Y[k+1];
        }


        // n-2 Ecuaciones: Continuidad de las primeras derivadas en los puntos interiores (X_{k+1})
        for(int k = 0; k < n-2; k++) {
            int row = 2*(n-1) + k;
            // S'_k(X_{k+1})
            for(int j = 0; j <= 2; j++) {
                A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j); // Derivada: 3*x^2, 2*x^1, 1
            }
            // - S'_{k+1}(X_{k+1})
            for(int j = 0; j <= 2; j++) {
                A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
            }
            b[row] = 0.0;
        }

        // n-2 Ecuaciones: Continuidad de las segundas derivadas en los puntos interiores (X_{k+1})
        for(int k = 0; k < n-2; k++) {
            int row = 2*(n-1) + (n-2) + k;
            // S''_k(X_{k+1})
            A[row][4*k] = 6 * X[k+1];  // Segunda derivada: 6*x
            A[row][4*k+1] = 2;         // Segunda derivada: 2
            // - S''_{k+1}(X_{k+1})
            A[row][4*(k+1)] = -6 * X[k+1];
            A[row][4*(k+1)+1] = -2;
            b[row] = 0.0;
        }

        // 2 Ecuaciones (Condiciones de Contorno): Spline Natural (Segunda derivada = 0 en los extremos)
        int row1 = 4*(n-1) - 2;
        int row2 = 4*(n-1) - 1;

        // S''_0(X_0) = 0 (primer Spline)
        A[row1][0] = 6 * X[0];
        A[row1][1] = 2;
        b[row1] = 0.0;

        // S''_{n-2}(X_{n-1}) = 0 (ultimo Spline)
        A[row2][4*(n-2)] = 6 * X[n-1];
        A[row2][4*(n-2)+1] = 2;
        b[row2] = 0.0;

        // 3. Resolver el sistema de ecuaciones para encontrar los coeficientes del Spline
        gauss_elimination(4*(n-1), A, b, solution);

        // Guardamos los coeficientes y los puntos para usarlos en la función f(x)
        for (int i = 0; i < n; i++) {
            spline_X[i] = X[i];
        }
        for (int i = 0; i < 4*(n-1); i++) {
            spline_solution[i] = solution[i];
        }
        spline_n = n; // Indica que se usara el Spline en f(x)
    }

    // --- INTEGRACION POR CUADRATURA DE GAUSS-LEGENDRE ---

    printf("\nInserte los limites de integracion (a b):\n");
    scanf("%lf %lf", &a_limit, &b_limit);

    printf("Inserte el numero de puntos de Gauss-Legendre (entre 2 y 6):\n");
    scanf("%d", &number_of_points);

    // Transformacion de variable: integral en [a, b] a integral en [-1, 1]
    // Argumento transformado para f(x): ((b-a) * x_i + (b+a)) / 2
    // Factor de escalamiento: (b-a) / 2
    double scale_factor = (b_limit - a_limit) / 2.0;
    double offset = (b_limit + a_limit) / 2.0;

    switch(number_of_points) {
        case 2: 
            // 2-puntos: Pesos (c_i) y Raices (x_i)
            c0 = 1.0; c1 = 1.0;
            x0 = -0.577350269; x1 = 0.577350269;
            integral = scale_factor * (
                c0 * f(scale_factor * x0 + offset) + 
                c1 * f(scale_factor * x1 + offset)
            );
            break;
        case 3:
            // 3-puntos
            c0 = 0.5555556; c1 = 0.8888889; c2 = 0.5555556;
            x0 = -0.774596669; x1 = 0.0; x2 = 0.774596669;
            integral = scale_factor * (
                c0 * f(scale_factor * x0 + offset) + 
                c1 * f(scale_factor * x1 + offset) +
                c2 * f(scale_factor * x2 + offset)
            );
            break;
        case 4:
            // 4-puntos
            c0 = 0.3478548; c1 = 0.6521452; c2 = 0.6521452; c3 = 0.3478548;
            x0 = -0.861136312; x1 = -0.339981044; x2 = 0.339981044; x3 = 0.861136312;
            integral = scale_factor * (
                c0 * f(scale_factor * x0 + offset) + 
                c1 * f(scale_factor * x1 + offset) + 
                c2 * f(scale_factor * x2 + offset) +
                c3 * f(scale_factor * x3 + offset)
            );
            break;
        case 5:
            // 5-puntos
            c0 = 0.2369269; c1 = 0.4786287; c2 = 0.5688889; c3 = 0.4786287; c4 = 0.2369269;
            x0 = -0.906179846; x1 = -0.538469310; x2 = 0.0; x3 = 0.538469310; x4 = 0.906179846;
            integral = scale_factor * (
                c0 * f(scale_factor * x0 + offset) + 
                c1 * f(scale_factor * x1 + offset) + 
                c2 * f(scale_factor * x2 + offset) + 
                c3 * f(scale_factor * x3 + offset) + 
                c4 * f(scale_factor * x4 + offset)
            );
            break;
        case 6:
            // 6-puntos
            c0 = 0.1713245; c1 = 0.3607616; c2 = 0.4679139; c3 = 0.4679139; c4 = 0.3607616; c5 = 0.1713245;
            x0 = -0.932469514; x1 = -0.661209386; x2 = -0.238619186; x3 = 0.238619186; x4 = 0.661209386; x5 = 0.932469514;
            integral = scale_factor * (
                c0 * f(scale_factor * x0 + offset) + 
                c1 * f(scale_factor * x1 + offset) + 
                c2 * f(scale_factor * x2 + offset) + 
                c3 * f(scale_factor * x3 + offset) + 
                c4 * f(scale_factor * x4 + offset) + 
                c5 * f(scale_factor * x5 + offset)
            );
            break;
        default:
            printf("Error: El numero de puntos debe estar entre 2 y 6\n");
            continue;
    }

    printf("\nEl valor aproximado de la integral es: %lf\n", integral);

    printf("\n--------------------------------------------------------\n");
    printf("Desea realizar otro calculo? (s/n): ");
    scanf(" %c", &continuar);
    printf("\n");

    } while (continuar == 's' || continuar == 'S');

    printf("Gracias por usar el programa!\n");
    printf("========================================================\n");
    return 0;
}


double f(double x) {
    // Si se cargo un Spline (spline_n > 0), evaluamos el Spline en x.
    if(spline_n > 0) {
        return evaluate_spline(spline_X, spline_solution, spline_n, x);
    }
    // Si no se cargo un Spline (opcion 2), por defecto usamos una funcion de prueba.
    // Aca se ingresa la funcion a evaluar
    return 2.0 * (x * x * x); 
}


int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: No se pudo abrir el archivo '%s'\n", filename);
        return 0;
    }
    
    printf("Archivo '%s' abierto correctamente\n", filename);
    
    // Leer el numero de puntos de datos
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: No se pudo leer el numero de puntos de datos\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Numero de puntos invalido (%d). Maximo permitido: %d\n", *n, MAX_POINTS);
        fclose(fp);
        return 0;
    }
    
    // Leer los pares de datos
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el par de datos %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Se leyeron %d puntos de datos correctamente\n\n", *n);
    return 1;
}

void print_data_points(double X[], double Y[], int n) {
    printf("Puntos de Datos:\n");
    printf("=============\n");
    printf("   i   |       Xi      |       Yi      \n");
    printf("-------|---------------|---------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d   | %12.6f  | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Buscar el intervalo de Spline correcto [X_k, X_{k+1}] donde se encuentra x
    if (x <= X[0]) {
        k = 0;      // Usar el primer Spline para valores <= X[0]
    } else if (x >= X[n-1]) {
        k = n - 2;  // Usar el último Spline para valores >= X[n-1]
    } else {
        // Búsqueda lineal del intervalo
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluar S_k(x) = a*x^3 + b*x^2 + c*x + d
    // Los coeficientes (a, b, c, d) para el Spline k están en solution[4*k] a solution[4*k+3]
    double y = solution[4*k] * pow(x, 3) +
               solution[4*k+1] * pow(x, 2) +
               solution[4*k+2] * x +
               solution[4*k+3];

    return y;
}