#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constantes para el tamano maximo de los arreglos
#define MAX_POINTS 20
#define MAX_SIZE 100

// Se asume que este archivo existe y contiene la funcion gauss_elimination
#include "gauss.h" 

/**
 * Funcion para leer pares de datos Xi, Yi desde un archivo
 * Formato esperado: la primera linea contiene el numero de puntos, luego pares xi yi
 * @param filename Nombre del archivo a leer
 * @param X Arreglo para almacenar los valores de X
 * @param Y Arreglo para almacenar los valores de Y
 * @param n Puntero para almacenar el numero de puntos de datos leidos
 * @return 1 si tiene exito, 0 en caso contrario
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * Funcion para mostrar los puntos de datos
 * @param X Arreglo de valores X
 * @param Y Arreglo de valores Y
 * @param n Numero de puntos de datos
 */
void print_data_points(double X[], double Y[], int n);

/**
 * Funcion para definir la funcion f(x) a integrar
 * @param x El punto en el que evaluar la funcion 
 * @return El valor de la funcion en x
 */
double f(double x);

/**
 * Funcion para evaluar el spline cubico en un punto dado x
 * @param X Arreglo de valores X (puntos de datos)
 * @param solution Arreglo de coeficientes del spline
 * @param n Numero de puntos de datos
 * @param x El valor x a evaluar
 * @return El valor y interpolado
 */
double evaluate_spline(double X[], double solution[], int n, double x);


/** * Segunda derivada usando aproximacion de diferencias finitas
 * @param func Puntero a la funcion
 * @param x El punto en el que evaluar la segunda derivada
 * @param h Un valor pequeno para la aproximacion de diferencias finitas (por defecto es 1e-5)
 * @return La segunda derivada de la funcion en el punto x
 * */ 
double second_derivative(double (*func)(double), double x, double h = 1e-5);

/** * Cuarta derivada usando aproximacion de diferencias finitas
 * @param func Puntero a la funcion
 * @param x El punto en el que evaluar la cuarta derivada
 * @param h Un valor pequeno para la aproximacion de diferencias finitas (por defecto es 1e-5)
 * @return La cuarta derivada de la funcion en el punto x
 * */
double fourth_derivative(double (*func)(double), double x, double h = 0.01);

int main(int argc, char const *argv[]) {
    int opcion, subintervalos;
    char continuar;
    // Limites de integracion
    double a, b;
    // Integral calculada (suma)
    double sum = 0.0;
    // Valores intermedios para dividir el intervalo de integracion [a,b]
    double x = 0.0;
    // Distancia entre dos puntos consecutivos (paso h)
    double h = 0.0;

    printf("========================================================\n");
    printf("    METODOS DE INTEGRACION NUMERICA - SIMPSON\n");
    printf("========================================================\n");
    printf("Funcion actual: f(x) = x^2 + 1\n\n");

    do {
        printf("Elija una opcion:\n");
        printf("1. Simpson Compuesto\n");
        printf("2. Simpson 1/3\n");
        printf("--------------------------------------------------------\n");
        printf("Ingrese su eleccion: ");
    scanf("%d", &opcion);

    if(opcion == 1) {
        printf("\nTiene una funcion o una tabla de datos?\n");
        printf("1. Tengo una funcion\n");
        printf("2. Tengo una tabla de datos\n");
        printf("Ingrese su eleccion: ");
        scanf("%d", &opcion);
        
        switch(opcion) {
            case 1: {
                // Valores para dividir el intervalo de integracion [a,b]
                double x_arr[MAX_POINTS];
                printf("\nInserte los limites de integracion (a b): ");
                scanf("%lf %lf", &a, &b);
                printf("Por favor ingrese el numero de subintervalos (debe ser un numero par): ");
                scanf("%d", &subintervalos);
                
                if(subintervalos % 2 != 0) {
                    printf("ERROR: El numero de subintervalos debe ser PAR.\n");
                    continue;
                }
                
                // Calcular I con la regla de Simpson Compuesto
                sum = f(a) + f(b);
                h = (b - a) / subintervalos;        
                for(int i = 1; i < subintervalos; i++) {
                    x_arr[i] = a + (i * h);
                    if (i % 2 == 0) {
                        sum += 2 * f(x_arr[i]);  // índices pares (múltiplos de 2)
                    } else {
                        sum += 4 * f(x_arr[i]);  // índices impares
                    }
                }
                sum = (h/3) * sum;

                // Imprimir integral usando Simpson Compuesto
                printf("La integral es: %lf\n", sum);
                break;
            }
            case 2: {
                // Arreglos para el calculo de coeficientes del polinomio a usar para el Spline
                double A[MAX_SIZE+1][MAX_SIZE+1], b_vec[MAX_SIZE+1], solution[MAX_SIZE+1];
                // Puntos de datos a leer desde el archivo de texto
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Puntos de datos para calcular la integral
                double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                // Numero de puntos
                int n;

                // Leer puntos de datos desde archivo
                if (!read_data_points("../data.txt", X, Y, &n)) {
                    printf("Fallo al leer datos del archivo.\n");
                    continue;
                }
                // Imprimir los puntos de datos
                print_data_points(X, Y, n);

                // Inicializar A y b
                for (int i = 0; i < 4*(n-1); i++) {
                    b_vec[i] = 0.0;
                    for (int j = 0; j < 4*(n-1); j++) {
                        A[i][j] = 0.0;
                    }
                }

                // Determinar el número de subintervalos. Debe ser par para Simpson.
                // Si n-1 (el número de intervalos iniciales) es par, lo usamos.
                // Si n-1 es impar, usamos n (crearemos n puntos igualmente espaciados, resultando en n-1 intervalos)
                if((n-1) % 2 == 0) {
                    subintervalos = n-1;
                } else {
                    subintervalos = n;
                }

                // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)]
                // (Mismas ecuaciones de contorno del Spline que en el caso 1)
                for(int k = 0; k < n-1; k++) {
                    for(int j = 0; j <= 3; j++) {
                        A[2*k][4*k+j] = pow(X[k], 3-j);
                    }
                    b_vec[2*k] = Y[k];
                
                    for(int j = 0; j <= 3; j++) {
                        A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                    }
                    b_vec[2*k+1] = Y[k+1];
                }

                // Siguientes n-2 ecuaciones: Continuidad de las primeras derivadas
                for(int k = 0; k < n-2; k++) {
                    int row = 2*(n-1) + k;
                    for(int j = 0; j <= 2; j++) {
                        A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                    }
                    for(int j = 0; j <= 2; j++) {
                        A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                    }
                    b_vec[row] = 0.0;
                }

                // Siguientes n-2 ecuaciones: Continuidad de las segundas derivadas
                for(int k = 0; k < n-2; k++) {
                    int row = 2*(n-1) + (n-2) + k;
                    A[row][4*k] = 6 * X[k+1];
                    A[row][4*k+1] = 2;
                    A[row][4*(k+1)] = -6 * X[k+1];
                    A[row][4*(k+1)+1] = -2;
                    b_vec[row] = 0.0;
                }

                // Dos condiciones de contorno: Spline natural
                int row1 = 4*(n-1) - 2;
                int row2 = 4*(n-1) - 1;

                // Segunda derivada = 0 en X[0]
                A[row1][0] = 6 * X[0];
                A[row1][1] = 2;
                b_vec[row1] = 0.0;

                // Segunda derivada = 0 en X[n-1]
                A[row2][4*(n-2)] = 6 * X[n-1];
                A[row2][4*(n-2)+1] = 2;
                b_vec[row2] = 0.0;

                // Resolver el sistema
                gauss_elimination(4*(n-1), A, b_vec, solution);


                // Dividir el intervalo [X[0], X[n-1]] en 'subintervalos' de igual longitud
                h = (X[n - 1] - X[0]) / subintervalos;

                // Calcular X[i] y Y[i] igualmente espaciados
                for(int i = 0; i <= subintervalos; i++) {
                    new_X[i] = X[0] + i * h;
                    new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                }

                // Calcular I aplicando la regla de Simpson Compuesto
                sum = new_Y[0] + new_Y[subintervalos];
                for(int i = 1; i < subintervalos; i++) {
                    if(i % 2 == 0) {
                        sum += 2 * new_Y[i]; // Coeficiente 2 para índices pares
                    } else {
                        sum += 4 * new_Y[i]; // Coeficiente 4 para índices impares
                    }
                }
                sum = (h/3) * sum;

                // Imprimir integral usando Simpson Compuesto con puntos igualmente espaciados creados
                printf("La integral es: %lf\n", sum);
                break;
            }
            default:
                printf("Opcion invalida\n");
                continue;
        }
    } else if(opcion == 2) {
        // Regla de Simpson 1/3 Simple. Funciona perfectamente con polinomios de grado 3 o menor
        // Punto para calcular el error
        double c;
        // Integral de Simpson 1/3 aproximada
        double Iaprox = 0.0;
        // Error aproximado de Simpson 1/3
        double aprox_error = 0.0;
        // Error exacto y error porcentual
        double exact_error = 0.0, porcentual_error = 0.0;
        // Integral Exacta (necesita ser calculada manualmente)
        double Iexact = 0.0;

        printf("\nInserte los limites de integracion (a b): ");
        scanf("%lf %lf", &a, &b);
        printf("Inserte un valor dentro del intervalo [a,b] para calcular el error (c): ");
        scanf("%lf", &c);
        printf("Inserte el valor exacto de la integral para calcular el error exacto: ");
        scanf("%lf", &Iexact);

        // Calcular I
        Iaprox = ((b-a)/6.0) * (f(a) + 4.0*f((a+b)/2.0) + f(b));
        aprox_error = fabs(-(1.0/2880.0) * pow(b-a, 5) * fourth_derivative(f, c));
        exact_error = fabs(Iexact - Iaprox);
        porcentual_error = (fabs(Iexact - Iaprox) / fabs(Iexact)) * 100.0;

        // Imprimir resultados
        printf("La integral aproximada es: %lf\n", Iaprox);
        printf("El error aproximado es: %lf\n", aprox_error);
        printf("El error exacto es: %lf\n", exact_error);
        printf("El error porcentual es: %lf%%\n", porcentual_error);
    } else {
        printf("ERROR: Debe insertar un numero en el intervalo [1, 2]\n");
        continue;
    }
    
    // Preguntar si desea realizar otro calculo
    printf("\n--------------------------------------------------------\n");
    printf("Desea realizar otro calculo? (s/n): ");
    scanf(" %c", &continuar);
    printf("\n");
    
    } while (continuar == 's' || continuar == 'S');
    
    printf("Gracias por usar el programa!\n");
    printf("========================================================\n");
    
    return 0;
}

// Funcion para definir f(x)
double f(double x) {
    // Ejemplo de funcion: sin(2x) * e^(-x)
    return (x * x) + 1.0;
}

// Implementacion de la segunda derivada por diferencia central
double second_derivative(double (*func)(double), double x, double h) {
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
}

// Implementacion de la cuarta derivada por diferencia central
double fourth_derivative(double (*func)(double), double x, double h) {
    return (func(x - 2*h) - 4*func(x - h) + 6*func(x) - 4*func(x + h) + func(x + 2*h)) / (pow(h, 4));
}


// Implementacion de la lectura de puntos
int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", filename);
        return 0;
    }
    
    printf("Archivo '%s' abierto exitosamente\n", filename);
    
    // Leer el numero de puntos de datos
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: No se puede leer el numero de puntos de datos\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Numero de puntos invalido (%d)\n", *n);
        fclose(fp);
        return 0;
    }
    
    // Leer los puntos de datos
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se puede leer el punto de datos %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Se leyeron %d puntos de datos exitosamente\n\n", *n);
    return 1;
}

// Implementacion de la impresion de puntos
void print_data_points(double X[], double Y[], int n) {
    printf("Puntos de Datos:\n");
    printf("=============\n");
    printf("   i  |      Xi      |      Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

// Implementacion de la evaluacion del spline
double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Encontrar el intervalo de spline correcto
    if (x <= X[0]) {
        k = 0;  // Usar el primer spline para x menor o igual al primer punto
    } else if (x >= X[n-1]) {
        k = n - 2; // Usar el ultimo spline para x mayor o igual al ultimo punto
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluar S_k(x) = a*x^3 + b*x^2 + c*x + d
    double y = solution[4*k] * pow(x, 3) +
                solution[4*k+1] * pow(x, 2) +
                solution[4*k+2] * x +
                solution[4*k+3];

    return y;
}