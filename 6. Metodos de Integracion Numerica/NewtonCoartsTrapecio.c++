#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constantes para el tamano maximo de los arreglos
#define MAX_POINTS 20
#define MAX_SIZE 100

// Se asume que este archivo existe y contiene la funcion gauss_elimination
#include "gauss.h" 

/**
// Funcion para leer pares de datos Xi, Yi desde un archivo
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
    // Limites de integracion
    double a, b;
    // Integral calculada (suma)
    double sum = 0.0;
    // Valores intermedios para dividir el intervalo de integracion [a,b]
    double x = 0.0;
    // Distancia entre dos puntos consecutivos (paso h)
    double h = 0.0;
    char continuar;

    printf("========================================================\n");
    printf("    METODOS DE INTEGRACION NUMERICA - TRAPECIO\n");
    printf("========================================================\n");
    printf("Funcion actual: f(x) = x^2 + 1\n\n");

    do {
        printf("Elija una opcion:\n");
        printf("1. Trapecio Compuesto\n");
        printf("2. Trapecio Simple\n");
        printf("--------------------------------------------------------\n");
        printf("Ingrese su eleccion: ");
        scanf("%d", &opcion);

    if (opcion == 1) {
        printf("\nTiene una funcion o una tabla de datos?\n");
        printf("1. Tengo una funcion\n");
        printf("2. Tengo una tabla de datos\n");
        printf("Ingrese su eleccion: ");
        scanf("%d", &opcion);
        
        switch (opcion) {
            case 1:
                printf("\nInserte los limites de integracion:\n");
                scanf("%lf %lf", &a, &b);
                printf("Por favor ingrese el numero de subintervalos:\n");
                scanf("%d", &subintervalos);

                // ============================================================
                // TRAPECIO COMPUESTO CON FUNCION
                // Formula:
                // I = h/2 * [ f(a) + f(b) + 2*(f(x1) + ... + f(xn-1)) ]
                // con h = (b-a)/n
                // ============================================================

                // PASO 1 (extremos): suma f(a) + f(b)
                sum = f(a) + f(b);

                // PASO 2 (paso): h = (b-a)/n
                h = (b - a) / subintervalos;

                // PASO 3 (sumatoria interna): recorre los puntos x1..x(n-1)
                for(int i = 1; i <= subintervalos-1; i++) {
                    // x_i = a + i*h: punto interno del intervalo
                    x = a + i * h;
                    // Acumula 2*f(x_i) en la suma
                    sum += 2 * f(x);
                }

                // PASO 4 (integral): I = h/2 * suma
                // (b-a)/(2n) es equivalente a h/2
                sum = ((b-a)/(2*subintervalos)) * sum;
                
                // Imprimir el resultado de la integral
                printf("La integral es: %lf\n", sum);
                break;

            case 2: {
                // Arreglos para el cálculo de coeficientes del polinomio a usar para el Spline
                double A[MAX_SIZE+1][MAX_SIZE+1], b_vec[MAX_SIZE+1], solution[MAX_SIZE+1];
                // Puntos de datos a leer desde el archivo de texto
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Puntos de datos para calcular la integral
                double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                // Número de puntos
                int n;
                // Bandera para saber si los puntos están equiespaciados
                int equiespaciados = 1;

                // PASO 1 (lectura): leer X[i], Y[i] desde el archivo
                // El archivo tiene: primera linea = n, luego n pares "x y"
                if (!read_data_points("../data.txt", X, Y, &n)) {
                    printf("Fallo al leer datos del archivo. Saliendo.\n");
                    continue;
                }
                // Imprimir los puntos de datos
                print_data_points(X, Y, n);

                // Se necesitan al menos 2 puntos para poder integrar
                if (n < 2) {
                    printf("Se necesitan al menos 2 puntos para integrar.\n");
                    break;
                }

                // PASO 2 (paso promedio): h = (x[n-1] - x[0]) / (n-1)
                h = (X[n - 1] - X[0]) / (n - 1);

                // PASO 3 (verificar equiespaciado):
                // Se compara cada diferencia X[i+1]-X[i] contra el paso h.
                // Si alguna difiere, los puntos NO estan equiespaciados.
                for (int i = 0; i < n - 1; i++) {
                    double hi = X[i + 1] - X[i];
                    if (fabs(hi - h) > 1e-8 * (1.0 + fabs(h))) {
                        printf("Los puntos no están equiespaciados. Se interpolará con Spline.\n");
                        equiespaciados = 0;
                        break;
                    }
                }

                if (equiespaciados) {
                    // PASO 4a (equiespaciados): usar directamente la tabla
                    for (int i = 0; i < n; i++) {
                        new_X[i] = X[i];
                        new_Y[i] = Y[i];
                    }
                } else {
                    // ============================================================
                    // PASO 4b (NO equiespaciados): interpolar con Spline cubico natural
                    // Se arma el sistema A*solution = b_vec y se resuelve con Gauss.
                    // Luego se evalua el spline en puntos equiespaciados.
                    // ============================================================

                    // Inicializar A y b para el Spline
                    for (int i = 0; i < 4*(n-1); i++) {
                        b_vec[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)] para el Spline
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos puntos finales
                    for(int k = 0; k < n-1; k++) {
                        // El Spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b_vec[2*k] = Y[k];
                    
                        // El Spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b_vec[2*k+1] = Y[k+1];
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de las primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada del spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada del spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b_vec[row] = 0.0;
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de las segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada del spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada del spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b_vec[row] = 0.0;
                    }

                    // Dos condiciones de contorno: Spline natural (segundas derivadas = 0 en los puntos finales)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primer spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b_vec[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (último spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b_vec[row2] = 0.0;

                    // Resolver el sistema lineal A*solution = b_vec con eliminacion Gaussiana
                    gauss_elimination(4*(n-1), A, b_vec, solution);

                    // PASO 5 (nuevos puntos): xp[i] = X[0] + i*h
                    // yp[i] = Spline(xp[i]): se evalua el spline en cada punto nuevo
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }
                }

                // ============================================================
                // PASO 6 (integral): Trapecio Compuesto sobre los puntos obtenidos
                // I = h/2 * [ yp[0] + yp[n-1] + 2*(yp[1] + ... + yp[n-2]) ]
                // ============================================================

                // Suma de los extremos
                sum = new_Y[0] + new_Y[n - 1];
                // Sumatoria de los puntos internos, multiplicados por 2
                for(int i = 1; i <= n - 2; i++) {
                    sum += 2 * new_Y[i];
                }
                // Multiplicacion final por h/2
                sum = (h / 2.0) * sum;

                printf("La integral es: %lf\n", sum);
                break;
            }
            default:
                printf("Opcion invalida\n");
                break;
        }
    } else if(opcion == 2) {
        // ============================================================
        // TRAPECIO SIMPLE + CALCULO DE ERRORES
        // Integral: I = (b-a) * (f(a) + f(b)) / 2
        // Error aproximado: |-(1/12) * f''(c) * (b-a)^3|
        // ============================================================
        // Punto para calcular el error
        double c;
        // Integral del Trapecio Simple
        double Iaprox = 0.0;
        // Error aproximado del Trapecio Simple
        double aprox_error = 0.0;
        // Error exacto y error porcentual
        double exact_error = 0.0, porcentual_error = 0.0;
        // Integral Exacta (necesita ser calculada manualmente)
        double Iexact = 0.0;

        printf("\nInserte los limites de integracion:\n");
        scanf("%lf %lf", &a, &b);
        printf("Inserte un valor dentro del intervalo [a,b] para calcular el error (c):\n");
        scanf("%lf", &c);
        printf("Inserte el valor exacto de la integral para calcular el error exacto:\n");
        scanf("%lf", &Iexact);

        // Integral del trapecio simple: (b-a)*(f(a)+f(b))/2
        Iaprox = (b-a) * ((f(b) + f(a)) / 2.0);
        // Error aproximado: |-(1/12)*f''(c)*(b-a)^3|
        aprox_error = fabs(-(1.0/12.0) * second_derivative(f, c) * pow((b - a), 3));
        // Error exacto: |Iexacta - Iaprox|
        exact_error = fabs(Iexact - Iaprox);
        // Error porcentual: |Iexacta - Iaprox| / |Iexacta| * 100
        porcentual_error = (fabs(Iexact - Iaprox) / fabs(Iexact)) * 100.0;

        // Imprimir resultados
        printf("La integral aproximada es: %lf\n", Iaprox);
        printf("El error aproximado es: %lf\n", aprox_error);
        printf("El error exacto es: %lf\n", exact_error);
        printf("El error porcentual es: %lf%%\n", porcentual_error);

    } else {
        printf("Opcion invalida\n");
    }

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
// Aqui se define la curva que se va a integrar.
// Cambiar esta linea si se quiere integrar otra funcion.
double f(double x) {
    return (x * x) + 1.0;
}

// Implementación de la segunda derivada por diferencia central
// Formula: f''(x) ~= [f(x+h) - 2*f(x) + f(x-h)] / h^2
double second_derivative(double (*func)(double), double x, double h) {
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
}

// Implementación de la cuarta derivada por diferencia central
// Formula: f''''(x) ~= [f(x-2h) - 4f(x-h) + 6f(x) - 4f(x+h) + f(x+2h)] / h^4
// Nota: en este archivo no se usa, queda disponible para otros metodos.
double fourth_derivative(double (*func)(double), double x, double h) {
    return (func(x - 2*h) - 4*func(x - h) + 6*func(x) - 4*func(x + h) + func(x + 2*h)) / (pow(h, 4));
}


// Implementacion de la lectura de puntos
// Formato del archivo:
//   primera linea: n (cantidad de puntos)
//   luego n lineas: x y
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
    
    // Leer los puntos de datos (bucle de lectura de X[i], Y[i])
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

    // PASO 1: encontrar el intervalo [X[k], X[k+1]] que contiene a x
    if (x <= X[0]) {
        k = 0;  // Usar el primer spline para x menor o igual al primer punto
    } else if (x >= X[n-1]) {
        k = n - 2; // Usar el último spline para x mayor o igual al último punto
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // PASO 2: evaluar el polinomio del tramo k
    // S_k(x) = a*x^3 + b*x^2 + c*x + d
    // solution[4k], solution[4k+1], solution[4k+2], solution[4k+3] son a, b, c, d
    double y = solution[4*k] * pow(x, 3) +
                solution[4*k+1] * pow(x, 2) +
                solution[4*k+2] * x +
                solution[4*k+3];

    return y;
}
