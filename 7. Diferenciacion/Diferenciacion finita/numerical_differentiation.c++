#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definiciones de constantes
#define MAX_POINTS 20 // Número máximo de puntos de datos que puede manejar el programa
#define MAX_SIZE 100  // Tamaño máximo general para arreglos
#define PI 3.14159265358979323846

// Se asume que este archivo contiene la implementación de métodos numéricos
// como la Eliminación Gaussiana, necesaria para calcular los coeficientes del spline.
#include "gauss.h"

/**
 * Función para leer pares de datos Xi, Yi desde un archivo
 * Formato esperado: la primera línea contiene el número de puntos, luego pares xi yi
 * @param nombre_archivo Nombre del archivo a leer
 * @param X Arreglo para almacenar los valores X (Puntos nodales)
 * @param Y Arreglo para almacenar los valores Y (Valores de la función en los nodos)
 * @param n Puntero para almacenar el número de puntos de datos leídos
 * @return 1 si tiene éxito, 0 en caso contrario
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * Función para mostrar los puntos de datos
 * @param X Arreglo de valores X
 * @param Y Arreglo de valores Y
 * @param n Número de puntos de datos
 */
void print_data_points(double X[], double Y[], int n);

/**
 * Función para definir la función f(x)
 * @param x El punto en el que se evalúa la función
 * @return El valor de la función en x
 */
double f(double x);

/**
 * Función para guardar los valores calculados en un archivo de texto
 * @param x Arreglo de valores x (posiciones)
 * @param derivada_primera Arreglo de valores de la primera derivada (fp significa f prime o f prima)
 * @param n Número de subintervalos
 */
void save_in_txt(double x[], double fp[], int n);
/**
 * Función para evaluar el spline cúbico en un punto x dado
 * @param X_nodos Arreglo de los valores X (puntos nodales)
 * @param solucion_spline Arreglo de coeficientes del spline (a, b, c, d para cada intervalo)
 * @param n Número de puntos de datos
 * @param x El valor x a evaluar
 * @return El valor y interpolado por el spline en x
 */
double evaluate_spline(double X[], double solution[], int n, double x);

int main(int argc, char const *argv[]) {
    // Interval extremes
    double a, b;    
    // Número de subintervalos
    int n;
    // Distancia entre puntos
    double h;
    // Valor de la primera derivada
    double fp[MAX_SIZE + 1];
    // Valor de la segunda derivada
    double fpp[MAX_SIZE + 1];
    // Valor de la tercera derivada
    double fppp[MAX_SIZE + 1];
    // Valores del dominio
    double x[MAX_SIZE + 1];
    // Tipo de error O(h) o O(h²). Recuerda que los operadores de diferencia finita centrados tienen error O(h²) para O(h) y tienen error O(h⁴) para O(h²)
    int error_type;
    int number_of_derivative;
    // Variable para almacenar si es una función o una tabla de datos
    int function_data_table;
    
    printf("Insertar la derivada que quieres calcular: 1. Primera derivada 2. Segunda derivada 3. Tercera derivada\n");
    scanf("%d", &number_of_derivative);
    
    printf("Insertar el error que quieres usar: 1. O(h) 2. O(h²)\n");
    scanf("%d", &error_type);

    printf("¿Tienes una función o una tabla de datos?\n");
    printf("1. Tengo una funcion\n");
    printf("2. Tengo una tabla de datos\n");
    scanf("%d", &function_data_table);

    if(error_type == 1) {
        switch(number_of_derivative) {
            case 1: 
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar número de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular la distancia entre puntos
                    h = (b - a) / n;
                    // Calcular la primera derivada en los puntos extremos 
                    fp[0] = ((f(a+h) - f(a))/h); // This is O(h)
                    fp[n] = ((f(b) - f(b-h))/h); // This is O(h)
                
                    x[0] = a;
                    x[n] = b;
                
                    // Bucle para los puntos internos
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fp[i] = (f(x[i]+h) - f(x[i]-h)) / (2*h); // This is O(h²)
                    }
    
                    // Print results
                    printf("x\t\tf'(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fp[i]);
                    }
                
                    // Guardar x[i] y fp[i] en un archivo de texto
                    save_in_txt(x, fp, n);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Matrices para el cálculo de coeficientes polinómicos para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Número de puntos
                    int n;

                    // Read data points from file
                     if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("No se pudo leer los datos del archivo. Saliendo.\n");
                        return 1;
                    }
                    // Print the data points
                    print_data_points(X, Y, n);
                    /*
                    // 2. Inicializar A y b
                    for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // We calculate A[4(n-1)][4(n-1)] and b[4(n-1)]
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos extremos
                    for(int k = 0; k < n-1; k++) {
                        // La spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // Spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
                    }


                    // Siguientes n-2 ecuaciones: Continuidad de primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada de la spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada de la spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b[row] = 0.0;
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada de la spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada de la spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b[row] = 0.0;
                    }

                    // Dos condiciones de frontera: Spline natural (segundas derivadas = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primera spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (última spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Usamos la función de gauss.h para resolver el sistema con eliminación gaussiana
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Dividir el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1);

                    // Calcular X[n] y Y[n] igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }

                    // Inicializar arreglo x con los nuevos puntos igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    } */

                    h = X[1] - X[0]; // Asumiendo un espaciado uniforme en los datos originales

                    // Calcular primera derivada en los puntos extremos 
                    // fp[0] = ((f(a+h) - f(a))/h); // This is O(h)
                    fp[0] = ((Y[1] - Y[0])/h);
                    // fp[n-1] = ((f(b) - f(b-h))/h); // This is O(h)
                    fp[n-1] = ((Y[n-1] - Y[n-2])/h);
                
                
                    // Para los puntos internos
                    for(int i = 1; i <= n-2; i++) {
                        // x[i] = a + i*h;
                        // fp[i] = (f(x[i]+h) - f(x[i]-h)) / (2*h); // This is O(h²)
                        fp[i] = (Y[i+1] - Y[i-1]) / (2*h); // This is O(h²)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'(x)\n");
                    for(int i = 0; i < n; i++) {  
                        printf("%lf\t%lf\n", X[i], fp[i]);
                    }
                
                    // Guardar x[i] y fp[i] en un archivo de texto
                    save_in_txt(X, fp, n-1); //We have n points
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                }
                break; 
            case 2:
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar número de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular distancia entre puntos
                    h = (b - a) / n;

                    // Calcular segunda derivada en los puntos extremos
                    fpp[0] = (f(a+2*h) - 2*f(a+h) + f(a))/(h*h); // This is O(h)
                    fpp[n] = (f(b) - 2*f(b-h) + f(b-2*h))/(h*h); // This is O(h)
    
                    x[0] = a;
                    x[n] = b;
                
                    // Para los puntos internos
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fpp[i] = (f(x[i]+h) - 2*f(x[i]) + f(x[i]-h))/(h*h); // This is O(h²)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf''(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fpp[i]);
                    }
    
                    // Guardar x[i] y fpp[i] en un archivo de texto
                    save_in_txt(x, fpp, n);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Arreglos para el cálculo de coeficientes polinomiales para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde el archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Number of points
                    int n;

                    // Leer puntos de datos desde el archivo
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("No se pudo leer los datos del archivo. Saliendo.\n");
                        return 1;
                    }
                    // Imprimir los puntos de datos
                    print_data_points(X, Y, n);

                    // 2. Inicializar A y b
                    /* for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)].
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos extremos
                    for(int k = 0; k < n-1; k++) {
                        // Spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // Spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
                    }


                    // Next n-2 ecuaciones: Continuidad de primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada de spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada de spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b[row] = 0.0;
                    }

                    // Next n-2 ecuaciones: Continuidad de segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada de spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada de spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b[row] = 0.0;
                    }

                    // Dos condiciones de frontera: Spline natural (segundas derivadas = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primera spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (última spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Utilizamos la función de gauss.h para resolver el sistema mediante eliminación gaussiana.
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Dividir el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1);

                    // Calcular X[n] y Y[n] equidistantes
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }

                    // Inicializar la matriz x con los nuevos puntos equidistantes
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    } */

                    h = X[1] - X[0]; // Asumiendo espaciamiento uniforme en los datos originales

                    // Calcular segunda derivada en los puntos extremos
                    // fpp[0] = (f(a+2*h) - 2*f(a+h) + f(a))/(h*h); // This is O(h)
                    fpp[0] = (Y[2] - 2*Y[1] + Y[0])/(h*h); // This is O(h)
                    // fpp[n] = (f(b) - 2*f(b-h) + f(b-2*h))/(h*h); // This is O(h)
                    fpp[n-1] = (Y[n-1] - 2*Y[n-2] + Y[n-3])/(h*h); // This is O(h)
    
                    // Para los puntos internos (evitando acceso fuera de rango)
                    for(int i = 1; i <= n-2; i++) {  // Changed n-1 to n-2
                        // x[i] = a + i*h;
                        // fpp[i] = (f(x[i]+h) - 2*f(x[i]) + f(x[i]-h))/(h*h); // This is O(h²)
                        fpp[i] = (Y[i+1] - 2*Y[i] + Y[i-1])/(h*h); // This is O(h²)
                    }

                    // Imprimir resultados
                    printf("x\t\tf''(x)\n");
                    for(int i = 0; i < n; i++) {
                        printf("%lf\t%lf\n", X[i], fpp[i]);
                    }
    
                    // Guardar x[i] y fpp[i] en un archivo de texto
                    save_in_txt(X, fpp, n-1);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                }
                break; 
            case 3:
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar número de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular distancia entre puntos
                    h = (b - a) / n;

                    // Calcular tercera derivada en los puntos extremos
                    fppp[0] = (f(a+3*h) - 3*f(a+2*h) + 3*f(a+h) - f(a))/(h*h*h); // This is O(h)
                    fppp[n] = (f(b) - 3*f(b-h) + 3*f(b-2*h) -f(b-3*h) )/(h*h*h); // This is O(h)
    
                    x[0] = a;
                    x[n] = b;
                
                    // Para los puntos internos
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fppp[i] = (f(x[i]+2*h) - 2*f(x[i]+h) + 2*f(x[i]-h) - f(x[i]-2*h)) / (2*(h*h*h)); // This is O(h²)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'''(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fppp[i]);
                    }
                
                    // Guardar x[i] y fppp[i] en un archivo de texto
                    save_in_txt(x, fppp, n);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Arreglos para el cálculo de coeficientes polinomiales para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde el archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Number of points
                    int n;

                    // Leer puntos de datos desde el archivo
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("No se pudieron leer los datos del archivo. Saliendo..\n");
                        return 1;
                    }
                    // Imprimir los puntos de datos
                    print_data_points(X, Y, n);

                    // 2. Inicializar A y b
                    /* for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)].
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos extremos
                    for(int k = 0; k < n-1; k++) {
                        // La spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // Spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
                    }


                    // Siguientes n-2 ecuaciones: Continuidad de primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada de la spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada de la spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b[row] = 0.0;
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada de la spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada de la spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b[row] = 0.0;
                    }

                    // Dos condiciones de frontera: Spline natural (segundas derivadas = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primera spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (última spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Usamos la función de gauss.h para resolver el sistema con eliminación gaussiana
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Dividir el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1);

                    // Calcular X[n] y Y[n] igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }

                    // Inicializar el arreglo x con los nuevos puntos igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    } */

                    // Usar el código de arriba si los datos originales no están uniformemente espaciados
                    h = X[1] - X[0]; // Suponiendo un espaciado uniforme en los datos originales

                    // Calcular la tercera derivada en los puntos extremos
                    // fppp[0] = (f(a+3*h) - 3*f(a+2*h) + 3*f(a+h) - f(a))/(h*h*h); // This is O(h)
                    fppp[0] = (Y[3] - 3*Y[2] + 3*Y[1] - Y[0])/(h*h*h); // This is O(h)
                    // fppp[n] = (f(b) - 3*f(b-h) + 3*f(b-2*h) -f(b-3*h) )/(h*h*h); // This is O(h)
                    fppp[n-1] = (Y[n-1] - 3*Y[n-2] + 3*Y[n-3] - Y[n-4])/(h*h*h); // This is O(h)

                    // Para los puntos interiores
                    for(int i = 1; i <= n-2; i++) {
                        // x[i] = a + i*h;
                        // fppp[i] = (f(x[i]+2*h) - 2*f(x[i]+h) + 2*f(x[i]-h) - f(x[i]-2*h)) / (2*(h*h*h)); // This is O(h²)
                        fppp[i] = (Y[i+2] - 2*Y[i+1] + 2*Y[i-1] - Y[i-2]) / (2*(h*h*h)); // This is O(h²)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'''(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", X[i], fppp[i]);
                    }
                
                    // Guardar x[i] y fppp[i] en un archivo de texto
                    save_in_txt(X, fppp, n);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                }

        }
    } else if(error_type == 2) {
        switch(number_of_derivative) {
            case 1: 
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar numero de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular la distancia entre puntos
                    h = (b - a) / n;

                    // Calcular la primera derivada en los puntos extremos
                    fp[0] = (-f(a+2*h) + 4*f(a+h) - 3*f(a))/(2*h); // This is O(h²)
                    fp[n] = (3*f(b) - 4*f(b-h) + f(b-2*h))/(2*h); // This is O(h²)
                
                    x[0] = a;
                    x[n] = b;
                
                    // Para los puntos interiores
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fp[i] = (-f(x[i]+2*h) + 8*f(x[i]+h) - 8*f(x[i]-h) + f(x[i]-2*h))/(12*h); // This is O(h⁴)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fp[i]);
                    }
                
                    // Guardar x[i] y fp[i] en un archivo de texto
                    save_in_txt(x, fp, n);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Arreglos para el cálculo de coeficientes polinomiales para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // número de puntos
                    int n;

                    // Leer puntos de datos desde el archivo
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("No se pudo leer los datos del archivo. Saliendo.\n");
                        return 1;
                    }
                    // Imprime los puntos de datos
                    print_data_points(X, Y, n);

                    // 2. Inicializar A y b
                    /* for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)].
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos extremos
                    for(int k = 0; k < n-1; k++) {
                        // La spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // La spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
                    }


                    // Siguientes n-2 ecuaciones: Continuidad de las primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada de la spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada de la spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b[row] = 0.0;
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de las segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada de la spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada de la spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b[row] = 0.0;
                    }

                    // Dos condiciones de contorno: Spline natural (segunda derivada = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primera spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (última spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Usamos la función de gauss.h para resolver el sistema con eliminación gaussiana
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Dividir el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1); 
                    

                    // Calcular X[n] y Y[n] igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }

                    // Inicializar el arreglo x con los nuevos puntos igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    }

                    // Inicializar el arreglo x con los puntos de datos originales
                    for(int i = 0; i < n; i++) {
                        x[i] = X[i];
                    } */

                    // Utilice el código anterior si los datos originales no tienen un espaciado uniforme.
                    h = X[1] - X[0]; // Asumiendo espaciamiento uniforme en los datos originales
                    
                    // Calcular la primera derivada con error O(h²)
                    // Diferencia hacia adelante O(h²) para el primer punto (necesita 3 puntos)
                    fp[0] = (-Y[2] + 4*Y[1] - 3*Y[0])/(2*h);
                    
                    // Diferencia hacia atrás O(h²) para el último punto (necesita 3 puntos)
                    fp[n-1] = (3*Y[n-1] - 4*Y[n-2] + Y[n-3])/(2*h);
                
                    // Calcular las derivadas de todos los puntos interiores con las fórmulas apropiadas.
                    for(int i = 1; i <= n-2; i++) {
                        // if(i >= 2 && i <= n-3) {
                            // Diferencia central O(h⁴) cuando tenemos 5 puntos disponibles
                        fp[i] = (-Y[i+2] + 8*Y[i+1] - 8*Y[i-1] + Y[i-2])/(12*h);
                        // } else {
                            // Diferencia central O(h²) cuando solo tenemos 3 puntos
                            // fp[i] = (Y[i+1] - Y[i-1])/(2*h);
                        // }
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'(x)\n");
                    for(int i = 0; i < n; i++) {
                        printf("%lf\t%lf\n", X[i], fp[i]);
                    }
                
                    // Guardar x[i] y fp[i] en un archivo de texto
                    save_in_txt(X, fp, n-1);
                
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                }
            case 2:
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar numero de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular la distancia entre puntos
                    h = (b - a) / n;

                    // Calcular la segunda derivada en puntos extremos
                    fpp[0] = (-f(a+3*h) + 4*f(a+2*h) - 5*f(a+h) + 2*f(a))/(h*h); // This is O(h²)
                    fpp[n] = (2*f(b) - 5*f(b-h) + 4*f(b-2*h) - f(b-3*h))/(h*h); // This is O(h²)
    
                    x[0] = a;
                    x[n] = b;
                
                    // Para los puntos interiores
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fpp[i] = (-f(x[i]+2*h) + 16*f(x[i]+h) - 30*f(x[i]) + 16*f(x[i]-h) - f(x[i]-2*h))/(12*h*h); // This is O(h⁴)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf''(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fpp[i]);
                    }
                
                    // Guardar x[i] y fpp[i] en un archivo de texto
                    save_in_txt(x, fpp, n);
    
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Arreglos para el cálculo de coeficientes polinomiales para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Número de puntos
                    int n;

                    // Leer puntos de datos desde archivo
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("Error al leer datos desde el archivo. Saliendo.\n");
                        return 1;
                    }
                    // Imprime los puntos de datos
                    print_data_points(X, Y, n);

                    // 2. Inicializar A y b
                    /* for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)].
                    // First 2(n-1) ecuaciones: Cada spline pasa por sus dos puntos finales
                    for(int k = 0; k < n-1; k++) {
                        // El spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // El spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
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
                        b[row] = 0.0;
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
                        b[row] = 0.0;
                    }

                    // Dos condiciones de contorno: Spline natural (segunda derivada = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primer spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (último spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Usamos la función de gauss.h para resolver el sistema con eliminación gaussiana
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Dividir el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1); 
                    

                    // Calcular X[n] y Y[n] igualmente espaciados
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }

                    // Inicializar el arreglo x con los nuevos puntos igualmente espaciados (consistente con new_Y)
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    } */

                    // Utilice el código anterior si los datos originales no tienen un espaciado uniforme.
                    h = X[1] - X[0]; // Asumiendo un espaciado uniforme en los datos originales

                    // Calcular la segunda derivada con un error de O(h²) utilizando datos spline.
                    // Diferencia hacia adelante O(h²) para el primer punto (necesita 4 puntos)
                    fpp[0] = (-Y[3] + 4*Y[2] - 5*Y[1] + 2*Y[0])/(h*h);
                    
                    // Diferencia hacia atrás O(h²) para el último punto (necesita 4 puntos)
                    fpp[n-1] = (2*Y[n-1] - 5*Y[n-2] + 4*Y[n-3] - Y[n-4])/(h*h);
    
                    // Diferencia central para los puntos internos
                    for(int i = 1; i <= n-2; i++) {
                        // if(i >= 2 && i <= n-3) {
                            // Diferencia central de 5 puntos O(h⁴) cuando es posible
                        fpp[i] = (-Y[i+2] + 16*Y[i+1] - 30*Y[i] + 16*Y[i-1] - Y[i-2])/(12*h*h);
                        // } else {
                            // Diferencia central de 3 puntos O(h²) cuando no hay suficientes puntos
                            // fpp[i] = (Y[i+1] - 2*Y[i] + Y[i-1])/(h*h);
                        // }
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf''(x)\n");
                    for(int i = 0; i < n; i++) {
                        printf("%lf\t%lf\n", X[i], fpp[i]);
                    }
                
                    // Guardar x[i] y fpp[i] en un archivo de texto
                    save_in_txt(X, fpp, n-1);
    
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados.
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                }
            case 3:
                if(function_data_table == 1) {
                    printf("Insertar extremos del intervalo (a,b)\n");
                    scanf("%lf %lf", &a, &b);
                                
                    printf("Insertar numero de subintervalos\n");
                    scanf("%d", &n);
                                
                    // Calcular la distancia entre puntos
                    h = (b - a) / n;

                    // Calcular la tercera derivada en los puntos extremos
                    fppp[0] = (-3*f(a+4*h) + 14*f(a+3*h) - 24*f(a+2*h) + 18*f(a+h) - 5*f(a))/(2*h*h*h); // This is O(h²)
                    fppp[n] = (5*f(b) - 18*f(b-h) + 24*f(b-2*h) - 14*f(b-3*h) + 3*f(b-4*h))/(2*h*h*h); // This is O(h²)
    
                    x[0] = a;
                    x[n] = b;
                
                    // Para los puntos internos
                    for(int i = 1; i <= n-1; i++) {
                        x[i] = a + i*h;
                        fppp[i] = (-f(x[i]+3*h) + 8*f(x[i]+2*h) - 13*f(x[i]+h) + 13*f(x[i]-h) - 8*f(x[i]-2*h) + f(x[i]-3*h))/(8*h*h*h);// This is O(h⁴)
                    }
    
                    // Imprimir resultados
                    printf("x\t\tf'''(x)\n");
                    for(int i = 0; i <= n; i++) {
                        printf("%lf\t%lf\n", x[i], fppp[i]);
                    }
    
                    // Guardar x[i] y fppp[i] en un archivo de texto
                    save_in_txt(x, fppp, n);
    
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados.
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                } else if(function_data_table == 2) {
                    // Arreglos para el cálculo de coeficientes polinomiales para usar en Spline
                    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                    // Puntos de datos para leer desde archivo de texto
                    double X[MAX_POINTS], Y[MAX_POINTS];
                    // Puntos de datos para calcular la integral
                    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                    // Número de puntos
                    int n;

                    // Leer puntos de datos desde archivo
                    if (!read_data_points("data.txt", X, Y, &n)) {
                        printf("No se pudo leer los datos del archivo. Saliendo.\n");
                        return 1;
                    }
                    // Imprimir los puntos de datos
                    print_data_points(X, Y, n);

                    // 2. Inicializar A y b
                    /* for (int i = 0; i < 4*(n-1); i++) {
                        b[i] = 0.0;
                        for (int j = 0; j < 4*(n-1); j++) {
                            A[i][j] = 0.0;
                        }
                    }

                    // Calculamos A[4(n-1)][4(n-1)] y b[4(n-1)].
                    // Primeras 2(n-1) ecuaciones: Cada spline pasa por sus dos puntos extremos.
                    for(int k = 0; k < n-1; k++) {
                        // La spline k pasa por el punto (X[k], Y[k])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k][4*k+j] = pow(X[k], 3-j);
                        }
                        b[2*k] = Y[k];
                    
                        // La spline k pasa por el punto (X[k+1], Y[k+1])
                        for(int j = 0; j <= 3; j++) {
                            A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                        }
                        b[2*k+1] = Y[k+1];
                    }


                    // Siguientes n-2 ecuaciones: Continuidad de las primeras derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + k;
                        // Primera derivada de la spline k en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                        }
                        // Primera derivada de la spline k+1 en X[k+1]
                        for(int j = 0; j <= 2; j++) {
                            A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                        }
                        b[row] = 0.0;
                    }

                    // Siguientes n-2 ecuaciones: Continuidad de las segundas derivadas
                    for(int k = 0; k < n-2; k++) {
                        int row = 2*(n-1) + (n-2) + k;
                        // Segunda derivada de la spline k en X[k+1]
                        A[row][4*k] = 6 * X[k+1];
                        A[row][4*k+1] = 2;
                        // Segunda derivada de la spline k+1 en X[k+1]
                        A[row][4*(k+1)] = -6 * X[k+1];
                        A[row][4*(k+1)+1] = -2;
                        b[row] = 0.0;
                    }

                    // Dos condiciones de contorno: Spline natural (segunda derivada = 0 en los extremos)
                    int row1 = 4*(n-1) - 2;
                    int row2 = 4*(n-1) - 1;

                    // Segunda derivada = 0 en X[0] (primera spline)
                    A[row1][0] = 6 * X[0];
                    A[row1][1] = 2;
                    b[row1] = 0.0;

                    // Segunda derivada = 0 en X[n-1] (última spline)
                    A[row2][4*(n-2)] = 6 * X[n-1];
                    A[row2][4*(n-2)+1] = 2;
                    b[row2] = 0.0;

                    // Usamos la función de gauss.h para resolver el sistema con eliminación gaussiana
                    gauss_elimination(4*(n-1), A, b, solution);

                    // Divide el intervalo [X[0], X[n-1]] en n-1 subintervalos de igual longitud
                    h = (X[n - 1] - X[0]) / (n - 1); 
                    

                    // Calcula X[n] e Y[n] equidistantes.
                    for(int i = 0; i < n; i++) {
                        new_X[i] = X[0] + i * h;
                        new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                    }
                    
                    // Inicializar el array x con los nuevos puntos equidistantes (de acuerdo con new_Y)
                    for(int i = 0; i < n; i++) {
                        x[i] = new_X[i];
                    }
                    */ 
                   
                    // Usar el código anterior si los datos originales no tienen un espaciado uniforme
                    h = X[1] - X[0]; // Se asume un espaciado uniforme en los datos originales

                    // Verificar el número mínimo de puntos necesarios para la tercera derivada O(h²)
                    if(n < 7) {
                        printf("Error: Se necesitan al menos 7 puntos para la tercera derivada con precisión O(h²).\n");
                        return 1;
                    }
                    
                    // Diferencia hacia adelante O(h²) para el primer punto (necesita 5 puntos)
                    fppp[0] = (-3*Y[4] + 14*Y[3] - 24*Y[2] + 18*Y[1] - 5*Y[0])/(2*h*h*h);
                    
                    // Diferencia hacia atrás O(h²) para el último punto (necesita 5 puntos)
                    fppp[n-1] = (5*Y[n-1] - 18*Y[n-2] + 24*Y[n-3] - 14*Y[n-4] + 3*Y[n-5])/(2*h*h*h);
                
                    // Diferencia central O(h⁴) para puntos interiores
                    for(int i = 1; i <= n-2; i++) {
                        fppp[i] = (-Y[i+3] + 8*Y[i+2] - 13*Y[i+1] + 13*Y[i-1] - 8*Y[i-2] + Y[i-3])/(8*h*h*h);
                    }
    
                    // Imprimo resultados
                    printf("x\t\tf'''(x)\n");
                    for(int i = 0; i < n; i++) {
                        printf("%lf\t%lf\n", X[i], fppp[i]);
                    }
    
                    // Guarda x[i] y fppp[i] en un archivo de texto.
                    save_in_txt(X, fppp, n-1);
    
                    // Finalmente, imprimimos el archivo results.txt en un gráfico usando Python para visualizar los resultados.
                    system("C:\\Users\\sofim\\AppData\\Local\\Programs\\Python\\Python313\\python.exe ..\\graph_points.py");
                    break;
                    
                }
        }
    }

    return 0;
}

double f(double t) {
    return 10 * exp(-t/10) * sin(2*t);
}

void save_in_txt(double x[], double fp[], int n) {
    FILE *archivo = fopen("results.txt", "w");
    if (archivo == NULL) {
        printf("Error: No se puede crear el archivo.\n");
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", x[i], fp[i]);
    }

    fclose(archivo);
}

int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", filename);
        return 0;
    }
    
    printf("El archivo '%s' se abrió correctamente.\n", filename);
    
    // Read number of data points
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: No se puede leer el numero de puntos de datos\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Numero de puntos no valido (%d)\n", *n);
        fclose(fp);
        return 0;
    }
    
    // Read data points
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se puede leer el punto de datos %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Se leyeron correctamente %d puntos de datos\n\n", *n);
    return 1;
}

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

/**
 * Función para evaluar el spline cúbico en un punto x dado
 * @param X Array de valores X (puntos de datos)
 * @param solution Array de coeficientes del spline
 * @param n Número de puntos de datos
 * @param x El valor de x a evaluar
 * @return El valor de y interpolado
 */
double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Encontrar el intervalo de spline correcto
    if (x <= X[0]) {
        k = 0;  // Usar el primer spline para x menor que el primer punto
    } else if (x >= X[n-1]) {
        k = n - 2; // Usar el último spline para x mayor que el último punto
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

double X(double tita) {
    // x(tita) = R*cos(tita) + (L² - R² * sin²(tita))^(1/2)
    // Donde R = 90 mm = 0.09 m y L = 2.5*R = 0.225 m
    double R = 0.09;
    double L = 0.225;
    return (R * cos(tita)) + (sqrt(pow(L, 2) - pow(R, 2) * pow(sin(tita), 2)));
}

void points_generator_deg(const char *file_name) {
    FILE *fp = fopen(file_name, "w");
    if (!fp) {
        perror("Error al abrir el archivo");
        return;
    }

    double tita_deg, tita_rad, x;
    // double h = 2 * PI / (N - 1);

    for (int i = 0; i <= 36; i++) {
        tita_deg = i * 5; // De 0 a 180 grados en pasos de 5 grados
        tita_rad = tita_deg * (PI / 180.0); // Convertir a radianes
        x = X(tita_rad);
        fprintf(fp, "%.6f %.6f\n", tita_deg, x);
    }

    fclose(fp);
    // Imprime el mensaje: "Archivo '%s' generado correctamente."
    printf("Archivo '%s' generado correctamente.\n", file_name);
}

void points_generator_rad(const char *file_name) {
    FILE *fp = fopen(file_name, "w");
    if (!fp) {
        // Imprime el error: "Error al abrir el archivo"
        perror("Error al abrir el archivo");
        return;
    }

    double tita_deg, tita_rad, x;
    // double h = 2 * PI / (N - 1);

    for (int i = 0; i <= 36; i++) {
        // De 0 a 180 grados en pasos de 5 grados
        tita_deg = i * 5; 
        tita_rad = tita_deg * (PI / 180.0); // Convertir a radianes

        x = X(tita_rad);
        fprintf(fp, "%.6f %.6f\n", tita_rad, x); 
    }

    fclose(fp);
    // Imprime el mensaje: "Archivo '%s' generado correctamente."
    printf("Archivo '%s' generado correctamente.\n", file_name);
}
