#include <cstdio>
#include <stdlib.h>
#include <cmath>

/* Este archivo contiene correcciones de errores realizadas con Github Copilot en los métodos de Gauss-Seidel y Relajación.
   Además, se ha corregido la función de impresión de resultados para mostrar números positivos cuando el resultado es -0.0
*/

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
bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_SIZE+1], double b[], int* n);

/**
 * Verifica si la matriz es diagonalmente dominante y que no haya ceros en la diagonal.
 * Retorna 0 si todo está bien, de lo contrario retorna 1
 * @param a Matriz de coeficientes
 * @param n Tamaño de la matriz
 * @return 0 si todo está bien, 1 si hay un cero en la diagonal
 * o una advertencia si la matriz no es diagonalmente dominante
 */
int verificarDominanciaDiagonal(double a[][MAX_SIZE+1], int n);

/** Función para inicializar el vector de estimación con ceros
 * @param Xv Vector a inicializar == X Anterior (Iteración Previa)
 * @param n Tamaño del vector
 * * */ 
void inicializarEstimacion(double Xv[], int n);

/** Función para calcular el error entre Xn y Xv
 * @param Xn X Nuevo = Solución de la matriz
 * @param Xv X Anterior = Iteración previa de la solución
 * @param n Tamaño de los vectores
 * @return El error calculado (Norma Euclidiana)
 * */
double calcularError(double Xn[], double Xv[], int n);


/** Función para imprimir la solución
 * @param nombreMetodo Nombre del método utilizado
 * @param Xn Vector solución
 * @param n Tamaño del vector
 * @param iteraciones Número de iteraciones tomadas para converger
 * @param error Error final
 * */ 
void imprimirSolucion(const char* nombreMetodo, double Xn[], int n, int iteraciones, double error);


/**
 * Implementación del método de Jacobi para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la iteración actual
 * @param n Tamaño de los vectores
 * */
void metodoJacobi(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/**
 * Implementación del método de Gauss-Seidel para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodoGaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);


/**
 * Implementación del método de Relajación (SOR) para resolver sistemas lineales
 * Es similar a Gauss-Seidel pero con un factor de relajación omega y una línea de código adicional
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodoRelajacion(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/* Matriz a resolver (Ejemplo):
    3x1 - 2x2 + x3 + 0x4 + 0x5 + 0x6 = 10
    -2x1 + 4x2 - 2x3 + x4 + 0x5 + 0x6 = -8
    x1 - 2x2 + 4x3 - 2x4 + x5 + 0x6 = 10
    0x1 + 1x2 - 2x3 + 4x4 - 2x5 + x6 = 10
    0x1 + 0x2 + x3 - 2x4 + 4x5 - 2x6 = -8
    0x1 + 0x2 + 0x3 + 1x4 - 2x5 + 3x6 = 10

    Encontramos el siguiente sistema de ecuaciones a resolver:

    En otros términos: 
    A = | 3 -2 1 0 0 0|
        | -2 4 -2 1 0 0|
        | 1 -2 4 -2 1 0|
        | 0 1 -2 4 -2 1|
        | 0 0 1 -2 4 -2|
        | 0 0 0 1 -2 3|


    b = | 10 |
        | -8 |
        | 10 |
        | 10 |
        | -8 |
        | 10 |

    Tolerancia = 1e-11
    

    Aplicaremos el Método de Relajación para resolver el sistema de ecuaciones.

    --> omega = 1.1

    ⚠️ Advertencia: La matriz no es diagonalmente dominante en la fila 2.
    ⚠️ Advertencia: La matriz no es diagonalmente dominante en la fila 3.
    ⚠️ Advertencia: La matriz no es diagonalmente dominante en la fila 4.
    ⚠️ Advertencia: La matriz no es diagonalmente dominante en la fila 5.

    La solución es:
        
    x1 = 2
    x2 = 0
    x3 = 4
    x4 = 4
    x5 = 0
    x6 = 2

    El método convergió en 39 iteraciones con un error de 0.0

    Volveremos a aplicar el mismo método pero usando omega = 1.0 para resolver el sistema de ecuaciones anterior.

    --> omega = 1.0

    La solución es:
    
    x1 = 2
    x2 = 0
    x3 = 4
    x4 = 4
    x5 = 0
    x6 = 2

    El método convergió en 43 iteraciones con un error de 0.0

    IMPORTANTE: 
    --> Observa que la solución con omega > 1 (sobre-relajación) converge más rápido que con omega = 1 (Gauss-Seidel).
    --> Observa que la solución al sistema lineal se encuentra aun cuando no es diagonalmente dominante. Esto demuestra que la condición de dominancia diagonal
        es solo una condición suficiente, no necesaria para que exista la convergencia.

    si omega = 1 tenemos Gauss-Seidel
    si omega 0 <= omega < 1 tenemos sub-relajación --> Se utiliza para hacer converger un sistema que de otra manera no convergería usando Gauss-Seidel 
    si omega 1 < omega <= 2 tenemos sobre-relajación --> Se utiliza para acelerar la convergencia de un sistema que ya converge con Gauss-Seidel
    */

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma, auxiliar;
    double Xv[MAX_SIZE+1], Xn[MAX_SIZE+1]; // X Anterior, X Nuevo
    double tolerancia, error_anterior, error_nuevo;
    int iteraciones;
    double omega;

    // Definiendo arreglos usando la constante global MAX_SIZE
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Leer matriz desde el archivo usando la función
    if(!leer_archivo_matriz("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia4\\Ejercicio5\\data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Sistema de ecuaciones original:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);        // 10 espacios totales, 6 decimales
    }
    printf("\n");

    // Verificación de diagonal dominante
    if (verificarDominanciaDiagonal(a, n) != 0) {
        printf("El método no puede continuar. La matriz tiene ceros en la diagonal. El programa termina.");
        return 1;
    }
    printf("Verificación completada.\n");

    printf("Elija un método para resolver el sistema:\n");
    printf("1. Método de Jacobi\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Método de Relajación\n");
    printf("Ingrese opción: ");

    int opcion;
    scanf("%d", &opcion);

    switch(opcion) {
        case 1:
            metodoJacobi(a, b, Xv, Xn, n);
            break;
        case 2:
            metodoGaussSeidel(a, b, Xv, Xn, n);
            break;
        case 3:
            metodoRelajacion(a, b, Xv, Xn, n);
            break;
        default:
            printf("❌ Opción inválida.\n");
    }


    return 0;
}

int verificarDominanciaDiagonal(double a[][MAX_SIZE+1], int n) {
    for (int i = 1; i <= n; i++) {
        double suma = 0.0;

        // Primero verificamos si hay un cero en la diagonal
        if (fabs(a[i][i]) == 0.0) {
            printf("❌ Error: Elemento cero en la diagonal en la posición a[%d][%d].\n", i, i);
            return 1;
        }

        // Suma de elementos fuera de la diagonal
        for (int j = 1; j <= n; j++) {
            if (j != i) {
                suma += fabs(a[i][j]);
            }
        }

        // Verificamos la condición de dominancia
        if (fabs(a[i][i]) < suma) {
            printf("⚠️ Advertencia: La matriz no es diagonalmente dominante en la fila %d.\n", i);
        }
    }

    return 0; // Todo está OK
}

void inicializarEstimacion(double Xv[], int n) {
    for(int i = 1; i <= n; i++) {
        Xv[i] = 0.0;
    }
}

double calcularError(double Xn[], double Xv[], int n) {
    double error = 0.0;
    for(int i = 1; i <= n; i++) {
        error += pow(Xn[i] - Xv[i], 2);
    }
    return sqrt(error);
}

void imprimirSolucion(const char* nombreMetodo, double Xn[], int n, int iteraciones, double error) {
    printf("------------------SOLUCIÓN DEL %s------------------\n", nombreMetodo);
    printf("La solución del sistema es:\n");
    for(int i = 1; i <= n; i++) {
        // Corrección para el problema de visualización de -0.0: si el valor es muy cercano a cero, mostrar como 0.0
        double valor = (fabs(Xn[i]) < 1e-10) ? 0.0 : Xn[i];
        printf("Xn[%d] = %10.6lf\n", i, valor);
    }
    printf("El método convergió en %d iteraciones con un error de %10.6lf\n", iteraciones, error);
}


void metodoJacobi(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double suma, tolerancia, error_anterior, error_nuevo;
    int iteraciones;

    inicializarEstimacion(Xv, n);

    printf("Ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    error_anterior = 1000;
    iteraciones = 0;

    do {
        iteraciones++;
        for(int i = 1; i <= n; i++) {
            suma = 0.0;
            for(int j = 1; j <= n; j++) {
                if(j != i) {
                    suma += a[i][j] * Xv[j];
                }
            }
            Xn[i] = (b[i] - suma) / a[i][i];
        }

        error_nuevo = calcularError(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimirSolucion("MÉTODO DE JACOBI", Xn, n, iteraciones, error_nuevo);
}

void metodoGaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double suma, tolerancia, error_anterior, error_nuevo;
    int iteraciones;

    inicializarEstimacion(Xv, n);

    printf("Ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    error_anterior = 1000;
    iteraciones = 0;

    do {
        iteraciones++;
        for(int i = 1; i <= n; i++) {
            suma = 0.0;  // Reiniciar la suma para cada fila
            
            // Suma de elementos antes de la diagonal (usando valores NUEVOS Xn)
            for(int j = 1; j <= i-1; j++) {
                suma += a[i][j] * Xn[j];
            }
            
            // Suma de elementos después de la diagonal (usando valores ANTERIORES Xv)
            for(int j = i+1; j <= n; j++) {
                suma += a[i][j] * Xv[j];
            }
            
            Xn[i] = (b[i] - suma) / a[i][i];
        }

        error_nuevo = calcularError(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimirSolucion("GAUSS-SEIDEL", Xn, n, iteraciones, error_nuevo);
}

void metodoRelajacion(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double suma, tolerancia, error_anterior, error_nuevo, omega;
    int iteraciones;

    inicializarEstimacion(Xv, n);

    printf("Ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    printf("Ingrese el factor de relajación (0 < omega < 2):");
    scanf("%lf", &omega);

    error_anterior = 1000;
    iteraciones = 0;

    do {
        iteraciones++;
        for(int i = 1; i <= n; i++) {
            suma = 0.0;  // Reiniciar la suma para cada fila
            
            // Suma de elementos antes de la diagonal (usando valores NUEVOS Xn)
            for(int j = 1; j <= i-1; j++) {
                suma += a[i][j] * Xn[j];
            }
            
            // Suma de elementos después de la diagonal (usando valores ANTERIORES Xv)
            for(int j = i+1; j <= n; j++) {
                suma += a[i][j] * Xv[j];
            }
            
            // Calcular el paso de Gauss-Seidel
            double gauss_seidel = (b[i] - suma) / a[i][i];
            
            // Aplicar factor de relajación (SOR)
            Xn[i] = omega * gauss_seidel + (1.0 - omega) * Xv[i];
        }

        error_nuevo = calcularError(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimirSolucion("MÉTODO DE RELAJACIÓN", Xn, n, iteraciones, error_nuevo);
}


bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_SIZE+1], double b[], int* n) {
    FILE *p_archivo;
    char c;
    
    // Abrir archivo de datos
    p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }
    
    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar filas en el archivo (Lógica para determinar n)
    int filas = 0;
    while(!feof(p_archivo)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), p_archivo) != NULL) {
            filas++;
        }
    }
    
    // El sistema debe ser cuadrado, entonces n = filas
    *n = filas;
    printf("📊 Tamaño del sistema: %d x %d\n", *n, *n);

    // Cerrar y reabrir el archivo para reiniciar el puntero
    fclose(p_archivo);
    p_archivo = fopen(nombre_archivo, "r");
    
    // Verificar el tamaño máximo usando la constante global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_SIZE);
        fclose(p_archivo);
        return false;
    }

    // Leer matriz aumentada desde el archivo
    // Formato esperado: cada fila contiene n coeficientes + 1 término independiente
    // Ejemplo para 3x3: a11 a12 a13 b1
    //                   a21 a22 a23 b2  
    //                   a31 a32 a33 b3
    
    int i, j;
    for(i = 1; i <= *n; i++) {
        // Leyendo los coeficientes de la matriz
        for(j = 1; j <= *n; j++) {
            if(fscanf(p_archivo, "%lf", &a[i][j]) != 1) {
                printf("❌ Error al leer el elemento a[%d][%d]\n", i, j);
                fclose(p_archivo);
                return false;
            }
        }
        // Leer el término independiente
        if(fscanf(p_archivo, "%lf", &b[i]) != 1) {
            printf("❌ Error al leer el término independiente b[%d]\n", i);
            fclose(p_archivo);
            return false;
        }
    }
    
    fclose(p_archivo);
    printf("✅ Matriz leída exitosamente del archivo\n\n");
    
    return true;
}