#include <cstdio>  // Incluye las funciones de I/O de C (printf, scanf, FILE)
#include <stdlib.h>
#include <cmath>   // Incluye las funciones matemáticas (fabs, pow, sqrt)

/* Este archivo fue creado utilizando el conocimiento de mi profesor de la universidad */

// Ahora puedes manejar matrices de hasta 50x50
#define MAX_TAMANO 50  // Tamaño máximo de la matriz y vectores

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
 * @return true si el archivo fue leído exitosamente, false en caso contrario
 * */
bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_TAMANO+1], double b[], int* n);

/**
 * Comprueba si la matriz es diagonalmente dominante y que no haya ceros en la diagonal.
 * Devuelve 0 si todo está bien, de lo contrario devuelve 1
 * @param a Matriz de coeficientes
 * @param n Tamaño de la matriz
 * @return 0 si todo está bien, 1 si hay un cero en la diagonal
 * o una advertencia si la matriz no es diagonalmente dominante
 */
int verificar_dominancia_diagonal(double a[][MAX_TAMANO+1], int n);

/** * Función para inicializar el vector de estimación (guess) con ceros
 * @param Xv Vector a inicializar == X Anterior
 * @param n Tamaño del vector
 * */
void inicializar_estimacion(double Xv[], int n);

/** * Función para calcular el error entre Xn y Xv
 * @param Xn X Nuevo = Solución de la matriz
 * @param Xv X Anterior = Solución de la iteración previa
 * @param n Tamaño de los vectores
 * @return El error calculado (norma euclidiana)
 * */
double calcular_error(double Xn[], double Xv[], int n);


/** * Función para imprimir la solución
 * @param nombre_metodo Nombre del método utilizado
 * @param Xn Vector solución
 * @param n Tamaño del vector
 * @param iteraciones Número de iteraciones tomadas para converger
 * @param error Error final
 * */
void imprimir_solucion(const char* nombre_metodo, double Xn[], int n, int iteraciones, double error);


/**
 * Implementación del método de Jacobi para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Solución de la iteración previa
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodo_jacobi(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n);

/**
 * Implementación del método de Gauss-Seidel para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Solución de la iteración previa
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodo_gauss_seidel(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n);


/**
 * Implementación del método de Relajación (SOR) para resolver sistemas lineales
 * Es similar a Gauss-Seidel pero con un factor de relajación omega y una línea de código adicional
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Anterior = Solución de la iteración previa
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodo_relajacion(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n);

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma_temp, aux;
    double Xv[MAX_TAMANO+1], Xn[MAX_TAMANO+1]; // X Anterior, X Nuevo
    double tolerancia, error_anterior, error_nuevo;
    int iteraciones;
    double omega;

    // Definición de arreglos usando el MAX_TAMANO global
    double a[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], X[MAX_TAMANO+1];

    // Leer arreglo desde el archivo usando la función
    if(!leer_archivo_matriz("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia4\\data.dat", a, b, &n)) {
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

    // Verificación de diagonal dominante
    if (verificar_dominancia_diagonal(a, n) != 0) {
        printf("No se puede continuar con el método. La matriz tiene ceros en la diagonal o no es dominante. Se cierra el programa\n");
        return 1;
    }
    printf("Verificación completada.\n");

    printf("Elija un método para resolver el sistema:\n");
    printf("1. Método de Jacobi\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Método de Relajación (SOR)\n");
    printf("Ingrese opción: ");

    int opcion;
    scanf("%d", &opcion);

    switch(opcion) {
        case 1:
            metodo_jacobi(a, b, Xv, Xn, n);
            break;
        case 2:
            metodo_gauss_seidel(a, b, Xv, Xn, n);
            break;
        case 3:
            metodo_relajacion(a, b, Xv, Xn, n);
            break;
        default:
            printf("❌ Opción inválida.\n");
    }


    return 0;
}

int verificar_dominancia_diagonal(double a[][MAX_TAMANO+1], int n) {
    for (int i = 1; i <= n; i++) {
        double suma_fila = 0.0;

        // Primero verificamos si hay cero en la diagonal
        if (fabs(a[i][i]) == 0.0) {
            printf("❌ Error: Elemento cero en la diagonal en la posición a[%d][%d].\n", i, i);
            return 1;
        }

        // Suma de los elementos fuera de la diagonal
        for (int j = 1; j <= n; j++) {
            if (j != i) {
                suma_fila += fabs(a[i][j]);
            }
        }

        // Verificamos la condición de dominancia
        if (fabs(a[i][i]) < suma_fila) {
            printf("⚠️  Advertencia: La matriz no es estrictamente diagonalmente dominante en la fila %d.\n", i);
        }
    }

    return 0; // Todo está OK
}

void inicializar_estimacion(double Xv[], int n) {
    for(int i = 1; i <= n; i++) {
        Xv[i] = 0.0;
    }
}

double calcular_error(double Xn[], double Xv[], int n) {
    double error = 0.0;
    for(int i = 1; i <= n; i++) {
        error += pow(Xn[i] - Xv[i], 2);
    }
    return sqrt(error); // Norma euclidiana
}

void imprimir_solucion(const char* nombre_metodo, double Xn[], int n, int iteraciones, double error) {
    printf("------------------SOLUCIÓN DEL %s------------------\n", nombre_metodo);
    printf("La solución del sistema es:\n");
    for(int i = 1; i <= n; i++) {
        printf("X[%d] = %10.6lf\n", i, Xn[i]);
    }
    printf("El método convergió en %d iteraciones con un error de %10.6lf\n", iteraciones, error);
}


void metodo_jacobi(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n) {
    double suma_temp, tolerancia, error_anterior, error_nuevo;
    int iteraciones;

    inicializar_estimacion(Xv, n);

    printf("Por favor, ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    error_anterior = 1000.0;
    iteraciones = 0;

    do {
        iteraciones++;
        // Cálculo de las nuevas variables Xn usando SOLO los valores de Xv
        for(int i = 1; i <= n; i++) {
            suma_temp = 0.0;
            for(int j = 1; j <= n; j++) {
                if(j != i) {
                    suma_temp += a[i][j] * Xv[j]; // Usa X anterior (Xv)
                }
            }
            Xn[i] = (b[i] - suma_temp) / a[i][i];
        }

        error_nuevo = calcular_error(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        // Actualizar Xv para la siguiente iteración (Xv = Xn)
        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimir_solucion("MÉTODO DE JACOBI", Xn, n, iteraciones, error_nuevo);
}

void metodo_gauss_seidel(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n) {
    double suma_temp, tolerancia, error_anterior, error_nuevo;
    int iteraciones;

    inicializar_estimacion(Xv, n);
    // Inicializar Xn con la estimación inicial también
    for(int i = 1; i <= n; i++) {
        Xn[i] = Xv[i];
    }

    printf("Por favor, ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    error_anterior = 1000.0;
    iteraciones = 0;

    do {
        iteraciones++;
        // Cálculo de las nuevas variables Xn
        for(int i = 1; i <= n; i++) {
            suma_temp = 0.0;
            
            // Suma con los valores de Xn ya calculados (j < i)
            for(int j = 1; j <= i-1; j++) {
                suma_temp += a[i][j] * Xn[j];
            }
            
            // Suma con los valores de Xv (X anterior) para los que aún no se han calculado (j > i)
            for(int j = i+1; j <= n; j++) {
                suma_temp += a[i][j] * Xv[j];
            }
            
            // El cálculo solo se almacena en Xn
            Xn[i] = (b[i] - suma_temp) / a[i][i];
        }

        error_nuevo = calcular_error(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        // Actualizar Xv para la siguiente iteración (Xv = Xn)
        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimir_solucion("GAUSS-SEIDEL", Xn, n, iteraciones, error_nuevo);
}

void metodo_relajacion(double a[][MAX_TAMANO+1], double b[], double Xv[], double Xn[], int n) {
    double suma_temp, tolerancia, error_anterior, error_nuevo, omega;
    int iteraciones;

    inicializar_estimacion(Xv, n);
    // Inicializar Xn con la estimación inicial también
    for(int i = 1; i <= n; i++) {
        Xn[i] = Xv[i];
    }

    printf("Por favor, ingrese la tolerancia:");
    scanf("%lf", &tolerancia);

    printf("Por favor, ingrese el factor de relajación (0 < omega < 2):");
    scanf("%lf", &omega);

    error_anterior = 1000.0;
    iteraciones = 0;

    do {
        iteraciones++;
        // Cálculo de las nuevas variables Xn
        for(int i = 1; i <= n; i++) {
            suma_temp = 0.0;
            
            // Suma con los valores de Xn ya calculados (j < i)
            for(int j = 1; j <= i-1; j++) {
                suma_temp += a[i][j] * Xn[j];
            }
            
            // Suma con los valores de Xv (X anterior) para los que aún no se han calculado (j > i)
            for(int j = i+1; j <= n; j++) {
                suma_temp += a[i][j] * Xv[j];
            }
            
            // Se calcula el valor "temporal" (similar a Gauss-Seidel)
            double X_gs_temp = (b[i] - suma_temp) / a[i][i];

            // Paso de relajación (SOR):
            // Xn[i] = omega * X_gauss_seidel + (1 - omega) * X_anterior
            Xn[i] = (omega * X_gs_temp) + ((1.0 - omega) * Xv[i]);
        }

        error_nuevo = calcular_error(Xn, Xv, n);

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        // Actualizar Xv para la siguiente iteración (Xv = Xn)
        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimir_solucion("MÉTODO DE RELAJACIÓN (SOR)", Xn, n, iteraciones, error_nuevo);
}


bool leer_archivo_matriz(const char* nombre_archivo, double a[][MAX_TAMANO+1], double b[], int* n) {
    FILE *p_archivo;
    char c;
    
    // Abrir archivo de datos
    p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }
    
    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar filas en el archivo (mejor alternativa)
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
    
    // Verificar el tamaño máximo usando MAX_TAMANO global
    if(*n > MAX_TAMANO) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_TAMANO);
        fclose(p_archivo);
        return false;
    }

    // Leer la matriz aumentada del archivo
    // Formato esperado: cada fila contiene n coeficientes + 1 término independiente
    
    int i, j;
    for(i = 1; i <= *n; i++) {
        // Leer los coeficientes de la matriz
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