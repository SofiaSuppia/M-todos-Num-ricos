#include <cstdio>
#include <stdlib.h>
#include <cmath>

/* Como se menciona, X**(k+1) - X**(k) solo mide el cambio entre iteraciones, mientras que la norma
de A*X**(k) - b mide qué tan bien la aproximación satisface el sistema original. Esto lo hace más confiable. */

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
 * Verifica si la matriz es estrictamente diagonalmente dominante y que no haya ceros en la diagonal.
 * Retorna 0 si todo está bien, de lo contrario retorna 1
 * @param a Matriz de coeficientes
 * @param n Tamaño de la matriz
 * @return 0 si todo está bien, 1 si hay un cero en la diagonal
 * o una advertencia si la matriz no es diagonalmente dominante
 */
int verificarDominanciaDiagonal(double a[][MAX_SIZE+1], int n);

/** Función para inicializar el vector de estimación con ceros
 * @param Xv Vector a inicializar == X Antiguo
 * @param n Tamaño del vector
 * * */ 
void inicializarEstimacion(double Xv[], int n);

/** Función para calcular el residual ||AX - b||
 * Es otra forma de medir el error de la solución
 * Generalmente es un criterio más apropiado al decidir cuándo detener el proceso iterativo del método.
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * @return El residual calculado (norma euclidiana)
 */
double calcularResidual(double a[][MAX_SIZE+1], double b[], double Xn[], int n);

/** Función para imprimir la solución
 * @param nombreMetodo Nombre del método utilizado
 * @param Xn Vector solución
 * @param n Tamaño del vector
 * @param iteraciones Número de iteraciones requeridas para converger
 * @param error Error final (Norma del residual)
 * */ 
void imprimirSolucion(const char* nombreMetodo, double Xn[], int n, int iteraciones, double error);


/**
 * Implementación del método de Jacobi para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Antiguo = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodoJacobi(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/**
 * Implementación del método de Gauss-Seidel para resolver sistemas lineales
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Antiguo = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodoGaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);


/**
 * Implementación del método de Relajación (SOR) para resolver sistemas lineales
 * Es similar a Gauss-Seidel pero con un factor de relajación omega y una línea de código adicional
 * @param a Matriz de coeficientes
 * @param b Vector de términos independientes
 * @param Xv X Antiguo = Iteración previa de la solución
 * @param Xn X Nuevo = Solución de la matriz
 * @param n Tamaño de los vectores
 * */
void metodoRelajacion(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/* Inciso a)
Matriz a resolver:
    Encontramos el siguiente sistema de ecuaciones a resolver:

    En otros términos: 
    A = | 3 1 1 |
        | 2 6 1 |
        | 1 1 4 |

    b = | 5 |
        | 9 |
        | 6 |

    Tolerancia = 1e-11
    

    Aplicaremos Jacobi para resolver el sistema de ecuaciones.

    La solución es:
        x1 = 1
        x2 = 1
        x3 = 1

    El método convergió en 48 iteraciones con un error de 0.0

    Aplicaremos Gauss-Seidel para resolver el sistema de ecuaciones anterior.

    La solución es:
        x1 = 1
        x2 = 1
        x3 = 1

    El método convergió en 14 iteraciones con un error de   0.0

Inciso b)
Matriz a resolver:
Encontramos el siguiente sistema de ecuaciones a resolver:

    En otros términos: 
    A = | 5 7 6 5  |
        | 7 10 8 7 |
        | 5 7 9 10 |

    b = | 23 |
        | 32 |
        | 33 |
        | 31 |

    Tolerancia = 1e-4

    ⚠️  Advertencia: La matriz no es diagonalmente dominante en la fila 1.
    ⚠️  Advertencia: La matriz no es diagonalmente dominante en la fila 2.
    ⚠️  Advertencia: La matriz no es diagonalmente dominante en la fila 3.
    ⚠️  Advertencia: La matriz no es diagonalmente dominante en la fila 4.
    

    Aplicaremos Jacobi para resolver el sistema de ecuaciones. 
    El método no converge, detenemos el proceso.


    Aplicaremos Gauss-Seidel para resolver el sistema de ecuaciones anterior.
    El método no converge, detenemos el proceso.

*/


int main(int argc, char const *argv[]) {
    int n, p;
    double factor, producto, suma, aux;
    double Xv[MAX_SIZE+1], Xn[MAX_SIZE+1]; // X Antiguo, X Nuevo
    double tolerancia, error_anterior, error_nuevo;
    int iteraciones;
    double omega;

    // Definiendo arreglos usando el tamaño máximo global
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Leer arreglo desde archivo usando la función
    if(!leer_archivo_matriz("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia4\\Ejercicio3\\data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Sistema de ecuaciones original:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 espacios totales, 6 decimales
        }
        printf("| %10.6lf\n", b[i]);       // 10 espacios totales, 6 decimales
    }
    printf("\n");

    // Verificación de diagonal dominante
    if (verificarDominanciaDiagonal(a, n) != 0) {
        printf("El método no puede continuar. La matriz tiene ceros en la diagonal. El programa termina.\n");
        return 1;
    }
    printf("Verificación completada.\n");

    printf("Elige un método para resolver el sistema:\n");
    printf("1. Método de Jacobi\n");
    printf("2. Método de Gauss-Seidel\n");
    printf("3. Método de Relajación (SOR)\n");
    printf("Introduce la opción: ");

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

        // Suma de los elementos fuera de la diagonal
        for (int j = 1; j <= n; j++) {
            if (j != i) {
                suma += fabs(a[i][j]);
            }
        }

        // Verificamos la condición de dominancia
        if (fabs(a[i][i]) < suma) {
            printf("⚠️  Advertencia: La matriz no es diagonalmente dominante en la fila %d.\n", i);
        }
    }

    return 0; // Todo está OK
}

void inicializarEstimacion(double Xv[], int n) {
    for(int i = 1; i <= n; i++) {
        Xv[i] = 0.0;
    }
}

/* double calcularError(double Xn[], double Xv[], int n) {
    double error = 0.0;
    for(int i = 1; i <= n; i++) {
        error += pow(Xn[i] - Xv[i], 2);
    }
    return sqrt(error);
} */

double calcularResidual(double a[][MAX_SIZE+1], double b[], double Xn[], int n) {
    double residual[MAX_SIZE+1];
    for (int i = 1; i <= n; i++) {
        double suma_producto = 0.0;
        for (int j = 1; j <= n; j++) {
            suma_producto += a[i][j] * Xn[j];
        }
        residual[i] = suma_producto - b[i]; 
    }
    // Norma euclidiana del residual
    double norma = 0.0;
    for (int i = 1; i <= n; i++) {
        norma += residual[i] * residual[i];
    }
    return sqrt(norma);
}


void imprimirSolucion(const char* nombreMetodo, double Xn[], int n, int iteraciones, double error) {
    printf("------------------SOLUCIÓN DEL %s------------------\n", nombreMetodo);
    printf("La solución del sistema es:\n");
    for(int i = 1; i <= n; i++) {
        printf("Xn[%d] = %10.6lf\n", i, Xn[i]);
    }
    printf("El método convergió en %d iteraciones con una norma del residual de %10.6lf\n", iteraciones, error);
}


void metodoJacobi(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double suma, tolerancia, error_anterior, error_nuevo;
    int iteraciones;

    inicializarEstimacion(Xv, n);

    printf("Por favor, introduce la tolerancia:");
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

        // error_nuevo = calcularError(Xn, Xv, n);
        error_nuevo = calcularResidual(a, b, Xn, n); // Usando el residual como error

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

    printf("Por favor, introduce la tolerancia:");
    scanf("%lf", &tolerancia);

    error_anterior = 1000;
    iteraciones = 0;

    do {
        iteraciones++;
        for(int i = 1; i <= n; i++) {
            suma = 0.0;  // Reiniciar suma para cada fila
            
            // Sumar elementos antes de la diagonal (usando NUEVOS valores Xn)
            for(int j = 1; j <= i-1; j++) {
                suma += a[i][j] * Xn[j];
            }
            
            // Sumar elementos después de la diagonal (usando ANTIGUOS valores Xv)
            for(int j = i+1; j <= n; j++) {
                suma += a[i][j] * Xv[j];
            }
            
            Xn[i] = (b[i] - suma) / a[i][i];
        }

        // error_nuevo = calcularError(Xn, Xv, n);
        error_nuevo = calcularResidual(a, b, Xn, n); // Usando el residual como error

        if(error_nuevo > error_anterior) {
            printf("El método no converge, detenemos el proceso.\n");
            return;
        }

        error_anterior = error_nuevo;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(error_nuevo > tolerancia);

    imprimirSolucion("MÉTODO DE GAUSS-SEIDEL", Xn, n, iteraciones, error_nuevo);
}

void metodoRelajacion(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double suma, tolerancia, error_anterior, error_nuevo, omega;
    int iteraciones;

    inicializarEstimacion(Xv, n);

    printf("Por favor, introduce la tolerancia:");
    scanf("%lf", &tolerancia);

    printf("Por favor, introduce el factor de relajación (0 < omega < 2):");
    scanf("%lf", &omega);

    error_anterior = 1000;
    iteraciones = 0;

    do {
        iteraciones++;
        for(int i = 1; i <= n; i++) {
            suma = 0.0;  // Reiniciar suma para cada fila
            
            // Sumar elementos antes de la diagonal (usando NUEVOS valores Xn)
            for(int j = 1; j <= i-1; j++) {
                suma += a[i][j] * Xn[j];
            }
            
            // Sumar elementos después de la diagonal (usando ANTIGUOS valores Xv)
            for(int j = i+1; j <= n; j++) {
                suma += a[i][j] * Xv[j];
            }
            
            // Calcular el paso de Gauss-Seidel
            double gauss_seidel = (b[i] - suma) / a[i][i];
            
            // Aplicar el factor de relajación (SOR)
            Xn[i] = omega * gauss_seidel + (1.0 - omega) * Xv[i];
        }

        // error_nuevo = calcularError(Xn, Xv, n);
        error_nuevo = calcularResidual(a, b, Xn, n); // Usando el residual como error

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
    FILE *fp;
    char c;
    
    // Abrir archivo de datos
    fp = fopen(nombre_archivo, "r");
    if (fp == NULL) {
        printf("❌ Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return false;
    }
    
    printf("✅ Archivo '%s' abierto\n\n", nombre_archivo);

    // Contar filas en el archivo
    /* int filas = 0;
    while((c = fgetc(fp)) != EOF) {
        if(c == '\n') {
            filas++;
        }
    } */

    // Esta alternativa funciona mejor que la anterior
    int filas = 0;
    while(!feof(fp)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), fp) != NULL) {
            filas++;
        }
    }
    
    // El sistema debe ser cuadrado, por lo tanto n = filas
    *n = filas;
    printf("📊 Tamaño del sistema: %d x %d\n", *n, *n);

    // Cerrar y volver a abrir el archivo para reiniciar el puntero
    fclose(fp);
    fp = fopen(nombre_archivo, "r");
    
    // Verificar el tamaño máximo usando la constante global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: Sistema demasiado grande (%d). Máximo permitido: %d\n", *n, MAX_SIZE);
        fclose(fp);
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
            if(fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("❌ Error al leer el elemento a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        // Leer el término independiente
        if(fscanf(fp, "%lf", &b[i]) != 1) {
            printf("❌ Error al leer el término independiente b[%d]\n", i);
            fclose(fp);
            return false;
        }
    }
    
    fclose(fp);
    printf("✅ Arreglo leído correctamente desde el archivo\n\n");
    
    return true;
}