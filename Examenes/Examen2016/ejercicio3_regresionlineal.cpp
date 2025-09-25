#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50
#define MAX_TAMANO 100 // Usado para la matriz de las Ecuaciones Normales
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

/* a = 0.999996 = 1
    b = -1.999974 = -2
    f(x) = exp(x²) - 2
*/

int main(int argc, char const *argv[]) {
    // Arreglos para los puntos de datos a leer del archivo de texto
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para el cálculo de los coeficientes polinomiales (o de ajuste)
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion_sistema[MAX_TAMANO+1];
    // Solución (coeficientes) del polinomio de Mínimos Cuadrados
    double a[MAX_TAMANO+1];
    // n = número de puntos de datos y grado es el grado del polinomio
    int n, grado;
    // Sumas para los cálculos de mínimos cuadrados
    double suma_xy, suma_x, suma_y, distancia, suma_e, suma_ye;
    double Sr, St, r, promedio_y, f_xi;

    printf("Ingrese el grado del polinomio (grado >= 1): ");
    scanf("%d", &grado); // Esta variable 'grado' se usa de forma limitada en el código actual (se asume ajuste lineal/no lineal de 2 parámetros)

    // Leer puntos de datos del archivo
    if (!leer_puntos_datos("exercise3.txt", X, Y, &n)) {
        printf("Fallo al leer los datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    // Verificar si tenemos suficientes puntos de datos
    if(n < grado) {
        printf("Error: No hay suficientes puntos de datos para el grado del polinomio elegido.\n");
        return 1;
    }

    // Construir la matriz A y el vector b para las ecuaciones normales
    // Nota: El siguiente bucle parece estar codificado específicamente para un ajuste
    // de 2 parámetros (grado + 1 = 2), ignorando la variable 'grado' para la matriz.
    for(int l = 0; l <= 1; l++) {
        // Calcular los términos independientes (b[l])
        suma_ye = 0.0;
        for(int i = 0; i < n; i++) {
            // Esto implica Y_i * g_l(X_i) donde g_l(X_i) = exp((1-l) * X_i^2)
            suma_ye += Y[i] * exp((1 - l) * pow(X[i], 2));
        }
        b[l] = suma_ye;
        
        // Calcular la Matriz A[l][m]
        for(int m = 0; m <= 1; m++) {
            suma_e = 0.0;
            for(int i = 0; i < n; i++) {
                // Esto implica g_l(X_i) * g_m(X_i) donde g_k(X_i) = exp((1-k) * X_i^2)
                suma_e += exp((2-l-m) * pow(X[i], 2));
            }
            A[l][m] = suma_e;
        }
    }


    // Usamos la función de gauss.h para resolver el sistema con Eliminación Gaussiana
    // Se usa 'grado + 1' como tamaño del sistema (2x2), aunque se haya ignorado 'grado' arriba.
    gauss_elimination(grado + 1, A, b, solucion_sistema);

    // Copiamos la solución a 'a[i]' para darle relevancia a nuestro contexto
    for(int i = 0; i <= grado; i++) {
        a[i] = solucion_sistema[i];
    }

    // Imprimir los coeficientes del polinomio/ajuste
    printf("------------------SOLUCIÓN DEL SISTEMA------------------\n");
    printf("La solución del sistema es:\n");
    for(int i = 0; i <= grado; i++) {
        printf("a[%d] = %lf\n", i, a[i]);
    }
    
    // NOTA: 'imprimir_polinomio' asume un ajuste polinomial estándar (a0 + a1*x + ...)
    printf("\n------------------POLINOMIO INTERPOLADOR/DE AJUSTE------------------\n");
    imprimir_polinomio(a, n);

    // Calcular Sr, St, r 
    suma_y = 0.0;
    St = 0.0;
    Sr = 0.0;

    // Calcular el promedio de Y
    for(int i = 0; i < n; i++) {
        suma_y += Y[i];
    }
    promedio_y = suma_y / n;

    // Calcular St - suma total de cuadrados (error total)
    for(int i = 0; i < n; i++) {
        St += pow((Y[i] - promedio_y), 2); // El orden (avg - Y) o (Y - avg) no afecta al ser elevado al cuadrado
    }

    // Calcular Sr - suma de cuadrados de los residuos (error del modelo)
    for(int i = 0; i < n; i++) {
        f_xi = 0.0;
        // La siguiente línea asume que el ajuste es un polinomio estándar:
        // P(x) = a0 + a1*x + a2*x^2 + ...
        for(int k = 0; k <= grado; k++) {
            f_xi += a[k] * pow(X[i], k);
        }
        Sr += pow((Y[i] - f_xi), 2); // El orden (f_xi - Y_i) o (Y_i - f_xi) no afecta
    }
    
    // Calcular r - coeficiente de correlación
    r = sqrt((St - Sr) / St);
    
    // Imprimir coeficiente de correlación 
    printf("\nEl coeficiente de correlación r es: %lf\n", r);

    // Verificar la bondad del ajuste
    distancia = fabs(r - 1.0);

    // Interpretación de la bondad del ajuste
    if(distancia < 0.1) {
        printf("El ajuste es muy bueno (r está cerca de 1)\n");
    } else if(distancia < 0.25) {
        printf("El ajuste es bueno\n");
    } else if(distancia < 0.5) {
        printf("El ajuste es aceptable\n");
    } else {
        printf("El ajuste es pobre\n");
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
    
    // Leer puntos de datos
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
    // Nota: Esta función imprime la forma polinomial estándar P(x) = a0 + a1*x + ...
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