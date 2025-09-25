#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50 // Máximo número de puntos
#define MAX_TAMANO 100 // Máximo tamaño para la matriz de las ecuaciones normales
#include "gauss.h" // Incluye la cabecera para la eliminación gaussiana

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
 * Función para imprimir el polinomio de ajuste Pn(x)
 * @param a Arreglo de coeficientes (a0, a1, a2, ...)
 * @param n Grado del polinomio + 1 (para iterar sobre los coeficientes a[0] hasta a[grado])
 */
void imprimir_polinomio(double a[], int n);

int main(int argc, char const *argv[]) {
    // Arreglos para los puntos de datos a leer del archivo de texto
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    // Arreglos para el cálculo de coeficientes: Matriz A, vector b y solución
    double A[MAX_TAMANO+1][MAX_TAMANO+1], b[MAX_TAMANO+1], solucion[MAX_TAMANO+1];
    // Solución del polinomio de Mínimos Cuadrados (coeficientes a0, a1, ...)
    double a[MAX_TAMANO+1];
    // n = número de puntos de datos, y grado es el grado del polinomio de ajuste
    int n, grado;
    // Sumatorias para los cálculos de mínimos cuadrados
    double suma_xy, suma_x, suma_y, distancia;
    // Sr: Suma de los cuadrados de los residuos, St: Suma total de cuadrados, r: Coeficiente de correlación
    double Sr, St, r, promedio_y, f_xi;

    printf("Ingrese el grado del polinomio de ajuste (grado >= 1): ");
    scanf("%d", &grado);

    // Leer puntos de datos desde el archivo
    if (!leer_puntos_datos("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Guias\\Guia5\\data.txt", X, Y, &n)) {
        printf("Error al leer datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Imprimir los puntos de datos
    imprimir_puntos_datos(X, Y, n);

    // Verificar si tenemos suficientes puntos de datos (n debe ser mayor que el grado del polinomio)
    if(n < grado + 1) {
        printf("Error: No hay suficientes puntos de datos (%d) para el grado de polinomio elegido (%d).\n", n, grado);
        return 1;
    }

    // Construir la matriz A y el vector b para el sistema de Ecuaciones Normales
    for(int l = 0; l <= grado; l++) {
        suma_xy = 0.0;
        for(int i = 0; i < n; i++) {
            // b[l] = Sum(Yi * X_i^l)
            suma_xy += Y[i] * pow(X[i], l);
        }
        b[l] = suma_xy;
        for(int m = 0; m <= grado; m++) {
            suma_x = 0.0;
            for(int i = 0; i < n; i++) {
                // A[l][m] = Sum(X_i^(l+m))
                suma_x += pow(X[i], l + m); 
            }
            A[l][m] = suma_x;
        }
    }

    // Usamos la función de gauss.h para resolver el sistema con Eliminación Gaussiana
    // El tamaño del sistema es (grado + 1) x (grado + 1)
    eliminacion_gaussiana(grado + 1, A, b, solucion);

    // Copiamos la solución a a[i] para nuestro contexto de coeficientes (a0, a1, ...)
    for(int i = 0; i <= grado; i++) {
        a[i] = solucion[i];
    }

    // Imprimir los coeficientes del polinomio
    printf("------------------SOLUCIÓN------------------\n");
    printf("La solución del sistema (coeficientes del polinomio) es:\n");
    for(int i = 0; i <= grado; i++) {
        printf("a[%d] = %lf\n", i, a[i]);
    }
    
    printf("\n------------------POLINOMIO DE AJUSTE------------------\n");
    // El segundo argumento es (grado + 1) para indicar cuántos coeficientes hay
    imprimir_polinomio(a, grado + 1);

    // Calcular Sr, St, r 
    suma_y = 0.0;
    St = 0.0;
    Sr = 0.0;

    // Calcular promedio de Y
    for(int i = 0; i < n; i++) {
        suma_y += Y[i];
    }
    promedio_y = suma_y / n;

    // Calcular St - Suma total de cuadrados
    for(int i = 0; i < n; i++) {
        St += pow((Y[i] - promedio_y), 2);
    }

    // Calcular Sr - Suma de los cuadrados de los residuos
    for(int i = 0; i < n; i++) {
        f_xi = 0.0;
        // Evaluar el polinomio de ajuste en X[i]
        for(int k = 0; k <= grado; k++) {
            f_xi += a[k] * pow(X[i], k);
        }
        // Sr = Sumatoria de (f(Xi) - Yi)^2
        Sr += pow((Y[i] - f_xi), 2);
    }
    
    // Calcular r - Coeficiente de correlación
    r = sqrt((St - Sr) / St);
    
    // Imprimir coeficiente de correlación 
    printf("\nEl coeficiente de correlación r es: %lf\n", r);

    // Verificar la bondad del ajuste
    distancia = fabs(r - 1.0);

    // Interpretación de la bondad del ajuste
    if(distancia < 0.05) { // Usando un margen más estricto
        printf("El ajuste es excelente (r es muy cercano a 1)\n");
    } else if(distancia < 0.15) {
        printf("El ajuste es muy bueno\n");
    } else if(distancia < 0.3) {
        printf("El ajuste es aceptable\n");
    } else {
        printf("El ajuste es pobre o inadecuado\n");
    }

    return 0;
}

// Implementación de funciones auxiliares (traducidas)

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
        printf("Error: Número de puntos inválido (%d). Máximo permitido: %d\n", *n, MAX_PUNTOS);
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
    printf("P(x) = ");
    
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
            // Nota: Se elimina la impresión de 1.0 si no es el término constante, pero se mantiene si el término es X^i
            if (fabs(coeficiente - 1.0) > 1e-10 || i == 0) {
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