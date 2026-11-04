#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define MAX_SIZE 200

#include "gauss.h"

/*
 * ============================================================================
 * SPLINES LINEALES Y SPLINES CUBICOS NATURALES
 * ============================================================================
 *
 * MAPA DE CALCULOS DE ESTE ARCHIVO
 * --------------------------------
 * 1) Lectura y validacion:
 *    - Se leen los puntos (X[i], Y[i]) desde "data.txt".
 *    - Se exige que haya al menos 2 puntos y que las X sean estrictamente
 *      crecientes. Si no, el programa termina con un mensaje de error.
 *
 * 2) Spline lineal (opcion 1):
 *    - Para cada intervalo [X[k], X[k+1]] se calcula la pendiente:
 *          m_k = (Y[k+1] - Y[k]) / (X[k+1] - X[k])
 *    - Esa pendiente es la derivada constante del spline lineal.
 *    - La evaluacion usa:
 *          y = Y[k] + m_k * (x - X[k])
 *
 * 3) Spline cubico natural (opcion 2):
 *    - Hay num_puntos puntos y num_intervalos = num_puntos - 1.
 *    - Hay 4 * num_intervalos incognitas:
 *          a_k, b_k, c_k, d_k  para k = 0, ..., num_intervalos - 1
 *    - El polinomio de cada intervalo es:
 *          f_k(x) = a_k*x^3 + b_k*x^2 + c_k*x + d_k
 *    - El sistema A*Z = b se arma con cuatro bloques:
 *        1. Interpolacion: cada polinomio pasa por los extremos de su intervalo.
 *        2. Continuidad de la primera derivada en los nodos interiores.
 *        3. Continuidad de la segunda derivada en los nodos interiores.
 *        4. Frontera natural: f''(x_0) = 0 y f''(x_ultimo) = 0.
 *
 * 4) DERIVADAS:
 *    - Se usan derivadas simbolicas del polinomio cubico para armar el sistema:
 *          f_k'(x)  = 3*a_k*x^2 + 2*b_k*x + c_k
 *          f_k''(x) = 6*a_k*x + 2*b_k
 *    - NO se calculan derivadas numericas ni se imprimen f'(x) o f''(x).
 *    - La unica "derivada" que se imprime como dato es la pendiente del
 *      spline lineal en la opcion 1.
 *
 * 5) SUMATORIAS:
 *    - En ESTE archivo no hay acumuladores del tipo "suma += ...".
 *      Solo hay bucles for que recorren puntos, intervalos y filas de la matriz.
 *    - La unica sumatoria real esta en gauss.h, en la sustitucion hacia atras.
 *
 * 6) INTEGRALES:
 *    - En este archivo NO se calculan integrales.
 *
 * 7) Resolucion del sistema:
 *    - Se llama a gauss_elimination() de gauss.h.
 *    - El vector solucion contiene los coeficientes:
 *          Z[4*k]     = a_k
 *          Z[4*k + 1] = b_k
 *          Z[4*k + 2] = c_k
 *          Z[4*k + 3] = d_k
 *
 * 8) Evaluacion:
 *    - Se busca el intervalo que contiene a x_gorro.
 *    - Si x_gorro esta fuera del rango [X[0], X[num_puntos-1]], NO se extrapola:
 *      se muestra un mensaje y no se evalua.
 *    - Si esta dentro, se evalua el polinomio correspondiente.
 * ============================================================================
 */

/**
 * Lee los puntos desde un archivo y valida los datos.
 *
 * Formato del archivo:
 *   primera linea: numero de puntos
 *   siguientes lineas: X[i] Y[i]
 *
 * Validaciones:
 *   - numero de puntos entre 2 y MAX_POINTS.
 *   - el sistema cubico debe caber en MAX_SIZE.
 *   - las X deben ser estrictamente crecientes.
 *
 * @return 1 si todo salio bien, 0 si hubo error.
 */
int leer_puntos_de_datos(const char *nombre_archivo, double X[], double Y[],
                         int *num_puntos);

/**
 * Muestra en pantalla los puntos leidos.
 * No realiza calculos, solo impresion.
 */
void imprimir_puntos_de_datos(double X[], double Y[], int num_puntos);

/**
 * Arma la matriz A y el vector b del spline cubico natural.
 *
 * Esta funcion contiene los cuatro bloques del sistema:
 *   1. interpolacion,
 *   2. continuidad de la primera derivada,
 *   3. continuidad de la segunda derivada,
 *   4. condiciones naturales.
 *
 * @return el tamano del sistema (4 * num_intervalos) o -1 si hay error.
 */
int construir_sistema_spline_cubico(double X[], double Y[], int num_puntos,
                                    double A[][MAX_SIZE + 1], double b[]);

/**
 * Imprime los coeficientes a_k, b_k, c_k, d_k de cada spline cubico.
 */
void imprimir_splines_cubicos(double X[], double solucion[], int num_puntos);

/**
 * Busca el intervalo k tal que X[k] <= x <= X[k+1].
 *
 * @return k si x esta dentro del rango, -1 si esta fuera.
 */
int encontrar_intervalo(double X[], int num_puntos, double x);

/**
 * Evalua el spline cubico en x.
 *
 * @param y Puntero donde se guarda el resultado.
 * @return 1 si x esta dentro del rango y se evaluo, 0 si esta fuera.
 */
int evaluar_spline(double X[], double solucion[], int num_puntos,
                   double x, double *y);

/**
 * Calcula e imprime la pendiente de cada tramo del spline lineal.
 *
 * La pendiente m_k es la derivada del spline lineal en ese intervalo.
 * No se calculan integrales.
 */
void spline_lineal(double X[], double Y[], int num_puntos);

int main(void) {
    /* Puntos de datos leidos desde el archivo. */
    double X[MAX_POINTS], Y[MAX_POINTS];

    /* Matriz A, vector b y vector solucion para el sistema cubico. */
    double A[MAX_SIZE + 1][MAX_SIZE + 1];
    double b[MAX_SIZE + 1];
    double solucion[MAX_SIZE + 1];

    int num_puntos;
    int opcion;

    /* ------------------------------------------------------------------
     * PASO 1: leer los puntos y validar el archivo.
     * ------------------------------------------------------------------ */
    if (!leer_puntos_de_datos("data.txt", X, Y, &num_puntos)) {
        printf("Fallo al leer los datos del archivo. Saliendo.\n");
        return 1;
    }

    imprimir_puntos_de_datos(X, Y, num_puntos);

    /* ------------------------------------------------------------------
     * PASO 2: elegir el tipo de interpolacion.
     * ------------------------------------------------------------------ */
    printf("Seleccione una opcion para la interpolacion por splines:\n");
    printf("1. Spline Lineal\n");
    printf("2. Spline Cubico\n");
    printf("Opcion: ");

    if (scanf("%d", &opcion) != 1) {
        printf("Entrada invalida.\n");
        return 1;
    }

    switch (opcion) {
        case 1:
            /* ----------------------------------------------------------
             * OPCION 1: SPLINE LINEAL
             *
             * Aqui se calcula una pendiente por intervalo y se evalua
             * la recta que une dos puntos consecutivos.
             * No se usa la matriz A ni el metodo de Gauss.
             * ---------------------------------------------------------- */
            spline_lineal(X, Y, num_puntos);

            printf("¿Desea evaluar un punto x_gorro? (1 para si, 0 para no): ");
            if (scanf("%d", &opcion) != 1) {
                break;
            }

            while (opcion != 0) {
                double x_gorro;
                double y_gorro;
                int k;

                printf("Ingrese el valor de x_gorro: ");
                if (scanf("%lf", &x_gorro) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                /* Localiza el intervalo que contiene a x_gorro. */
                k = encontrar_intervalo(X, num_puntos, x_gorro);

                if (k < 0) {
                    /* No se extrapola fuera del rango de datos. */
                    printf("x_gorro = %.6f esta fuera del rango [%.6f, %.6f].\n",
                           x_gorro, X[0], X[num_puntos - 1]);
                    printf("No se realiza extrapolacion.\n");
                } else {
                    /* Formula del spline lineal:
                     *   y = Y[k] + m_k * (x_gorro - X[k])
                     * donde m_k = (Y[k+1] - Y[k]) / (X[k+1] - X[k]).
                     */
                    y_gorro = Y[k] + (Y[k + 1] - Y[k]) *
                              (x_gorro - X[k]) / (X[k + 1] - X[k]);

                    printf("El valor interpolado en x_gorro = %.6f es y_gorro = %.6f\n",
                           x_gorro, y_gorro);
                }

                printf("¿Desea evaluar otro punto x_gorro? (1 para si, 0 para no): ");
                if (scanf("%d", &opcion) != 1) {
                    break;
                }
            }
            break;

        case 2: {
            /* ----------------------------------------------------------
             * OPCION 2: SPLINE CUBICO NATURAL
             *
             * Se arma el sistema A*Z = b, se resuelve con eliminacion
             * gaussiana y se evalua el polinomio del intervalo elegido.
             * ---------------------------------------------------------- */
            int tam = construir_sistema_spline_cubico(X, Y, num_puntos, A, b);

            if (tam < 0) {
                printf("Error: no se pudo construir el sistema del spline cubico.\n");
                break;
            }

            /* Resuelve el sistema. La solucion queda en "solucion". */
            gauss_elimination(tam, A, b, solucion);

            /* Muestra los coeficientes de cada polinomio cubico. */
            imprimir_splines_cubicos(X, solucion, num_puntos);

            printf("¿Desea evaluar un punto x_gorro? (1 para si, 0 para no): ");
            if (scanf("%d", &opcion) != 1) {
                break;
            }

            while (opcion != 0) {
                double x_gorro;
                double y_gorro;

                printf("Ingrese el valor de x_gorro: ");
                if (scanf("%lf", &x_gorro) != 1) {
                    printf("Entrada invalida.\n");
                    break;
                }

                if (evaluar_spline(X, solucion, num_puntos, x_gorro, &y_gorro)) {
                    printf("El valor interpolado en x_gorro = %.6f es y_gorro = %.6f\n",
                           x_gorro, y_gorro);
                } else {
                    printf("x_gorro = %.6f esta fuera del rango [%.6f, %.6f].\n",
                           x_gorro, X[0], X[num_puntos - 1]);
                    printf("No se realiza extrapolacion.\n");
                }

                printf("¿Desea evaluar otro punto x_gorro? (1 para si, 0 para no): ");
                if (scanf("%d", &opcion) != 1) {
                    break;
                }
            }
            break;
        }

        default:
            printf("Opcion invalida. Intente de nuevo.\n");
            break;
    }

    return 0;
}

int leer_puntos_de_datos(const char *nombre_archivo, double X[], double Y[],
                         int *num_puntos) {
    FILE *fp;

    /* Abre el archivo en modo lectura. */
    fp = fopen(nombre_archivo, "r");
    if (fp == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return 0;
    }

    printf("Archivo '%s' abierto exitosamente\n", nombre_archivo);

    /* Lee la cantidad de puntos. */
    if (fscanf(fp, "%d", num_puntos) != 1) {
        printf("Error: No se pudo leer el numero de puntos de datos\n");
        fclose(fp);
        return 0;
    }

    /* Valida que haya al menos 2 puntos y no mas de MAX_POINTS. */
    if (*num_puntos < 2 || *num_puntos > MAX_POINTS) {
        printf("Error: se requieren entre 2 y %d puntos (se leyeron %d)\n",
               MAX_POINTS, *num_puntos);
        fclose(fp);
        return 0;
    }

    /* Valida que el sistema cubico 4*(num_puntos-1) quepa en la matriz. */
    if (4 * (*num_puntos - 1) > MAX_SIZE) {
        printf("Error: el sistema necesita %d incognitas y MAX_SIZE es %d\n",
               4 * (*num_puntos - 1), MAX_SIZE);
        fclose(fp);
        return 0;
    }

    /* Lee los pares (X[i], Y[i]). */
    for (int i = 0; i < *num_puntos; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el punto de datos %d\n", i + 1);
            fclose(fp);
            return 0;
        }

        /* Las X deben ser estrictamente crecientes para definir intervalos. */
        if (i > 0 && X[i] <= X[i - 1]) {
            printf("Error: las X deben ser estrictamente crecientes.\n");
            printf("X[%d] = %.6f no es mayor que X[%d] = %.6f\n",
                   i, X[i], i - 1, X[i - 1]);
            fclose(fp);
            return 0;
        }
    }

    fclose(fp);
    printf("Se leyeron %d puntos de datos exitosamente\n\n", *num_puntos);
    return 1;
}

void imprimir_puntos_de_datos(double X[], double Y[], int num_puntos) {
    printf("Puntos de Datos:\n");
    printf("=============\n");
    printf("    i   |      Xi       |      Yi       \n");
    printf("------|--------------|-------------\n");

    for (int i = 0; i < num_puntos; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }

    printf("\n");
}

int construir_sistema_spline_cubico(double X[], double Y[], int num_puntos,
                                    double A[][MAX_SIZE + 1], double b[]) {
    int num_intervalos = num_puntos - 1;
    int tam = 4 * num_intervalos;

    if (num_puntos < 2 || tam > MAX_SIZE) {
        return -1;
    }

    /* Inicializa toda la matriz A y el vector b en cero.
     * Esto es necesario porque solo se escriben algunos coeficientes.
     */
    for (int i = 0; i < tam; i++) {
        b[i] = 0.0;
        for (int j = 0; j < tam; j++) {
            A[i][j] = 0.0;
        }
    }

    /* ------------------------------------------------------------------
     * BLOQUE 1: INTERPOLACION
     *
     * Para cada intervalo k se impone:
     *   f_k(X[k])   = Y[k]
     *   f_k(X[k+1]) = Y[k+1]
     *
     * Como f_k(x) = a_k*x^3 + b_k*x^2 + c_k*x + d_k, los coeficientes
     * se guardan en las columnas 4*k ... 4*k+3.
     * ------------------------------------------------------------------ */
    for (int k = 0; k < num_intervalos; k++) {
        for (int j = 0; j <= 3; j++) {
            /* j=0 -> x^3, j=1 -> x^2, j=2 -> x, j=3 -> 1 */
            A[2 * k][4 * k + j] = pow(X[k], 3 - j);
            A[2 * k + 1][4 * k + j] = pow(X[k + 1], 3 - j);
        }

        b[2 * k] = Y[k];
        b[2 * k + 1] = Y[k + 1];
    }

    /* ------------------------------------------------------------------
     * BLOQUE 2: CONTINUIDAD DE LA PRIMERA DERIVADA
     *
     * En cada nodo interior X[k+1] se impone:
     *   f_k'(X[k+1]) = f_{k+1}'(X[k+1])
     *
     * Con f'(x) = 3*a*x^2 + 2*b*x + c, los coeficientes son:
     *   j=0 -> 3*x^2, j=1 -> 2*x, j=2 -> 1
     *
     * Hay num_intervalos - 1 nodos interiores, por lo tanto hay
     * num_intervalos - 1 ecuaciones de este tipo.
     * ------------------------------------------------------------------ */
    for (int k = 0; k < num_intervalos - 1; k++) {
        int fila = 2 * num_intervalos + k;

        for (int j = 0; j <= 2; j++) {
            double termino = (3 - j) * pow(X[k + 1], 2 - j);

            /* Coeficientes del spline k. */
            A[fila][4 * k + j] = termino;

            /* Coeficientes del spline k+1, con signo negativo. */
            A[fila][4 * (k + 1) + j] = -termino;
        }

        b[fila] = 0.0;
    }

    /* ------------------------------------------------------------------
     * BLOQUE 3: CONTINUIDAD DE LA SEGUNDA DERIVADA
     *
     * En cada nodo interior X[k+1] se impone:
     *   f_k''(X[k+1]) = f_{k+1}''(X[k+1])
     *
     * Con f''(x) = 6*a*x + 2*b:
     *   coeficiente de a -> 6*x
     *   coeficiente de b -> 2
     *
     * Los coeficientes c y d no aparecen porque su derivada segunda es cero.
     * ------------------------------------------------------------------ */
    for (int k = 0; k < num_intervalos - 1; k++) {
        int fila = 3 * num_intervalos - 1 + k;

        /* Segunda derivada del spline k en X[k+1]. */
        A[fila][4 * k] = 6 * X[k + 1];
        A[fila][4 * k + 1] = 2;

        /* Segunda derivada del spline k+1 en X[k+1], con signo negativo. */
        A[fila][4 * (k + 1)] = -6 * X[k + 1];
        A[fila][4 * (k + 1) + 1] = -2;

        b[fila] = 0.0;
    }

    /* ------------------------------------------------------------------
     * BLOQUE 4: CONDICIONES NATURALES
     *
     * Se impone segunda derivada cero en los extremos:
     *   f_0''(X[0]) = 0
     *   f_{n-1}''(X[n]) = 0
     *
     * Dividiendo por 2, la ecuacion queda:
     *   3*a_0*X[0] + b_0 = 0
     *   3*a_ultimo*X[ultimo] + b_ultimo = 0
     * ------------------------------------------------------------------ */
    {
        int fila_izquierda = tam - 2;
        int fila_derecha = tam - 1;
        int ultimo_spline = 4 * (num_intervalos - 1);

        /* Extremo izquierdo: coeficientes a_0 y b_0. */
        A[fila_izquierda][0] = 3 * X[0];
        A[fila_izquierda][1] = 1;
        b[fila_izquierda] = 0.0;

        /* Extremo derecho: coeficientes a_ultimo y b_ultimo. */
        A[fila_derecha][ultimo_spline] = 3 * X[num_puntos - 1];
        A[fila_derecha][ultimo_spline + 1] = 1;
        b[fila_derecha] = 0.0;
    }

    return tam;
}

void imprimir_splines_cubicos(double X[], double solucion[], int num_puntos) {
    printf("------------------SOLUCION------------------\n");
    printf("Coeficientes del Spline Cubico:\n");

    for (int k = 0; k < num_puntos - 1; k++) {
        printf("Spline %d (desde X[%d]=%.3f hasta X[%d]=%.3f):\n",
               k + 1, k, X[k], k + 1, X[k + 1]);

        /* Coeficientes: a = solucion[4k], b = 4k+1, c = 4k+2, d = 4k+3. */
        printf("  S%d(x) = %.6fx³ + %.6fx² + %.6fx + %.6f\n",
               k + 1,
               solucion[4 * k],
               solucion[4 * k + 1],
               solucion[4 * k + 2],
               solucion[4 * k + 3]);
    }

    printf("\n");
}

int encontrar_intervalo(double X[], int num_puntos, double x) {
    /* Si x esta fuera del rango de datos, no se busca intervalo. */
    if (x < X[0] || x > X[num_puntos - 1]) {
        return -1;
    }

    /* Busca k tal que X[k] <= x <= X[k+1]. */
    for (int k = 0; k < num_puntos - 1; k++) {
        if (x >= X[k] && x <= X[k + 1]) {
            return k;
        }
    }

    return -1;
}

int evaluar_spline(double X[], double solucion[], int num_puntos,
                   double x, double *y) {
    int k = encontrar_intervalo(X, num_puntos, x);

    /* Si x esta fuera del rango, no se extrapola. */
    if (k < 0) {
        return 0;
    }

    /* Evalua el polinomio cubico del intervalo k:
     *   y = a*x^3 + b*x^2 + c*x + d
     */
    *y = solucion[4 * k] * pow(x, 3) +
         solucion[4 * k + 1] * pow(x, 2) +
         solucion[4 * k + 2] * x +
         solucion[4 * k + 3];

    return 1;
}

void spline_lineal(double X[], double Y[], int num_puntos) {
    printf("------------------SPLINES LINEALES------------------\n");

    for (int k = 0; k < num_puntos - 1; k++) {
        /* Pendiente del tramo k. Es la derivada del spline lineal. */
        double pendiente = (Y[k + 1] - Y[k]) / (X[k + 1] - X[k]);

        printf("Spline %d (desde X[%d]=%.3f hasta X[%d]=%.3f):\n",
               k + 1, k, X[k], k + 1, X[k + 1]);

        /* Ecuacion de la recta: y = Y[k] + m_k * (x - X[k]). */
        printf("  f_%d(X) = %.6f + %.6f * (X - %.6f)\n\n",
               k + 1, Y[k], pendiente, X[k]);
    }
}
