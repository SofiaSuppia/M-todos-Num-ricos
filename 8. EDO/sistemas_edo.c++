#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * ============================================================================
 * SISTEMAS DE ECUACIONES DIFERENCIALES DE PRIMER ORDEN
 * ============================================================================
 *
 * Este programa resuelve un sistema de N ecuaciones de primer orden:
 *
 *      y1' = f1(x, y1, y2, ..., yN)
 *      y2' = f2(x, y1, y2, ..., yN)
 *      ...
 *      yN' = fN(x, y1, y2, ..., yN)
 *
 * con condiciones iniciales:
 *
 *      y1(x0) = y10,  y2(x0) = y20, ..., yN(x0) = yN0
 *
 * MÉTODOS DISPONIBLES:
 *      1. Euler
 *      2. Heun
 *      3. Punto Medio
 *      4. Runge-Kutta de orden 4 (RK4)
 *
 * IMPORTANTE: no confundir N con n.
 *      N = cantidad de ecuaciones del sistema (orden de la EDO original).
 *      n = cantidad de subintervalos usados por el método numérico.
 *
 * ----------------------------------------------------------------------------
 * CÓMO USAR PARA N = 3 Y N = 4
 * ----------------------------------------------------------------------------
 * Más abajo hay bloques #if N == 2, #elif N == 3 y #elif N == 4.
 * Solamente hay que cambiar la línea:
 *
 *      #define N 2
 *
 * por:
 *
 *      #define N 3     o     #define N 4
 *
 * y el programa compila con el ejemplo correspondiente.
 * Si se quiere resolver otro sistema, se debe editar la función f() del bloque
 * activo y, si se conoce, la función exacta().
 *
 * Ejemplo para N = 3:
 *      Sistema: y1' = y2, y2' = y3, y3' = 6
 *      Equivale a la EDO de tercer orden: y''' = 6
 *      Condiciones iniciales: y1(0) = 0, y2(0) = 0, y3(0) = 0
 *      Solución exacta: y1 = x^3, y2 = 3x^2, y3 = 6x
 *      Datos a ingresar: x0 = 0, xf = 1, n = 10, y1(0)=0, y2(0)=0, y3(0)=0
 *
 * Ejemplo para N = 4:
 *      Sistema: y1' = y2, y2' = y3, y3' = y4, y4' = 24
 *      Equivale a la EDO de cuarto orden: y'''' = 24
 *      Condiciones iniciales: y1(0)=0, y2(0)=0, y3(0)=0, y4(0)=0
 *      Solución exacta: y1 = x^4, y2 = 4x^3, y3 = 12x^2, y4 = 24x
 *      Datos a ingresar: x0 = 0, xf = 1, n = 10, y1(0)=0, y2(0)=0, y3(0)=0, y4(0)=0
 *
 * ============================================================================
 */

#define N 2
#define MAX_PASOS 200

/* ==========================================================================
 * DEFINICIÓN DEL SISTEMA Y DE LA SOLUCIÓN EXACTA
 * ========================================================================== */

#if N == 2

/*
 * Ejemplo de sistema 2x2:
 *      y1' = 3x + y2
 *      y2' = x^2 - y1 - 1
 *
 * Condiciones iniciales:
 *      y1(0) = 0
 *      y2(0) = 1
 *
 * Solución exacta:
 *      y1(x) = x^2 + sin(x)
 *      y2(x) = cos(x) - x
 */
void f(double x, const double y[], double dy[]) {
    dy[0] = 3.0 * x + y[1];
    dy[1] = x * x - y[0] - 1.0;
}

void exacta(double x, double ye[]) {
    ye[0] = x * x + sin(x);
    ye[1] = cos(x) - x;
}

#elif N == 3

/*
 * Ejemplo de sistema 3x3 que representa una EDO de tercer orden:
 *      y''' = 6
 *
 * Se define:
 *      y1 = y
 *      y2 = y'
 *      y3 = y''
 *
 * Entonces:
 *      y1' = y2
 *      y2' = y3
 *      y3' = 6
 *
 * Condiciones iniciales:
 *      y1(0) = 0
 *      y2(0) = 0
 *      y3(0) = 0
 *
 * Solución exacta:
 *      y1 = x^3
 *      y2 = 3x^2
 *      y3 = 6x
 */
void f(double x, const double y[], double dy[]) {
    (void)x;
    dy[0] = y[1];
    dy[1] = y[2];
    dy[2] = 6.0;
}

void exacta(double x, double ye[]) {
    ye[0] = x * x * x;
    ye[1] = 3.0 * x * x;
    ye[2] = 6.0 * x;
}

#elif N == 4

/*
 * Ejemplo de sistema 4x4 que representa una EDO de cuarto orden:
 *      y'''' = 24
 *
 * Se define:
 *      y1 = y
 *      y2 = y'
 *      y3 = y''
 *      y4 = y'''
 *
 * Entonces:
 *      y1' = y2
 *      y2' = y3
 *      y3' = y4
 *      y4' = 24
 *
 * Condiciones iniciales:
 *      y1(0) = 0
 *      y2(0) = 0
 *      y3(0) = 0
 *      y4(0) = 0
 *
 * Solución exacta:
 *      y1 = x^4
 *      y2 = 4x^3
 *      y3 = 12x^2
 *      y4 = 24x
 */
void f(double x, const double y[], double dy[]) {
    (void)x;
    dy[0] = y[1];
    dy[1] = y[2];
    dy[2] = y[3];
    dy[3] = 24.0;
}

void exacta(double x, double ye[]) {
    ye[0] = x * x * x * x;
    ye[1] = 4.0 * x * x * x;
    ye[2] = 12.0 * x * x;
    ye[3] = 24.0 * x;
}

#else
#error "N debe ser 2, 3 o 4. Agregue un bloque #elif N == ... con su f y exacta."
#endif

/* ==========================================================================
 * UN PASO DE CADA MÉTODO
 * ========================================================================== */

/*
 * Calcula un solo paso del método elegido.
 *
 * Entrada:
 *      metodo : 1 = Euler, 2 = Heun, 3 = Punto Medio, 4 = RK4
 *      x      : valor actual de x
 *      y[]    : valores actuales de todas las variables
 *      h      : tamaño de paso
 *
 * Salida:
 *      ynext[] : valores en x + h
 */
void paso_metodo(int metodo, double x, const double y[], double h, double ynext[]) {
    switch (metodo) {
        case 1: {
            /*
             * EULER:
             *      k = f(x, y)
             *      ynext = y + h * k
             *
             * Es una suma de la pendiente actual multiplicada por h.
             */
            double k[N];
            f(x, y, k);
            for (int j = 0; j < N; j++) {
                ynext[j] = y[j] + h * k[j];
            }
            break;
        }

        case 2: {
            /*
             * HEUN (predictor-corrector):
             *      k1 = f(x, y)
             *      yp = y + h * k1
             *      k2 = f(x + h, yp)
             *      ynext = y + (h/2) * (k1 + k2)
             *
             * Se promedian la pendiente inicial y la pendiente final estimada.
             */
            double k1[N], k2[N], yp[N];
            f(x, y, k1);
            for (int j = 0; j < N; j++) {
                yp[j] = y[j] + h * k1[j];
            }
            f(x + h, yp, k2);
            for (int j = 0; j < N; j++) {
                ynext[j] = y[j] + (h / 2.0) * (k1[j] + k2[j]);
            }
            break;
        }

        case 3: {
            /*
             * PUNTO MEDIO:
             *      k1 = f(x, y)
             *      k2 = f(x + h/2, y + (h/2) * k1)
             *      ynext = y + h * k2
             *
             * Se evalúa la pendiente en el punto medio del intervalo.
             */
            double k1[N], k2[N], yp[N];
            f(x, y, k1);
            for (int j = 0; j < N; j++) {
                yp[j] = y[j] + (h / 2.0) * k1[j];
            }
            f(x + h / 2.0, yp, k2);
            for (int j = 0; j < N; j++) {
                ynext[j] = y[j] + h * k2[j];
            }
            break;
        }

        case 4: {
            /*
             * RUNGE-KUTTA DE ORDEN 4:
             *      k1 = f(x, y)
             *      k2 = f(x + h/2, y + (h/2)*k1)
             *      k3 = f(x + h/2, y + (h/2)*k2)
             *      k4 = f(x + h,   y + h*k3)
             *      ynext = y + (h/6) * (k1 + 2k2 + 2k3 + k4)
             *
             * La suma ponderada usa los coeficientes 1/6, 2/6, 2/6 y 1/6.
             */
            double k1[N], k2[N], k3[N], k4[N], yp[N];

            f(x, y, k1);
            for (int j = 0; j < N; j++) {
                yp[j] = y[j] + (h / 2.0) * k1[j];
            }

            f(x + h / 2.0, yp, k2);
            for (int j = 0; j < N; j++) {
                yp[j] = y[j] + (h / 2.0) * k2[j];
            }

            f(x + h / 2.0, yp, k3);
            for (int j = 0; j < N; j++) {
                yp[j] = y[j] + h * k3[j];
            }

            f(x + h, yp, k4);
            for (int j = 0; j < N; j++) {
                ynext[j] = y[j] + (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
            }
            break;
        }

        default: {
            /*
             * Si el método no es válido, se deja el valor anterior.
             * La validación real se hace en main().
             */
            for (int j = 0; j < N; j++) {
                ynext[j] = y[j];
            }
            break;
        }
    }
}

/* ==========================================================================
 * RESOLUCIÓN COMPLETA DEL SISTEMA
 * ========================================================================== */

/*
 * Aplica el método elegido desde x0 hasta x0 + n*h.
 *
 * X[i]    guarda el valor de x en el paso i.
 * Y[i][j] guarda la variable j-ésima en el paso i.
 */
void resolver(int metodo, double x0, double h, int n, const double y0[], double X[], double Y[][N]) {
    X[0] = x0;
    for (int j = 0; j < N; j++) {
        Y[0][j] = y0[j];
    }

    for (int i = 0; i < n; i++) {
        /*
         * Se calcula el paso i+1 a partir del paso i.
         * Luego se avanza x: x_{i+1} = x_i + h.
         */
        paso_metodo(metodo, X[i], Y[i], h, Y[i + 1]);
        X[i + 1] = X[i] + h;
    }
}

/* ==========================================================================
 * NOMBRES Y ÓRDENES DE LOS MÉTODOS
 * ========================================================================== */

const char *nombre_metodo(int metodo) {
    switch (metodo) {
        case 1: return "Euler";
        case 2: return "Heun";
        case 3: return "Punto Medio";
        case 4: return "RK4";
        default: return "Desconocido";
    }
}

const char *archivo_metodo(int metodo) {
    switch (metodo) {
        case 1: return "euler";
        case 2: return "heun";
        case 3: return "punto_medio";
        case 4: return "rk4";
        default: return "desconocido";
    }
}

int orden_metodo(int metodo) {
    switch (metodo) {
        case 1: return 1; // Euler es de orden 1
        case 2: return 2; // Heun es de orden 2
        case 3: return 2; // Punto Medio es de orden 2
        case 4: return 4; // RK4 es de orden 4
        default: return 0;
    }
}

/* ==========================================================================
 * GUARDADO DE RESULTADOS
 * ========================================================================== */

/*
 * Guarda una tabla con:
 *      x, y1, ..., yN, y1_exacta, ..., yN_exacta, error1, ..., errorN
 *
 * El error se calcula como:
 *      error_j = |y_j_exacta - y_j_numerica|
 */
void guardar_resultados(const char *nombre, const double X[], const double Y[][N], int n) {
    FILE *archivo = fopen(nombre, "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear %s\n", nombre);
        return;
    }

    fprintf(archivo, "x");
    for (int j = 0; j < N; j++) {
        fprintf(archivo, "\ty%d", j + 1);
    }
    for (int j = 0; j < N; j++) {
        fprintf(archivo, "\ty%d_exacta", j + 1);
    }
    for (int j = 0; j < N; j++) {
        fprintf(archivo, "\terror%d", j + 1);
    }
    fprintf(archivo, "\n");

    for (int i = 0; i <= n; i++) {
        double ye[N];
        exacta(X[i], ye);

        fprintf(archivo, "%lf", X[i]);
        for (int j = 0; j < N; j++) {
            fprintf(archivo, "\t%lf", Y[i][j]);
        }
        for (int j = 0; j < N; j++) {
            fprintf(archivo, "\t%lf", ye[j]);
        }
        for (int j = 0; j < N; j++) {
            fprintf(archivo, "\t%le", fabs(ye[j] - Y[i][j]));
        }
        fprintf(archivo, "\n");
    }

    fclose(archivo);
}

/* ==========================================================================
 * FACTOR DE CONVERGENCIA Q
 * ========================================================================== */

/*
 * Para estimar el orden de convergencia se resuelve el mismo sistema con:
 *      h, h/2 y h/4
 *
 * Luego, para cada componente j, se calcula:
 *
 *      Q = log( |Yh - Yh2| / |Yh2 - Yh4| ) / log(2)
 *
 * Interpretación:
 *      Q cercano a 1  -> método de orden 1 (Euler)
 *      Q cercano a 2  -> método de orden 2 (Heun, Punto Medio)
 *      Q cercano a 4  -> método de orden 4 (RK4)
 *
 * Nota: Yh, Yh2 e Yh4 son los valores numéricos en el mismo x.
 *       Para compararlos se usan los índices 2*i y 4*i.
 *       Si la diferencia es casi cero (por ejemplo, cuando el método es
 *       exacto para el ejemplo), Q se muestra como N/A.
 */
void factor_convergencia(int metodo, int n1, double h1, double x0, const double y0[]) {
    if (n1 > MAX_PASOS) {
        printf("Error: n es demasiado grande para el cálculo de convergencia.\n");
        return;
    }

    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    /*
     * Arreglos para las tres corridas:
     *      h   -> n1 pasos
     *      h/2 -> 2*n1 pasos
     *      h/4 -> 4*n1 pasos
     */
    double Xh[MAX_PASOS + 1], Yh[MAX_PASOS + 1][N];
    double Xh2[2 * MAX_PASOS + 1], Yh2[2 * MAX_PASOS + 1][N];
    double Xh4[4 * MAX_PASOS + 1], Yh4[4 * MAX_PASOS + 1][N];

    resolver(metodo, x0, h1, n1, y0, Xh, Yh);
    resolver(metodo, x0, h2, 2 * n1, y0, Xh2, Yh2);
    resolver(metodo, x0, h3, 4 * n1, y0, Xh4, Yh4);

    double Q[MAX_PASOS + 1][N];

    printf("\n--- Factor de convergencia (%s) ---\n", nombre_metodo(metodo));
    printf("%-8s %-15s", "i", "x_i");
    for (int j = 0; j < N; j++) {
        printf("Q%-13d", j + 1);
    }
    printf("\n");

    for (int i = 0; i <= n1; i++) {
        for (int j = 0; j < N; j++) {
            Q[i][j] = 0.0;
        }
    }

    for (int i = 1; i <= n1; i++) {
        /*
         * En el paso i con h, el mismo x corresponde a:
         *      2*i pasos con h/2
         *      4*i pasos con h/4
         */
        int idx2 = 2 * i;
        int idx4 = 4 * i;

        printf("%-8d %-15lf", i, Xh[i]);

        for (int j = 0; j < N; j++) {
            double num = fabs(Yh[i][j] - Yh2[idx2][j]);
            double den = fabs(Yh2[idx2][j] - Yh4[idx4][j]);

            /*
             * Si el numerador o el denominador son prácticamente cero,
             * la razón de errores no es confiable y se reporta N/A.
             */
            if (den > 1e-12 && num > 1e-12) {
                Q[i][j] = log(num / den) / log(2.0);
                printf("%-15lf", Q[i][j]);
            } else {
                Q[i][j] = 0.0;
                printf("%-15s", "N/A");
            }
        }
        printf("\n");
    }

    printf("Orden esperado del metodo: %d\n", orden_metodo(metodo));

    char nombre_archivo[64];
    snprintf(nombre_archivo, sizeof(nombre_archivo), "convergence_%s.txt", archivo_metodo(metodo));

    FILE *archivo = fopen(nombre_archivo, "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear %s\n", nombre_archivo);
        return;
    }

    fprintf(archivo, "i\tx");
    for (int j = 0; j < N; j++) {
        fprintf(archivo, "\tQ%d", j + 1);
    }
    fprintf(archivo, "\n");

    for (int i = 1; i <= n1; i++) {
        fprintf(archivo, "%d\t%lf", i, Xh[i]);
        for (int j = 0; j < N; j++) {
            fprintf(archivo, "\t%lf", Q[i][j]);
        }
        fprintf(archivo, "\n");
    }

    fclose(archivo);
}

/* ==========================================================================
 * PROGRAMA PRINCIPAL
 * ========================================================================== */

int main(void) {
    double x0, xf, h;
    int n;
    double y0[N];

    double X[MAX_PASOS + 1];
    double Y[MAX_PASOS + 1][N];

    char continuar = 'n';

    do {
        printf("\n========================================\n");
        printf("SISTEMA DE %d ECUACIONES DIFERENCIALES\n", N);
        printf("========================================\n\n");

        /*
         * 1) Intervalo de integración.
         */
        printf("Ingrese x0 y xf: ");
        if (scanf("%lf %lf", &x0, &xf) != 2 || xf <= x0) {
            printf("Error: debe cumplirse xf > x0.\n");
            return 1;
        }

        /*
         * 2) Cantidad de subintervalos n o tamaño de paso h.
         *    En ambos casos se termina calculando h = (xf - x0) / n.
         */
        printf("¿Desea ingresar n (1) o h (2)? ");
        int opc_paso;
        if (scanf("%d", &opc_paso) != 1) {
            printf("Error: opcion invalida.\n");
            return 1;
        }

        if (opc_paso == 1) {
            printf("Ingrese n (1 a %d): ", MAX_PASOS);
            if (scanf("%d", &n) != 1 || n <= 0 || n > MAX_PASOS) {
                printf("Error: n invalido.\n");
                return 1;
            }
            h = (xf - x0) / n;
        } else if (opc_paso == 2) {
            printf("Ingrese h: ");
            if (scanf("%lf", &h) != 1 || h <= 0) {
                printf("Error: h invalido.\n");
                return 1;
            }

            double n_real = (xf - x0) / h;
            if (n_real > (double)MAX_PASOS) {
                printf("Error: el n calculado excede el maximo permitido (%d).\n", MAX_PASOS);
                return 1;
            }

            n = (int)round(n_real);
            if (n <= 0) {
                printf("Error: h es demasiado grande para el intervalo.\n");
                return 1;
            }
        } else {
            printf("Error: opcion invalida.\n");
            return 1;
        }

        /*
         * 3) Condiciones iniciales: y1(x0), y2(x0), ..., yN(x0).
         */
        printf("Ingrese las condiciones iniciales:\n");
        for (int j = 0; j < N; j++) {
            printf("y%d(x0) = ", j + 1);
            if (scanf("%lf", &y0[j]) != 1) {
                printf("Error: condicion inicial invalida.\n");
                return 1;
            }
        }

        /*
         * 4) Selección del método numérico.
         */
        printf("Seleccione el metodo:\n");
        printf("1. Euler\n");
        printf("2. Heun\n");
        printf("3. Punto Medio\n");
        printf("4. RK4\n");
        printf("Opcion: ");

        int metodo;
        if (scanf("%d", &metodo) != 1 || metodo < 1 || metodo > 4) {
            printf("Error: metodo invalido.\n");
            return 1;
        }

        /*
         * 5) Factor de convergencia opcional.
         *    Se calcula antes de la corrida principal.
         */
        printf("¿Desea calcular el factor de convergencia? (1.Si 2.No): ");
        int opc_convergencia;
        if (scanf("%d", &opc_convergencia) != 1) {
            printf("Error: opcion invalida.\n");
            return 1;
        }

        if (opc_convergencia == 1) {
            factor_convergencia(metodo, n, h, x0, y0);
        }

        /*
         * 6) Resolución del sistema con el método elegido.
         */
        resolver(metodo, x0, h, n, y0, X, Y);

        /*
         * 7) Tabla de valores numéricos y exactos.
         */
        printf("\nh = %lf\n\n", h);
        printf("%-6s %-12s", "i", "x");

        char etiqueta[20];
        for (int j = 0; j < N; j++) {
            snprintf(etiqueta, sizeof(etiqueta), "y%d", j + 1);
            printf("%-12s", etiqueta);
        }
        for (int j = 0; j < N; j++) {
            snprintf(etiqueta, sizeof(etiqueta), "y%d_ex", j + 1);
            printf("%-12s", etiqueta);
        }
        printf("\n");

        for (int i = 0; i <= n; i++) {
            double ye[N];
            exacta(X[i], ye);

            printf("%-6d %-12.6lf", i, X[i]);
            for (int j = 0; j < N; j++) {
                printf("%-12.6lf", Y[i][j]);
            }
            for (int j = 0; j < N; j++) {
                printf("%-12.6lf", ye[j]);
            }
            printf("\n");
        }

        /*
         * 8) Tabla de errores absolutos.
         */
        printf("\nErrores absolutos:\n");
        printf("%-6s %-12s", "i", "x");
        for (int j = 0; j < N; j++) {
            snprintf(etiqueta, sizeof(etiqueta), "error%d", j + 1);
            printf("%-15s", etiqueta);
        }
        printf("\n");

        for (int i = 0; i <= n; i++) {
            double ye[N];
            exacta(X[i], ye);

            printf("%-6d %-12.6lf", i, X[i]);
            for (int j = 0; j < N; j++) {
                printf("%-15.6le", fabs(ye[j] - Y[i][j]));
            }
            printf("\n");
        }

        /*
         * 9) Resumen final en x = xf.
         */
        double ye_final[N];
        exacta(X[n], ye_final);

        printf("\n--- Resumen final en x = %lf ---\n", X[n]);
        for (int j = 0; j < N; j++) {
            double error_abs = fabs(ye_final[j] - Y[n][j]);
            double error_rel = 0.0;
            if (fabs(ye_final[j]) > 1e-14) {
                error_rel = (error_abs / fabs(ye_final[j])) * 100.0;
            }
            printf("y%d: numerica = %.10lf | exacta = %.10lf | error = %.6le | error rel = %.6lf %%\n",
                   j + 1, Y[n][j], ye_final[j], error_abs, error_rel);
        }

        /*
         * 10) Guardado de resultados en archivo.
         */
        guardar_resultados("resultados_sistema.txt", X, Y, n);
        printf("\nResultados guardados en resultados_sistema.txt\n");

        printf("\n¿Desea realizar otro calculo? (s/n): ");
        if (scanf(" %c", &continuar) != 1) {
            continuar = 'n';
        }
    } while (continuar == 's' || continuar == 'S');

    printf("\nPrograma finalizado.\n");
    return 0;
}
