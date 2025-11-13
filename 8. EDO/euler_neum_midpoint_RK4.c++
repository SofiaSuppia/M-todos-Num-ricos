#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // Necesario para 'rename' si no está en stdio.h

#define MAX_SIZE 200 // Tamaño máximo para los arreglos de puntos

/**
 * Función para guardar los valores calculados en un archivo de texto
 * @param X Arreglo de valores de X
 * @param Y Arreglo de valores de Y
 * @param n Número de subintervalos (número de puntos - 1)
 */
void guardar_en_txt(double X[], double Y[], int n);

/**
 * Función para definir la función f(x, y) de la EDO: y' = f(x, y)
 * @param x La variable independiente
 * @param y La variable dependiente
 * @return El valor de la función en (x, y)
 */
double f(double x, double y);

/**
 * Función para definir la solución exacta y(x) para fines de comparación
 * @param x La variable independiente
 * @return El valor exacto de y en x
 */
double y_exacta(double x);

/**
 * Función para definir la derivada total de f con respecto a x
 * Esta es equivalente a la segunda derivada de y (y'')
 * @param x La variable independiente
 * @param y La variable dependiente
 * @return El valor de la derivada total en (x, y)
 */
double f_prima(double x, double y);

/**
 * Función para definir la tercera derivada de y con respecto a x (y''')
 * Usado para estimar el error de truncamiento local de orden 2.
 * @param x La variable independiente
 * @param y La variable dependiente
 * @return El valor de la tercera derivada en (x, y)
 */
double y3_prima(double x, double y);

/**
 * Función para realizar un solo paso del método Runge-Kutta de 4to orden
 * @param x Valor actual de x
 * @param y Valor actual de y
 * @param paso Tamaño del paso (h)
 * @return El nuevo valor de y después del paso RK4
 */
double paso_rk4(double x, double y, double paso);

/**
 * Función para calcular el factor de convergencia para el método de Euler (RK1)
 * @param n Número de subintervalos (para el paso h)
 * @param h Tamaño del paso h
 * @param X0 Valor inicial de x
 * @param Xf Valor final de x
 * @param Y0 Valor inicial de y
 */
void factor_convergencia_euler(int n, double h, double X0, double Xf, double Y0);

/**
 * Función para calcular el factor de convergencia para el método de Heun (RK2)
 * @param n Número de subintervalos (para el paso h)
 * @param h Tamaño del paso h
 * @param X0 Valor inicial de x
 * @param Xf Valor final de x
 * @param Y0 Valor inicial de y
 */
void factor_convergencia_heun(int n, double h, double X0, double Xf, double Y0);

/**
 * Función para calcular el factor de convergencia para el método del Punto Medio (RK2)
 * @param n Número de subintervalos (para el paso h)
 * @param h Tamaño del paso h
 * @param X0 Valor inicial de x
 * @param Xf Valor final de x
 * @param Y0 Valor inicial de y
 */
void factor_convergencia_punto_medio(int n, double h, double X0, double Xf, double Y0);


/**
 * Función para calcular el factor de convergencia para el método Runge-Kutta de 4to orden (RK4)
 * @param n Número de subintervalos (para el paso h)
 * @param h Tamaño del paso h
 * @param X0 Valor inicial de x
 * @param Xf Valor final de x
 * @param Y0 Valor inicial de y
 * */
void factor_convergencia_rk4(int n, double h, double X0, double Xf, double Y0);


/*
    Resumen de Ventajas y Desventajas de los métodos
    -----------------------------------------------
    Método de Euler (Orden 1):
        Extrema simplicidad: Fácil de implementar y entender.
        Bajo costo computacional: 1 evaluación de f(x,y) por paso.
        Baja precisión: Error de truncamiento local O(h²).
        Inestabilidad: Puede divergir fácilmente con un 'h' grande.

    Método de Heun (Orden 2):
        Buena precisión: Error O(h³) vs. O(h²) de Euler.
        Más estable: Menos propenso a la divergencia.
        Equilibrado: Buen balance precisión-costo.
        Costo moderado: 2 evaluaciones de f(x,y) vs 1 de Euler.

    Método del Punto Medio (Orden 2):
        Alta precisión: Generalmente más preciso que Heun para el mismo orden.
        Mejor para simetrías: Excelente para problemas con comportamiento simétrico.
        Costo similar a Heun: 2 evaluaciones de f(x,y).
        Implementación ligeramente más compleja: Cálculo del punto medio.

    Método Runge-Kutta 4 (Orden 4):
        Muy alta precisión: Error de truncamiento local O(h⁵).
        Gran estabilidad.
        Método estándar para la mayoría de las aplicaciones.
        Alto costo: 4 evaluaciones de f(x,y) por paso.
*/

int main(int argc, char const *argv[]) {
    // Variables generales de la EDO
    double X0, Xf, Y0, h;
    // Variables temporales para cálculos (usadas en Heun/Punto Medio)
    double Xp, Yp;
    int n;
    // Arreglos para almacenar los resultados
    double X[MAX_SIZE + 1], Y[MAX_SIZE + 1];

    // Variables para selección de opciones
    int opc_convergencia; // Número para seleccionar el cálculo del factor de convergencia

    // Variables para el cálculo de errores
    double error_exacto, error_trunc_local;
    
    char continuar;
    
    do {
        printf("\n========================================\n");
        printf("SOLUCION NUMERICA DE ECUACIONES DIFERENCIALES ORDINARIAS\n");
        printf("========================================\n\n");

        printf("Inserte X0 y Xf (inicio y fin del intervalo):\n");
        if (scanf("%lf %lf", &X0, &Xf) != 2 || Xf <= X0) {
            printf("Error: Entrada inválida para X0 y Xf.\n");
            return 1;
        }

    printf("Inserte el dato inicial Y0 = Y(X0):\n");
    if (scanf("%lf", &Y0) != 1) {
        printf("Error: Entrada inválida para Y0.\n");
        return 1;
    }

    printf("¿Desea insertar el número de subintervalos (n) o el tamaño de paso (h)?\n");
    printf("1. Quiero insertar n\n");
    printf("2. Quiero insertar h\n");
    int opc_paso;
    if (scanf("%d", &opc_paso) != 1) {
        printf("Error: Entrada inválida.\n");
        return 1;
    }

    if(opc_paso == 1) {
        printf("Inserte el número de subintervalos n (entero):\n");
        if (scanf("%d", &n) != 1 || n <= 0) {
            printf("Error: Entrada inválida para n.\n");
            return 1;
        }
        // Calcular la distancia entre puntos (h)
        h = (Xf - X0) / n;
    } else if(opc_paso == 2) {
        printf("Inserte el tamaño de paso h:\n");
        if (scanf("%lf", &h) != 1 || h <= 0) {
            printf("Error: Entrada inválida para h.\n");
            return 1;
        }
        n = (int)round((Xf - X0) / h); // Redondear n al entero más cercano
        if (n > MAX_SIZE) {
            printf("Error: el número de subintervalos (%d) excede el tamaño máximo permitido (%d).\n", n, MAX_SIZE);
            return 1;
        } else if (n <= 0) {
             printf("Error: El paso h es demasiado grande o el intervalo es muy pequeño.\n");
            return 1;
        }
    } else {
        printf("Opción de paso no válida.\n");
        return 1;
    }

    // Inicializar la solución
    X[0] = X0;
    Y[0] = Y0;

    printf("Inserte el método a usar:\n");
    printf("1. Método de Euler (R.K. 1)\n");
    printf("2. Método de Heun (R.K. 2)\n");
    printf("3. Método del Punto Medio (R.K. 2)\n");
    printf("4. Runge-Kutta de orden 4\n");
    if (scanf("%d", &opc_paso) != 1) {
        printf("Error: Entrada inválida para el método.\n");
        return 1;
    }

    switch(opc_paso) {
        case 1:
            printf("¿Desea calcular el factor de convergencia para el método de Euler? (1.Sí 2.No)\n");
            if (scanf("%d", &opc_convergencia) != 1) return 1;

            if(opc_convergencia == 1) {
               factor_convergencia_euler(n, h, X0, Xf, Y0);
            }
            // Método de Euler
            for(int i = 1; i <= n; i++) {
                // Calcular el siguiente X de forma incremental para evitar errores de redondeo
                X[i] = X[i-1] + h;
                // Fórmula de Euler
                Y[i] = Y[i-1] + h * f(X[i-1], Y[i-1]);
            }

            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n",
                   "i", "X[i]", "Y Exacta", "Y Euler", "Error Exacto", "Error Trunc. Local");
            printf("-------------------------------------------------------------------------------------------\n");

            // Si local_trunc_error < 0 => El valor calculado por Euler es menor que el exacto (subestima).
            // Si local_trunc_error > 0 => El valor calculado por Euler es mayor que el exacto (sobreestima).
            for(int i = 0; i <= n; i++) {
                error_exacto = fabs(y_exacta(X[i]) - Y[i]);
                // Error de truncamiento local para Euler: O(h²)
                error_trunc_local = (h * h / 2.0) * f_prima(X[i], y_exacta(X[i]));
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n",
                       i, X[i], y_exacta(X[i]), Y[i], error_exacto, error_trunc_local);
            }
            break;

        case 2:
            printf("¿Desea calcular el factor de convergencia para el método de Heun? (1.Sí 2.No)\n");
            if (scanf("%d", &opc_convergencia) != 1) return 1;

            if(opc_convergencia == 1) {
               factor_convergencia_heun(n, h, X0, Xf, Y0);
            }
            for(int i = 1; i <= n; i++) {
                X[i] = X[i-1] + h;
                // Paso predictor (similar a Euler)
                Yp = Y[i-1] + h * f(X[i-1], Y[i-1]);
                // Paso corrector (promedio de pendientes)
                Y[i] = Y[i-1] + (h/2.0) * (f(X[i-1], Y[i-1]) + f(X[i], Yp));
            }

            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n",
                   "i", "X[i]", "Y Exacta", "Y Heun", "Error Exacto", "Error Trunc. Local");
            printf("-------------------------------------------------------------------------------------------\n");

            for(int i = 0; i <= n; i++) {
                error_exacto = fabs(y_exacta(X[i]) - Y[i]);
                // Error de truncamiento local para Heun: O(h³)
                error_trunc_local = (pow(h, 3) / 12.0) * y3_prima(X[i], y_exacta(X[i]));
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n",
                       i, X[i], y_exacta(X[i]), Y[i], error_exacto, error_trunc_local);
            }
            break;

        case 3:
            printf("¿Desea calcular el factor de convergencia para el método del Punto Medio? (1.Sí 2.No)\n");
            if (scanf("%d", &opc_convergencia) != 1) return 1;

            if(opc_convergencia == 1) {
               factor_convergencia_punto_medio(n, h, X0, Xf, Y0);
            }

            // Método del Punto Medio (Runge-Kutta de orden 2)
            for(int i = 1; i <= n; i++) {
                X[i] = X[i-1] + h;
                // Pendiente (k1) al inicio del subintervalo
                double k1 = f(X[i-1], Y[i-1]);
                // Pendiente (k2) en el punto medio (usando el predictor de Euler)
                double k2 = f(X[i-1] + h/2.0, Y[i-1] + (h/2.0) * k1);
                // Usar la pendiente del punto medio para avanzar
                Y[i] = Y[i-1] + h * k2;
            }

            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n",
                   "i", "X[i]", "Y Exacta", "Y Punto Medio", "Error Exacto", "Error Trunc. Local");
            printf("-------------------------------------------------------------------------------------------\n");

            for(int i = 0; i <= n; i++) {
                error_exacto = fabs(y_exacta(X[i]) - Y[i]);
                // Error de truncamiento local para Punto Medio: O(h³)
                error_trunc_local = (pow(h, 3) / 24.0) * y3_prima(X[i], y_exacta(X[i]));
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n",
                       i, X[i], y_exacta(X[i]), Y[i], error_exacto, error_trunc_local);
            }
            break;

        case 4:
            printf("¿Desea calcular el factor de convergencia para Runge Kutta de orden 4? (1.Sí 2.No)\n");
            if (scanf("%d", &opc_convergencia) != 1) return 1;

            if(opc_convergencia == 1) {
               factor_convergencia_rk4(n, h, X0, Xf, Y0);
            }

            double k1, k2, k3, k4;
            // Runge-Kutta de orden 4
            for(int i = 0; i <= n-1; i++) {
                k1 = f(X[i], Y[i]);
                k2 = f(X[i] + h/2.0, Y[i] + (h/2.0) * k1);
                k3 = f(X[i] + h/2.0, Y[i] + (h/2.0) * k2);
                k4 = f(X[i] + h, Y[i] + h * k3);
                Y[i+1] = Y[i] + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4);
                X[i+1] = X[i] + h;
            }

            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n",
                   "i", "X[i]", "Y Exacta", "Y RK4", "Error Exacto", "Error Trunc. Local");
            printf("-------------------------------------------------------------------------------------------\n");

            for(int i = 0; i <= n; i++) {
                error_exacto = fabs(y_exacta(X[i]) - Y[i]);

                // Estimación práctica del Error de Truncamiento Local (LTE) para RK4.
                // Se usa la comparación del paso (h) vs el doble de pasos (h/2)
                double error_local_rk4 = 0.0;
                if (i > 0) {
                    // Recalcular el paso [i-1] a [i] usando dos pasos h/2
                    double y_temp_prev = Y[i-1];
                    double x_temp_prev = X[i-1];

                    double y_mid = paso_rk4(x_temp_prev, y_temp_prev, h/2.0);
                    double y_half = paso_rk4(x_temp_prev + h/2.0, y_mid, h/2.0);

                    // Estimador de error: (Y_h/2 - Y_h) / (2^p - 1). Para RK4 (p=4), 2^4 - 1 = 15.
                    error_local_rk4 = fabs((y_half - Y[i]) / 15.0);
                }
                error_trunc_local = error_local_rk4; // Usamos el estimador

                printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                       i, X[i], y_exacta(X[i]), Y[i], error_exacto, error_trunc_local);
            }
            break;

        default:
            printf("Opción de método no válida.\n");
            return 1;
    }


    // Imprimir resultados
    printf("\nResultados finales para el archivo results.txt:\n");
    printf("X[i]\t\tY[i]\n");
    // Imprimir todos los puntos calculados
    for(int i = 0; i <= n; i++) {
        printf("%lf\t%lf\n", X[i], Y[i]);
    }

    // Guardar x[i] y Y[i] en results.txt
    guardar_en_txt(X, Y, n);

    // Finalmente, se intenta ejecutar un script de Python para graficar los resultados
    system("python ..\\graph_points.py");

    printf("\n¿Desea realizar otro calculo? (s/n): ");
    scanf(" %c", &continuar);
    
    } while(continuar == 's' || continuar == 'S');
    
    printf("\nPrograma finalizado. Gracias por usar el programa.\n");

    return 0;
}

// ------------------------------------------
// IMPLEMENTACIÓN DE FUNCIONES AUXILIARES
// ------------------------------------------

double f(double x, double y) {
    // La EDO: y' = f(x, y) = (2x+1)*sqrt(y)
    return (2.0 * x + 1.0) * sqrt(y);
}

double y_exacta(double x) {
    // Solución analítica: y(x) = ((x² + x + C)²)/4
    // Para y(0) = 1, C = 2, entonces y(x) = ((x² + x + 2)²)/4
    double temp = x * x + x + 2.0;
    return (temp * temp) / 4.0;
}


double f_prima(double x, double y) {
    // Derivada total de f con respecto a x: d(f)/dx = ∂f/∂x + ∂f/∂y * f(x, y)
    // d(y')/dx = y''
    // f(x,y) = (2x+1)*sqrt(y)
    // ∂f/∂x = 2*sqrt(y)
    // ∂f/∂y = (2x+1)/(2*sqrt(y))
    double df_dx = 2.0 * sqrt(y);
    double df_dy = (2.0 * x + 1.0) / (2.0 * sqrt(y));
    return df_dx + df_dy * f(x, y);
}


/**
 * Tercera derivada de y(x)
 * Se calcula simbólicamente o usando la fórmula general para O(h³)
 * Si la función exacta es y(x) = 1/(x²+1), entonces y''' es compleja.
 * El código original usaba la tercera derivada para y(x)=e^{-x^2}, que es:
 * y''' = 4xy(3 - 2x²)
 * Dejaré la función original por si es la que se desea usar, aunque no corresponde a y_exacta.
 */
double y3_prima(double x, double y) {
    // Tercera derivada aproximada para f(x,y) = (2x+1)*sqrt(y)
    // Cálculo simplificado basado en la cadena de derivadas
    double sqrt_y = sqrt(y);
    return 2.0 * sqrt_y + (2.0*x + 1.0)*(2.0*x + 1.0)/(4.0*sqrt_y);
}


double paso_rk4(double x, double y, double paso) {
    double k1 = f(x, y);
    double k2 = f(x + paso/2.0, y + (paso/2.0) * k1);
    double k3 = f(x + paso/2.0, y + (paso/2.0) * k2);
    double k4 = f(x + paso,      y + paso * k3);
    return y + (paso/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}


void guardar_en_txt(double X[], double Y[], int n) {
    FILE *archivo = fopen("results.txt", "w");
    if (archivo == NULL) {
        printf("Error: No se pudo crear el archivo.\n");
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", X[i], Y[i]);
    }

    fclose(archivo);
}

void guardar_en_txt_con_nombre(double X[], double Y[], int n, const char* filename) {
    FILE *archivo = fopen(filename, "w");
    if (archivo == NULL) {
        printf("Error: No se pudo crear el archivo %s.\n", filename);
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", X[i], Y[i]);
    }

    fclose(archivo);
}

// ------------------------------------------
// CÁLCULO DE FACTORES DE CONVERGENCIA (Q)
// ------------------------------------------

void factor_convergencia_euler(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;
    // La capacidad de los arrays debe ser validada para evitar desbordamiento
    if (4 * n1 + 1 > MAX_SIZE * 4) {
        printf("Error: n es demasiado grande para los cálculos de convergencia.\n");
        return;
    }

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    // Inicializaciones
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Euler con paso h
    for (int i = 0; i < n1; i++) {
        Xh[i+1] = Xh[i] + h1;
        Yh[i+1] = Yh[i] + h1 * f(Xh[i], Yh[i]);
    }

    // Euler con paso h/2
    for (int i = 0; i < 2*n1; i++) {
        Xh2[i+1] = Xh2[i] + h2;
        Yh2[i+1] = Yh2[i] + h2 * f(Xh2[i], Yh2[i]);
    }

    // Euler con paso h/4
    for (int i = 0; i < 4*n1; i++) {
        Xh4[i+1] = Xh4[i] + h3;
        Yh4[i+1] = Yh4[i] + h3 * f(Xh4[i], Yh4[i]);
    }

    printf("\n--- Factor de Convergencia (Euler) ---\n");
    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");

    // Cálculo del factor de convergencia Q en los puntos X del paso h
    Q[0] = 0.0;
    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;    // posición equivalente para h/2
        int idx4 = 4*i;    // posición equivalente para h/4
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) { // Evitar división por cero o errores de muy bajo valor
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A (Error nulo)");
        }
    }

    guardar_en_txt_con_nombre(Xh, Q, n1, "convergence_euler.txt");

    system("python ..\\graph_convergence.py");
}

void factor_convergencia_heun(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;
    if (4 * n1 + 1 > MAX_SIZE * 4) {
        printf("Error: n es demasiado grande para los cálculos de convergencia.\n");
        return;
    }

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Heun con paso h
    for (int i = 0; i < n1; i++) {
        double predictor = Yh[i] + h1 * f(Xh[i], Yh[i]);
        Yh[i+1] = Yh[i] + (h1/2.0)*(f(Xh[i], Yh[i]) + f(Xh[i]+h1, predictor));
        Xh[i+1] = Xh[i] + h1;
    }

    // Heun con paso h/2
    for (int i = 0; i < 2*n1; i++) {
        double predictor = Yh2[i] + h2 * f(Xh2[i], Yh2[i]);
        Yh2[i+1] = Yh2[i] + (h2/2.0)*(f(Xh2[i], Yh2[i]) + f(Xh2[i]+h2, predictor));
        Xh2[i+1] = Xh2[i] + h2;
    }

    // Heun con paso h/4
    for (int i = 0; i < 4*n1; i++) {
        double predictor = Yh4[i] + h3 * f(Xh4[i], Yh4[i]);
        Yh4[i+1] = Yh4[i] + (h3/2.0)*(f(Xh4[i], Yh4[i]) + f(Xh4[i]+h3, predictor));
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n--- Factor de Convergencia (Heun) ---\n");
    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A (Error nulo)");
        }
    }

    guardar_en_txt_con_nombre(Xh, Q, n1, "convergence_heun.txt");

    system("python ..\\graph_convergence.py");
}

void factor_convergencia_punto_medio(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;
    if (4 * n1 + 1 > MAX_SIZE * 4) {
        printf("Error: n es demasiado grande para los cálculos de convergencia.\n");
        return;
    }

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Punto Medio paso h
    for (int i = 0; i < n1; i++) {
        double k1 = f(Xh[i], Yh[i]);
        double k2 = f(Xh[i] + h1/2.0, Yh[i] + (h1/2.0)*k1);
        Yh[i+1] = Yh[i] + h1*k2;
        Xh[i+1] = Xh[i] + h1;
    }

    // Punto Medio paso h/2
    for (int i = 0; i < 2*n1; i++) {
        double k1 = f(Xh2[i], Yh2[i]);
        double k2 = f(Xh2[i] + h2/2.0, Yh2[i] + (h2/2.0)*k1);
        Yh2[i+1] = Yh2[i] + h2*k2;
        Xh2[i+1] = Xh2[i] + h2;
    }

    // Punto Medio paso h/4
    for (int i = 0; i < 4*n1; i++) {
        double k1 = f(Xh4[i], Yh4[i]);
        double k2 = f(Xh4[i] + h3/2.0, Yh4[i] + (h3/2.0)*k1);
        Yh4[i+1] = Yh4[i] + h3*k2;
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n--- Factor de Convergencia (Punto Medio) ---\n");
    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A (Error nulo)");
        }
    }

    guardar_en_txt_con_nombre(Xh, Q, n1, "convergence_midpoint.txt");

    system("python ..\\graph_convergence.py");
}

void factor_convergencia_rk4(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;
    if (4 * n1 + 1 > MAX_SIZE * 4) {
        printf("Error: n es demasiado grande para los cálculos de convergencia.\n");
        return;
    }

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // RK4 paso h
    for (int i = 0; i < n1; i++) {
        Yh[i+1] = paso_rk4(Xh[i], Yh[i], h1);
        Xh[i+1] = Xh[i] + h1;
    }

    // RK4 paso h/2
    for (int i = 0; i < 2*n1; i++) {
        Yh2[i+1] = paso_rk4(Xh2[i], Yh2[i], h2);
        Xh2[i+1] = Xh2[i] + h2;
    }

    // RK4 paso h/4
    for (int i = 0; i < 4*n1; i++) {
        Yh4[i+1] = paso_rk4(Xh4[i], Yh4[i], h3);
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n--- Factor de Convergencia (Runge-Kutta 4) ---\n");
    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A (Error nulo)");
        }
    }

    guardar_en_txt_con_nombre(Xh, Q, n1, "convergence_rk4.txt");

    system("python ..\\graph_convergence.py");
}