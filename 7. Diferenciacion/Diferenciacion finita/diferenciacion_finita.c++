#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * ============================================================================
 * DIFERENCIACIÓN NUMÉRICA CON DIFERENCIAS FINITAS
 * ============================================================================
 *
 * Este programa sigue el pseudocódigo:
 *
 *   Definir f(x)
 *   Ingresar x0 y xfinal
 *   Ingresar h o n (mejor n), con h = (xf - x0) / n
 *   Calcular f'(xi) con:
 *       - esquema hacia adelante
 *       - esquema hacia atrás
 *       - esquema centrado
 *
 * Además se agregó:
 *   - derivada exacta y errores absoluto/relativo,
 *   - segunda derivada centrada,
 *   - menú para elegir el esquema,
 *   - archivo comparativo results.txt,
 *   - tabla opcional del factor de convergencia Q.
 *
 * ============================================================================
 * MAPA DE CÁLCULOS
 * ============================================================================
 *
 * DERIVADAS:
 *   - f(x):                  función f(double x).
 *   - f'(x) exacta:          función f_prima_exacta(double x).
 *   - f''(x) exacta:         función f_biprima_exacta(double x).
 *   - f'(x) hacia adelante:  función derivada_adelante(double x, double h).
 *   - f'(x) hacia atrás:     función derivada_atras(double x, double h).
 *   - f'(x) centrada:        función derivada_centrada(double x, double h).
 *   - f''(x) centrada:       función derivada_segunda_centrada(double x, double h).
 *
 * ERRORES:
 *   - Error absoluto:        |valor_exacto - valor_numerico|.
 *   - Error relativo:        |valor_exacto - valor_numerico| / |valor_exacto| * 100.
 *   - Se calculan dentro de calcular_adelante(), calcular_atras(),
 *     calcular_centrada(), calcular_segunda() y
 *     escribir_resultados_comparativos().
 *
 * SUMATORIAS:
 *   - En este programa NO hay una sumatoria acumulada del tipo suma += ...
 *   - Los ciclos for recorren los puntos xi y calculan un valor independiente
 *     en cada iteración.
 *   - Las únicas sumas algebraicas están dentro de las fórmulas de diferencias
 *     finitas, por ejemplo: f(x+h) - f(x-h) o f(x+h) - 2f(x) + f(x-h).
 *
 * INTEGRALES:
 *   - En este programa NO se calculan integrales.
 *
 * FACTOR DE CONVERGENCIA Q:
 *   - Se calcula en tabla_convergencia() usando:
 *         Q = log(|D_h - D_h/2| / |D_h/2 - D_h/4|) / log(2)
 *   - No es una integral ni una sumatoria; es una razón logarítmica de errores.
 *
 * ============================================================================
 * ARCHIVOS GENERADOS
 * ============================================================================
 *
 *   forward.txt     -> primera derivada hacia adelante
 *   backward.txt    -> primera derivada hacia atrás
 *   centered.txt    -> primera derivada centrada
 *   second.txt      -> segunda derivada centrada
 *   results.txt     -> tabla comparativa con los esquemas seleccionados
 *
 *   convergence_forward.txt   -> factor Q del esquema hacia adelante
 *   convergence_backward.txt  -> factor Q del esquema hacia atrás
 *   convergence_centered.txt  -> factor Q del esquema centrado
 *   convergence_second.txt    -> factor Q de la segunda derivada centrada
 *
 * ============================================================================
 */

/*
 * Función a derivar.
 *
 * Ejemplo:
 *      f(x) = 10 * exp(-x/10) * sin(2x)
 *
 * Si se cambia esta función, también hay que cambiar f_prima_exacta() y
 * f_biprima_exacta() para que los errores sigan siendo válidos.
 */
double f(double x) {
    return 10.0 * exp(-x / 10.0) * sin(2.0 * x);
}

/*
 * Derivada exacta de la función f(x).
 *
 * Para f(x) = 10 * exp(-x/10) * sin(2x):
 *
 *      f'(x) = exp(-x/10) * [-sin(2x) + 20*cos(2x)]
 *
 * Esta función se usa solamente para calcular errores.
 */
double f_prima_exacta(double x) {
    return exp(-x / 10.0) * (-sin(2.0 * x) + 20.0 * cos(2.0 * x));
}

/*
 * Segunda derivada exacta de f(x).
 *
 * Para f(x) = 10 * exp(-x/10) * sin(2x):
 *
 *      f''(x) = exp(-x/10) * [-4*cos(2x) - 39.9*sin(2x)]
 *
 * Esta función se usa solamente para calcular errores de la segunda derivada.
 */
double f_biprima_exacta(double x) {
    return exp(-x / 10.0) * (-4.0 * cos(2.0 * x) - 39.9 * sin(2.0 * x));
}

/*
 * ============================================================================
 * FÓRMULAS DE DIFERENCIAS FINITAS
 * ============================================================================
 */

/*
 * Primera derivada hacia adelante.
 *
 * Fórmula:
 *      f'(x) ≈ [f(x + h) - f(x)] / h
 *
 * Error de truncamiento: O(h).
 */
double derivada_adelante(double x, double h) {
    return (f(x + h) - f(x)) / h;
}

/*
 * Primera derivada hacia atrás.
 *
 * Fórmula:
 *      f'(x) ≈ [f(x) - f(x - h)] / h
 *
 * Error de truncamiento: O(h).
 */
double derivada_atras(double x, double h) {
    return (f(x) - f(x - h)) / h;
}

/*
 * Primera derivada centrada.
 *
 * Fórmula:
 *      f'(x) ≈ [f(x + h) - f(x - h)] / (2h)
 *
 * Error de truncamiento: O(h²).
 */
double derivada_centrada(double x, double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

/*
 * Segunda derivada centrada.
 *
 * Fórmula:
 *      f''(x) ≈ [f(x + h) - 2f(x) + f(x - h)] / h²
 *
 * Error de truncamiento: O(h²).
 */
double derivada_segunda_centrada(double x, double h) {
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

/*
 * Error relativo porcentual.
 *
 * Fórmula:
 *      |exacto - aproximado| / |exacto| * 100
 *
 * Si el valor exacto es prácticamente cero, se devuelve 0 para evitar
 * una división por cero.
 */
double error_relativo_porcentual(double exacto, double aproximado) {
    if (fabs(exacto) < 1e-14) {
        return 0.0;
    }
    return (fabs(exacto - aproximado) / fabs(exacto)) * 100.0;
}

/*
 * ============================================================================
 * CÁLCULO Y GUARDADO DE CADA ESQUEMA
 * ============================================================================
 */

/*
 * Calcula la primera derivada hacia adelante en i = 0, ..., n-1.
 *
 * En cada punto:
 *      xi = x0 + i*h
 *      f'(xi) = [f(xi+h) - f(xi)] / h
 *
 * También calcula:
 *      error absoluto = |f'(xi) exacta - f'(xi) numérica|
 *      error relativo = error absoluto / |f'(xi) exacta| * 100
 *
 * Guarda en forward.txt y muestra en pantalla.
 */
void calcular_adelante(double x0, double h, int n) {
    FILE *archivo = fopen("forward.txt", "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear forward.txt\n");
        return;
    }

    fprintf(archivo, "x\tf_prima_adelante\tf_prima_exacta\terror_absoluto\terror_relativo\n");

    printf("\n--- Esquema hacia adelante ---\n");
    printf("%-12s %-16s %-16s %-16s %-16s\n",
           "x", "f'(x) numerica", "f'(x) exacta", "error absoluto", "error relativo %");
    printf("--------------------------------------------------------------------------------------------\n");

    /*
     * Este ciclo recorre n puntos. No es una sumatoria: en cada iteración
     * se calcula un valor independiente de la derivada.
     */
    for (int i = 0; i <= n - 1; i++) {
        double xi = x0 + i * h;

        /* Aquí se calcula la derivada numérica hacia adelante. */
        double derivada = derivada_adelante(xi, h);

        /* Aquí se calcula la derivada exacta para comparar. */
        double exacta = f_prima_exacta(xi);

        /* Aquí se calculan los errores absoluto y relativo. */
        double error_abs = fabs(exacta - derivada);
        double error_rel = error_relativo_porcentual(exacta, derivada);

        fprintf(archivo, "%lf\t%lf\t%lf\t%le\t%lf\n",
                xi, derivada, exacta, error_abs, error_rel);
        printf("%-12lf %-16lf %-16lf %-16le %-16lf\n",
               xi, derivada, exacta, error_abs, error_rel);
    }

    fclose(archivo);
}

/*
 * Calcula la primera derivada hacia atrás en i = 1, ..., n.
 *
 * En cada punto:
 *      xi = x0 + i*h
 *      f'(xi) = [f(xi) - f(xi-h)] / h
 *
 * Guarda en backward.txt y muestra en pantalla.
 */
void calcular_atras(double x0, double h, int n) {
    FILE *archivo = fopen("backward.txt", "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear backward.txt\n");
        return;
    }

    fprintf(archivo, "x\tf_prima_atras\tf_prima_exacta\terror_absoluto\terror_relativo\n");

    printf("\n--- Esquema hacia atras ---\n");
    printf("%-12s %-16s %-16s %-16s %-16s\n",
           "x", "f'(x) numerica", "f'(x) exacta", "error absoluto", "error relativo %");
    printf("--------------------------------------------------------------------------------------------\n");

    /*
     * Este ciclo recorre los puntos desde i = 1 hasta i = n.
     * No es una sumatoria; calcula un valor independiente por punto.
     */
    for (int i = 1; i <= n; i++) {
        double xi = x0 + i * h;

        /* Aquí se calcula la derivada numérica hacia atrás. */
        double derivada = derivada_atras(xi, h);

        /* Derivada exacta y errores. */
        double exacta = f_prima_exacta(xi);
        double error_abs = fabs(exacta - derivada);
        double error_rel = error_relativo_porcentual(exacta, derivada);

        fprintf(archivo, "%lf\t%lf\t%lf\t%le\t%lf\n",
                xi, derivada, exacta, error_abs, error_rel);
        printf("%-12lf %-16lf %-16lf %-16le %-16lf\n",
               xi, derivada, exacta, error_abs, error_rel);
    }

    fclose(archivo);
}

/*
 * Calcula la primera derivada centrada en i = 1, ..., n-1.
 *
 * En cada punto:
 *      xi = x0 + i*h
 *      f'(xi) = [f(xi+h) - f(xi-h)] / (2h)
 *
 * Guarda en centered.txt y muestra en pantalla.
 */
void calcular_centrada(double x0, double h, int n) {
    FILE *archivo = fopen("centered.txt", "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear centered.txt\n");
        return;
    }

    fprintf(archivo, "x\tf_prima_centrada\tf_prima_exacta\terror_absoluto\terror_relativo\n");

    printf("\n--- Esquema centrado ---\n");
    printf("%-12s %-16s %-16s %-16s %-16s\n",
           "x", "f'(x) numerica", "f'(x) exacta", "error absoluto", "error relativo %");
    printf("--------------------------------------------------------------------------------------------\n");

    /*
     * Se usan solamente los puntos interiores porque la fórmula centrada
     * necesita f(xi-h) y f(xi+h).
     */
    for (int i = 1; i <= n - 1; i++) {
        double xi = x0 + i * h;

        /* Aquí se calcula la derivada numérica centrada. */
        double derivada = derivada_centrada(xi, h);

        /* Derivada exacta y errores. */
        double exacta = f_prima_exacta(xi);
        double error_abs = fabs(exacta - derivada);
        double error_rel = error_relativo_porcentual(exacta, derivada);

        fprintf(archivo, "%lf\t%lf\t%lf\t%le\t%lf\n",
                xi, derivada, exacta, error_abs, error_rel);
        printf("%-12lf %-16lf %-16lf %-16le %-16lf\n",
               xi, derivada, exacta, error_abs, error_rel);
    }

    fclose(archivo);
}

/*
 * Calcula la segunda derivada centrada en i = 1, ..., n-1.
 *
 * En cada punto:
 *      xi = x0 + i*h
 *      f''(xi) = [f(xi+h) - 2f(xi) + f(xi-h)] / h²
 *
 * Guarda en second.txt y muestra en pantalla.
 */
void calcular_segunda(double x0, double h, int n) {
    FILE *archivo = fopen("second.txt", "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear second.txt\n");
        return;
    }

    fprintf(archivo, "x\tf_biprima_centrada\tf_biprima_exacta\terror_absoluto\terror_relativo\n");

    printf("\n--- Segunda derivada centrada ---\n");
    printf("%-12s %-18s %-18s %-16s %-16s\n",
           "x", "f''(x) numerica", "f''(x) exacta", "error absoluto", "error relativo %");
    printf("--------------------------------------------------------------------------------------------\n");

    for (int i = 1; i <= n - 1; i++) {
        double xi = x0 + i * h;

        /* Aquí se calcula la segunda derivada numérica centrada. */
        double derivada = derivada_segunda_centrada(xi, h);

        /* Segunda derivada exacta y errores. */
        double exacta = f_biprima_exacta(xi);
        double error_abs = fabs(exacta - derivada);
        double error_rel = error_relativo_porcentual(exacta, derivada);

        fprintf(archivo, "%lf\t%lf\t%lf\t%le\t%lf\n",
                xi, derivada, exacta, error_abs, error_rel);
        printf("%-12lf %-18lf %-18lf %-16le %-16lf\n",
               xi, derivada, exacta, error_abs, error_rel);
    }

    fclose(archivo);
}

/*
 * Escribe results.txt con una tabla comparativa de los esquemas elegidos.
 *
 * Si un esquema no está disponible en un punto (por ejemplo, el centrado en
 * los extremos), se escribe NA en sus columnas.
 *
 * No se calculan integrales ni sumatorias acumuladas aquí.
 */
void escribir_resultados_comparativos(double x0, double h, int n,
                                      int hacer_adelante,
                                      int hacer_atras,
                                      int hacer_centrada,
                                      int hacer_segunda) {
    FILE *archivo = fopen("results.txt", "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear results.txt\n");
        return;
    }

    /* Encabezado dinámico según los esquemas seleccionados. */
    fprintf(archivo, "x");

    if (hacer_adelante) {
        fprintf(archivo, "\tf_prima_adelante\tf_prima_exacta\terror_adelante\terror_rel_adelante");
    }
    if (hacer_atras) {
        fprintf(archivo, "\tf_prima_atras\tf_prima_exacta_atras\terror_atras\terror_rel_atras");
    }
    if (hacer_centrada) {
        fprintf(archivo, "\tf_prima_centrada\tf_prima_exacta_centrada\terror_centrada\terror_rel_centrada");
    }
    if (hacer_segunda) {
        fprintf(archivo, "\tf_biprima_centrada\tf_biprima_exacta\terror_segunda\terror_rel_segunda");
    }
    fprintf(archivo, "\n");

    /*
     * Se recorre todo el intervalo, de i = 0 a i = n.
     * Cada columna se llena solamente si el esquema está definido en ese punto.
     */
    for (int i = 0; i <= n; i++) {
        double xi = x0 + i * h;

        fprintf(archivo, "%lf", xi);

        /* Primera derivada hacia adelante: válida para i = 0,...,n-1. */
        if (hacer_adelante) {
            if (i <= n - 1) {
                double derivada = derivada_adelante(xi, h);
                double exacta = f_prima_exacta(xi);
                double error_abs = fabs(exacta - derivada);
                double error_rel = error_relativo_porcentual(exacta, derivada);
                fprintf(archivo, "\t%lf\t%lf\t%le\t%lf",
                        derivada, exacta, error_abs, error_rel);
            } else {
                fprintf(archivo, "\tNA\tNA\tNA\tNA");
            }
        }

        /* Primera derivada hacia atrás: válida para i = 1,...,n. */
        if (hacer_atras) {
            if (i >= 1) {
                double derivada = derivada_atras(xi, h);
                double exacta = f_prima_exacta(xi);
                double error_abs = fabs(exacta - derivada);
                double error_rel = error_relativo_porcentual(exacta, derivada);
                fprintf(archivo, "\t%lf\t%lf\t%le\t%lf",
                        derivada, exacta, error_abs, error_rel);
            } else {
                fprintf(archivo, "\tNA\tNA\tNA\tNA");
            }
        }

        /* Primera derivada centrada: válida para i = 1,...,n-1. */
        if (hacer_centrada) {
            if (i >= 1 && i <= n - 1) {
                double derivada = derivada_centrada(xi, h);
                double exacta = f_prima_exacta(xi);
                double error_abs = fabs(exacta - derivada);
                double error_rel = error_relativo_porcentual(exacta, derivada);
                fprintf(archivo, "\t%lf\t%lf\t%le\t%lf",
                        derivada, exacta, error_abs, error_rel);
            } else {
                fprintf(archivo, "\tNA\tNA\tNA\tNA");
            }
        }

        /* Segunda derivada centrada: válida para i = 1,...,n-1. */
        if (hacer_segunda) {
            if (i >= 1 && i <= n - 1) {
                double derivada = derivada_segunda_centrada(xi, h);
                double exacta = f_biprima_exacta(xi);
                double error_abs = fabs(exacta - derivada);
                double error_rel = error_relativo_porcentual(exacta, derivada);
                fprintf(archivo, "\t%lf\t%lf\t%le\t%lf",
                        derivada, exacta, error_abs, error_rel);
            } else {
                fprintf(archivo, "\tNA\tNA\tNA\tNA");
            }
        }

        fprintf(archivo, "\n");
    }

    fclose(archivo);
}

/*
 * ============================================================================
 * FACTOR DE CONVERGENCIA Q
 * ============================================================================
 */

/*
 * Devuelve el nombre del esquema para imprimir mensajes.
 */
const char *nombre_esquema(int tipo) {
    switch (tipo) {
        case 1: return "hacia adelante";
        case 2: return "hacia atras";
        case 3: return "centrado";
        case 4: return "segunda derivada centrada";
        default: return "desconocido";
    }
}

/*
 * Devuelve el orden esperado del error:
 *      adelante y atrás -> 1
 *      centrado y segunda derivada centrada -> 2
 */
int orden_esperado(int tipo) {
    switch (tipo) {
        case 1: return 1;
        case 2: return 1;
        case 3: return 2;
        case 4: return 2;
        default: return 0;
    }
}

/*
 * Evalúa la derivada numérica según el tipo de esquema.
 * Se usa para calcular la tabla de convergencia.
 */
double evaluar_derivada(int tipo, double x, double h) {
    switch (tipo) {
        case 1: return derivada_adelante(x, h);
        case 2: return derivada_atras(x, h);
        case 3: return derivada_centrada(x, h);
        case 4: return derivada_segunda_centrada(x, h);
        default: return 0.0;
    }
}

/*
 * Calcula el factor de convergencia Q para un esquema.
 *
 * Se evalúa la derivada en el mismo x con tres pasos: h, h/2 y h/4.
 *
 * Fórmula:
 *      Q = log(|D_h - D_h/2| / |D_h/2 - D_h/4|) / log(2)
 *
 * Interpretación:
 *      Q ≈ 1 -> orden 1
 *      Q ≈ 2 -> orden 2
 *
 * No es una integral ni una sumatoria: es una razón logarítmica de errores.
 */
void tabla_convergencia(const char *nombre_archivo, int tipo,
                        double x0, double h, int n) {
    if (n < 2) {
        printf("\nNo se puede calcular convergencia: se necesita n >= 2.\n");
        return;
    }

    FILE *archivo = fopen(nombre_archivo, "w");
    if (archivo == NULL) {
        printf("Error: no se pudo crear %s\n", nombre_archivo);
        return;
    }

    fprintf(archivo, "x\tQ\torden_esperado\n");

    printf("\n--- Convergencia: %s ---\n", nombre_esquema(tipo));
    printf("%-12s %-15s %-15s\n", "x", "Q", "orden esperado");
    printf("----------------------------------------\n");

    /*
     * Se usan los puntos interiores i = 1,...,n-1 para que existan
     * los valores f(x-h) y f(x+h).
     */
    for (int i = 1; i <= n - 1; i++) {
        double xi = x0 + i * h;

        /* Tres aproximaciones de la misma derivada con pasos h, h/2 y h/4. */
        double d_h = evaluar_derivada(tipo, xi, h);
        double d_h2 = evaluar_derivada(tipo, xi, h / 2.0);
        double d_h4 = evaluar_derivada(tipo, xi, h / 4.0);

        /*
         * Numerador y denominador de la razón de errores.
         * No son sumatorias: son diferencias entre aproximaciones.
         */
        double numerador = fabs(d_h - d_h2);
        double denominador = fabs(d_h2 - d_h4);

        int orden = orden_esperado(tipo);

        /*
         * Si las diferencias son casi cero, Q no es confiable.
         * Se escribe NA en el archivo.
         */
        if (numerador > 1e-12 && denominador > 1e-12) {
            double Q = log(numerador / denominador) / log(2.0);
            fprintf(archivo, "%lf\t%lf\t%d\n", xi, Q, orden);
            printf("%-12lf %-15lf %-15d\n", xi, Q, orden);
        } else {
            fprintf(archivo, "%lf\tNA\t%d\n", xi, orden);
            printf("%-12lf %-15s %-15d\n", xi, "NA", orden);
        }
    }

    fclose(archivo);
}

/*
 * ============================================================================
 * PROGRAMA PRINCIPAL
 * ============================================================================
 */

int main(void) {
    double x0, xf, h;
    int n;
    int opcion_paso;
    int orden_derivada;

    printf("========================================\n");
    printf("DIFERENCIACION NUMERICA - DIFERENCIAS FINITAS\n");
    printf("========================================\n\n");

    /*
     * 1) Intervalo de trabajo.
     */
    printf("Ingrese x0 y xfinal: ");
    if (scanf("%lf %lf", &x0, &xf) != 2 || xf <= x0) {
        printf("Error: debe cumplirse xf > x0.\n");
        return 1;
    }

    /*
     * 2) El usuario puede ingresar n o h.
     *    Se recomienda n porque así el último punto coincide exactamente con xf.
     */
    printf("¿Desea ingresar n (1) o h (2)? (se recomienda n): ");
    if (scanf("%d", &opcion_paso) != 1) {
        printf("Error: opcion invalida.\n");
        return 1;
    }

    if (opcion_paso == 1) {
        printf("Ingrese el numero de subintervalos n: ");
        if (scanf("%d", &n) != 1 || n <= 0) {
            printf("Error: n debe ser un entero positivo.\n");
            return 1;
        }

        /* Fórmula del pseudocódigo: h = (xf - x0) / n. */
        h = (xf - x0) / n;
    } else if (opcion_paso == 2) {
        printf("Ingrese el tamano de paso h: ");
        if (scanf("%lf", &h) != 1 || h <= 0) {
            printf("Error: h debe ser positivo.\n");
            return 1;
        }

        /*
         * Si se ingresa h, se calcula n = (xf - x0) / h.
         * Luego se recalcula h para que el intervalo termine exactamente en xf.
         */
        double n_real = (xf - x0) / h;
        n = (int)round(n_real);

        if (n <= 0) {
            printf("Error: h es demasiado grande para el intervalo.\n");
            return 1;
        }

        h = (xf - x0) / n;
    } else {
        printf("Error: opcion invalida.\n");
        return 1;
    }

    /*
     * 3) Elegir qué derivada se desea calcular.
     */
    printf("\n¿Que derivada desea calcular?\n");
    printf("1. Primera derivada\n");
    printf("2. Segunda derivada centrada\n");
    printf("3. Ambas\n");
    printf("Opcion: ");
    if (scanf("%d", &orden_derivada) != 1 ||
        orden_derivada < 1 || orden_derivada > 3) {
        printf("Error: opcion invalida.\n");
        return 1;
    }

    int hacer_adelante = 0;
    int hacer_atras = 0;
    int hacer_centrada = 0;
    int hacer_segunda = 0;

    if (orden_derivada == 1 || orden_derivada == 3) {
        int opcion_esquema;

        printf("\nEsquema para la primera derivada:\n");
        printf("1. Hacia adelante\n");
        printf("2. Hacia atras\n");
        printf("3. Centrado\n");
        printf("4. Los tres esquemas\n");
        printf("Opcion: ");
        if (scanf("%d", &opcion_esquema) != 1 ||
            opcion_esquema < 1 || opcion_esquema > 4) {
            printf("Error: opcion invalida.\n");
            return 1;
        }

        hacer_adelante = (opcion_esquema == 1 || opcion_esquema == 4);
        hacer_atras = (opcion_esquema == 2 || opcion_esquema == 4);
        hacer_centrada = (opcion_esquema == 3 || opcion_esquema == 4);
    }

    if (orden_derivada == 2 || orden_derivada == 3) {
        hacer_segunda = 1;
    }

    printf("\nDatos usados:\n");
    printf("x0 = %lf\n", x0);
    printf("xf = %lf\n", xf);
    printf("n  = %d\n", n);
    printf("h  = %lf\n", h);

    if (n < 2) {
        printf("\nAdvertencia: con n < 2 no hay puntos interiores para el esquema centrado ni para la segunda derivada.\n");
    }

    /*
     * 4) Cálculo de los esquemas seleccionados.
     *    Cada función calcula la derivada, los errores y guarda su archivo.
     */
    if (hacer_adelante) {
        calcular_adelante(x0, h, n);
    }

    if (hacer_atras) {
        calcular_atras(x0, h, n);
    }

    if (hacer_centrada) {
        calcular_centrada(x0, h, n);
    }

    if (hacer_segunda) {
        calcular_segunda(x0, h, n);
    }

    /*
     * 5) Archivo comparativo con todos los esquemas seleccionados.
     */
    escribir_resultados_comparativos(x0, h, n,
                                     hacer_adelante,
                                     hacer_atras,
                                     hacer_centrada,
                                     hacer_segunda);

    /*
     * 6) Tabla opcional del factor de convergencia Q.
     */
    printf("\n¿Desea calcular la tabla de convergencia Q? (1.Si 2.No): ");
    int opcion_convergencia;
    if (scanf("%d", &opcion_convergencia) != 1) {
        printf("Error: opcion invalida.\n");
        return 1;
    }

    if (opcion_convergencia == 1) {
        if (hacer_adelante) {
            tabla_convergencia("convergence_forward.txt", 1, x0, h, n);
        }
        if (hacer_atras) {
            tabla_convergencia("convergence_backward.txt", 2, x0, h, n);
        }
        if (hacer_centrada) {
            tabla_convergencia("convergence_centered.txt", 3, x0, h, n);
        }
        if (hacer_segunda) {
            tabla_convergencia("convergence_second.txt", 4, x0, h, n);
        }
    }

    printf("\nArchivos generados:\n");
    if (hacer_adelante) {
        printf("  forward.txt\n");
    }
    if (hacer_atras) {
        printf("  backward.txt\n");
    }
    if (hacer_centrada) {
        printf("  centered.txt\n");
    }
    if (hacer_segunda) {
        printf("  second.txt\n");
    }
    printf("  results.txt\n");

    if (opcion_convergencia == 1) {
        if (hacer_adelante) {
            printf("  convergence_forward.txt\n");
        }
        if (hacer_atras) {
            printf("  convergence_backward.txt\n");
        }
        if (hacer_centrada) {
            printf("  convergence_centered.txt\n");
        }
        if (hacer_segunda) {
            printf("  convergence_second.txt\n");
        }
    }

    printf("\nPrograma finalizado.\n");
    return 0;
}
