/*
 * ===========================================================================
 * PARCIAL 2 - 2017: EJERCICIO 2 - PASO 3
 * Integracion numerica de df/dx usando Trapecio y Simpson 1/3
 * ===========================================================================
 * 
 * PASO 3: Calcular I = integral de 1 a 2 de (df/dx) dx
 * Metodos: a) Trapecio Compuesto
 *          b) Simpson 1/3 Compuesto
 * 
 * Valor exacto: I = 2.0596753
 * Calcular: Error absoluto y error porcentual para ambos metodos
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_PUNTOS 11  // 11 puntos equiespaciados
#define VALOR_EXACTO 2.0596753

// ===========================================================================
// DATOS DE ENTRADA (RESULTADOS DEL PASO 2)
// ===========================================================================

// Valores de x (equiespaciados con h = 0.1)
double X[N_PUNTOS] = {
    1.0000, 1.1000, 1.2000, 1.3000, 1.4000,
    1.5000, 1.6000, 1.7000, 1.8000, 1.9000, 2.0000
};

// Valores de df/dx calculados en el PASO 2
double df_dx[N_PUNTOS] = {
   -1.4100000000,   // x = 1.0
   -1.5050000000,   // x = 1.1
   -1.6291666665,   // x = 1.2
   -1.6500000000,   // x = 1.3
   -1.6208333335,   // x = 1.4
   -1.4300000000,   // x = 1.5
   -1.2600000000,   // x = 1.6
   -0.9612500000,   // x = 1.7
   -0.5328571430,   // x = 1.8
   -0.2087500000,   // x = 1.9
   -0.0142857140    // x = 2.0
};

// ===========================================================================
// FUNCIONES DE INTEGRACION NUMERICA
// =========================================================================== 


/**
 * Metodo del Trapecio Compuesto
 * 
 * Formula: I = (h/2) * [f(x0) + 2*f(x1) + 2*f(x2) + ... + 2*f(xn-1) + f(xn)]
 * 
 * Para n subintervalos (n+1 puntos):
 * I = (h/2) * [f0 + 2*(f1 + f2 + ... + fn-1) + fn]
 */
double metodo_trapecio(double f[], int n, double h) {
    double suma = 0.0;
    
    printf("\n===============================================================\n");
    printf("METODO DEL TRAPECIO COMPUESTO\n");
    printf("===============================================================\n\n");
    
    printf("Formula: I = (h/2) * [f0 + 2*(f1+f2+...+fn-1) + fn]\n\n");
    printf("Paso h = %.4f\n", h);
    printf("Numero de subintervalos: %d\n\n", n-1);
    
    // Calcular suma de terminos intermedios
    printf("Calculando:\n");
    printf("  f(x0) = f(%.4f) = %.10f\n", X[0], f[0]);
    
    double suma_intermedios = 0.0;
    for (int i = 1; i < n-1; i++) {
        suma_intermedios += f[i];
    }
    printf("  2 * (f1 + f2 + ... + f%d) = 2 * %.10f = %.10f\n", 
           n-2, suma_intermedios, 2.0 * suma_intermedios);
    
    printf("  f(x%d) = f(%.4f) = %.10f\n\n", n-1, X[n-1], f[n-1]);
    
    // Formula del trapecio
    suma = f[0] + 2.0 * suma_intermedios + f[n-1];
    double integral = (h / 2.0) * suma;
    
    printf("Suma total = %.10f\n", suma);
    printf("I = (%.4f / 2) * %.10f\n", h, suma);
    printf("I = %.10f\n", integral);
    
    return integral;
}

/**
 * Metodo de Simpson 1/3 Compuesto
 * 
 * Formula: I = (h/3) * [f(x0) + 4*f(x1) + 2*f(x2) + 4*f(x3) + ... + 4*f(xn-1) + f(xn)]
 * 
 * Coeficientes: 1, 4, 2, 4, 2, 4, ..., 2, 4, 1
 * 
 * IMPORTANTE: Requiere n (numero de subintervalos) PAR
 */
double metodo_simpson(double f[], int n, double h) {
    printf("\n===============================================================\n");
    printf("METODO DE SIMPSON 1/3 COMPUESTO\n");
    printf("===============================================================\n\n");
    
    // Verificar que n-1 sea par (n puntos requiere n-1 subintervalos)
    if ((n-1) % 2 != 0) {
        printf("ERROR: Simpson 1/3 requiere un numero PAR de subintervalos.\n");
        printf("Tienes %d subintervalos (IMPAR). Ajustando...\n\n", n-1);
        // Usar n-1 puntos para tener n-2 subintervalos (par)
        n = n - 1;
    }
    
    printf("Formula: I = (h/3) * [f0 + 4*f1 + 2*f2 + 4*f3 + ... + 4*fn-1 + fn]\n");
    printf("Coeficientes: 1, 4, 2, 4, 2, 4, ..., 2, 4, 1\n\n");
    printf("Paso h = %.4f\n", h);
    printf("Numero de subintervalos: %d (PAR)\n\n", n-1);
    
    // Calcular suma con coeficientes de Simpson
    double suma = f[0] + f[n-1];  // Extremos con coeficiente 1
    
    printf("Calculando con coeficientes:\n");
    printf("  [1] f(x0) = %.10f\n", f[0]);
    
    double suma_impares = 0.0;
    double suma_pares = 0.0;
    
    for (int i = 1; i < n-1; i++) {
        if (i % 2 == 1) {  // Indices impares: coeficiente 4
            suma_impares += f[i];
            if (i <= 3) printf("  [4] f(x%d) = %.10f\n", i, f[i]);
        } else {  // Indices pares: coeficiente 2
            suma_pares += f[i];
            if (i <= 3) printf("  [2] f(x%d) = %.10f\n", i, f[i]);
        }
    }
    if (n > 5) printf("  ... (coeficientes similares para puntos intermedios)\n");
    printf("  [1] f(x%d) = %.10f\n\n", n-1, f[n-1]);
    
    suma += 4.0 * suma_impares + 2.0 * suma_pares;
    
    printf("Suma con coeficientes:\n");
    printf("  1*f0 = %.10f\n", f[0]);
    printf("  4*(suma impares) = 4 * %.10f = %.10f\n", suma_impares, 4.0*suma_impares);
    printf("  2*(suma pares) = 2 * %.10f = %.10f\n", suma_pares, 2.0*suma_pares);
    printf("  1*fn = %.10f\n", f[n-1]);
    printf("  Total = %.10f\n\n", suma);
    
    double integral = (h / 3.0) * suma;
    
    printf("I = (%.4f / 3) * %.10f\n", h, suma);
    printf("I = %.10f\n", integral);
    
    return integral;
}

/**
 * Calcula errores absoluto y porcentual
 */
void calcular_errores(double valor_calculado, double valor_exacto, 
                      double *error_abs, double *error_pct) {
    *error_abs = fabs(valor_exacto - valor_calculado);
    *error_pct = (*error_abs / fabs(valor_exacto)) * 100.0;
}

// ===========================================================================
// PROGRAMA PRINCIPAL
// ===========================================================================

int main() {
    double h = 0.1;  // Paso
    double I_trapecio, I_simpson;
    double error_abs_trap, error_pct_trap;
    double error_abs_simp, error_pct_simp;
    
    printf("===============================================================\n");
    printf("       PARCIAL 2 - 2017: EJERCICIO 2 - PASO 3              \n");
    printf("  Integracion Numerica: Trapecio y Simpson 1/3              \n");
    printf("===============================================================\n\n");

    // Mostrar datos de entrada
    printf("===============================================================\n");
    printf("DATOS DE ENTRADA (df/dx del PASO 2)\n");
    printf("===============================================================\n\n");
    
    printf("+-----+----------+-----------------+\n");
    printf("|  i  |    x     |     df/dx       |\n");
    printf("+-----+----------+-----------------+\n");
    for (int i = 0; i < N_PUNTOS; i++) {
        printf("| %3d | %8.4f | %15.10f |\n", i, X[i], df_dx[i]);
    }
    printf("+-----+----------+-----------------+\n\n");
    
    printf("OBJETIVO: Calcular I = integral de 1 a 2 de (df/dx) dx\n");
    printf("VALOR EXACTO: I = %.7f\n\n", VALOR_EXACTO);
    
    // ========================================================================
    // a) METODO DEL TRAPECIO
    // ========================================================================
    
    I_trapecio = metodo_trapecio(df_dx, N_PUNTOS, h);
    calcular_errores(I_trapecio, VALOR_EXACTO, &error_abs_trap, &error_pct_trap);
    
    printf("\n---------------------------------------------------------------\n");
    printf("RESULTADO TRAPECIO:\n");
    printf("  I_trapecio = %.10f\n", I_trapecio);
    printf("  I_exacto   = %.10f\n", VALOR_EXACTO);
    printf("  Error absoluto = |%.7f - %.7f| = %.10f\n", 
           VALOR_EXACTO, I_trapecio, error_abs_trap);
    printf("  Error porcentual = (%.10f / %.7f) * 100%% = %.6f%%\n",
           error_abs_trap, VALOR_EXACTO, error_pct_trap);
    printf("---------------------------------------------------------------\n");
    
    // ========================================================================
    // b) METODO DE SIMPSON 1/3
    // ========================================================================
    
    I_simpson = metodo_simpson(df_dx, N_PUNTOS, h);
    calcular_errores(I_simpson, VALOR_EXACTO, &error_abs_simp, &error_pct_simp);
    
    printf("\n---------------------------------------------------------------\n");
    printf("RESULTADO SIMPSON 1/3:\n");
    printf("  I_simpson = %.10f\n", I_simpson);
    printf("  I_exacto  = %.10f\n", VALOR_EXACTO);
    printf("  Error absoluto = |%.7f - %.7f| = %.10f\n",
           VALOR_EXACTO, I_simpson, error_abs_simp);
    printf("  Error porcentual = (%.10f / %.7f) * 100%% = %.6f%%\n",
           error_abs_simp, VALOR_EXACTO, error_pct_simp);
    printf("---------------------------------------------------------------\n");
    
    // ========================================================================
    // COMPARACION DE METODOS
    // ========================================================================
    
    printf("\n===============================================================\n");
    printf("COMPARACION DE METODOS\n");
    printf("===============================================================\n\n");
    
    printf("+------------------+------------------+------------------+------------------+\n");
    printf("|     METODO       |   I calculado    |  Error absoluto  | Error porcentual |\n");
    printf("+------------------+------------------+------------------+------------------+\n");
    printf("| Valor exacto     | %16.10f |        -         |        -         |\n", VALOR_EXACTO);
    printf("| Trapecio         | %16.10f | %16.10f | %15.6f%% |\n", 
           I_trapecio, error_abs_trap, error_pct_trap);
    printf("| Simpson 1/3      | %16.10f | %16.10f | %15.6f%% |\n",
           I_simpson, error_abs_simp, error_pct_simp);
    printf("+------------------+------------------+------------------+------------------+\n\n");
    
    // Determinar cual es mas preciso
    if (error_abs_trap < error_abs_simp) {
        printf(">>> METODO MAS PRECISO: TRAPECIO\n");
        printf("    (Error absoluto: %.10f < %.10f)\n", error_abs_trap, error_abs_simp);
    } else if (error_abs_simp < error_abs_trap) {
        printf(">>> METODO MAS PRECISO: SIMPSON 1/3\n");
        printf("    (Error absoluto: %.10f < %.10f)\n", error_abs_simp, error_abs_trap);
    } else {
        printf(">>> AMBOS METODOS TIENEN LA MISMA PRECISION\n");
    }
    printf("\n");
    
    // ========================================================================
    // ANALISIS TEORICO
    // ========================================================================
    
    printf("===============================================================\n");
    printf("ANALISIS TEORICO DE ERRORES\n");
    printf("===============================================================\n\n");
    
    printf("TRAPECIO COMPUESTO:\n");
    printf("  - Error teorico: O(h^2) = O(%.4f^2) = O(%.6f)\n", h, h*h);
    printf("  - Error real: %.10f\n", error_abs_trap);
    printf("  - Razon: Usa aproximacion lineal en cada subintervalo\n\n");
    
    printf("SIMPSON 1/3 COMPUESTO:\n");
    printf("  - Error teorico: O(h^4) = O(%.4f^4) = O(%.8f)\n", h, h*h*h*h);
    printf("  - Error real: %.10f\n", error_abs_simp);
    printf("  - Razon: Usa aproximacion parabolica (cuadratica)\n");
    printf("  - Simpson deberia ser mas preciso por tener error O(h^4)\n\n");
    
    // ========================================================================
    // RESUMEN FINAL
    // ========================================================================
    
    printf("===============================================================\n");
    printf("  RESUMEN DEL EJERCICIO COMPLETO                            \n");
    printf("===============================================================\n\n");
    
    printf("PASO 1: Interpolacion (9 puntos -> 11 puntos equiespaciados)\n");
    printf("  - Metodo: Interpolacion lineal\n");
    printf("  - Resultado: 11 puntos con h = 0.1\n\n");
    
    printf("PASO 2: Diferenciacion numerica\n");
    printf("  - Metodo: Diferencias finitas centradas\n");
    printf("  - Resultado: df/dx en 11 puntos\n\n");
    
    printf("PASO 3: Integracion numerica\n");
    printf("  - Metodo a) Trapecio: I = %.10f (Error: %.6f%%)\n", 
           I_trapecio, error_pct_trap);
    printf("  - Metodo b) Simpson:  I = %.10f (Error: %.6f%%)\n",
           I_simpson, error_pct_simp);
    printf("  - Valor exacto:       I = %.10f\n\n", VALOR_EXACTO);
    
    if (error_abs_simp < error_abs_trap) {
        printf("CONCLUSION: Simpson 1/3 es mas preciso que Trapecio\n");
        printf("            (como era de esperarse teoricamente)\n");
    } else {
        printf("NOTA: Trapecio resulto mas preciso en este caso particular\n");
        printf("      (puede deberse a la naturaleza de los datos interpolados)\n");
    }
    
    printf("\n===============================================================\n");
    printf("FIN DEL EJERCICIO 2 - PARCIAL 2 2017\n");
    printf("===============================================================\n\n");

    return 0;
}