/*
 * ===========================================================================
 * PARCIAL 2 - 2022: EJERCICIO 1
 * Calculo de tabla, diferenciacion e integracion numerica
 * ===========================================================================
 * 
 * FUNCION: f(x) = e^(sqrt(1+x)) * ln(1+2x^2)
 * INTERVALO: [0, 1]
 * 
 * ===========================================================================
 * PARTE a) GENERAR TABLA DE VALORES
 * ===========================================================================
 * 
 * OBJETIVO: Crear una tabla que muestre los valores de x y f(x) en el 
 *           intervalo dado [0,1], en valores equiespaciados de manera que 
 *           el intervalo de interes se divida en 6 sub-intervalos.
 * 
 * PROCEDIMIENTO:
 * 
 * 1. CALCULAR EL PASO (h):
 *    - Intervalo: [a, b] = [0, 1]
 *    - Numero de subintervalos: n = 6
 *    - Formula: h = (b - a) / n
 *    - Calculo: h = (1 - 0) / 6 = 1/6 = 0.16666666...
 * 
 * 2. GENERAR PUNTOS EQUIESPACIADOS:
 *    - Numero de puntos: n + 1 = 7 puntos (incluyendo extremos)
 *    - Formula general: x_i = a + i*h, donde i = 0, 1, 2, ..., n
 *    
 *    Puntos generados:
 *    x_0 = 0 + 0*(1/6) = 0.00000000  (extremo izquierdo)
 *    x_1 = 0 + 1*(1/6) = 0.16666666
 *    x_2 = 0 + 2*(1/6) = 0.33333333
 *    x_3 = 0 + 3*(1/6) = 0.50000000  (punto medio)
 *    x_4 = 0 + 4*(1/6) = 0.66666666
 *    x_5 = 0 + 5*(1/6) = 0.83333333
 *    x_6 = 0 + 6*(1/6) = 1.00000000  (extremo derecho)
 * 
 * 3. EVALUAR LA FUNCION f(x) EN CADA PUNTO:
 *    - Funcion: f(x) = e^(sqrt(1+x)) * ln(1+2x^2)
 *    
 *    Analisis de la funcion:
 *    - Es un PRODUCTO de una exponencial y un logaritmo
 *    - Termino 1: e^(sqrt(1+x)) - exponencial con raiz cuadrada
 *    - Termino 2: ln(1+2x^2) - logaritmo natural
 *    - Dominio: x >= 0 (para que sqrt(1+x) y ln(1+2x^2) esten definidos)
 *    
 *    Pasos para evaluar f(x_i):
 *    a) Calcular: 1 + x_i
 *    b) Calcular: sqrt(1 + x_i)
 *    c) Calcular: e^(sqrt(1 + x_i))
 *    d) Calcular: 2*x_i^2
 *    e) Calcular: 1 + 2*x_i^2
 *    f) Calcular: ln(1 + 2*x_i^2)
 *    g) Multiplicar: e^(sqrt(1+x_i)) * ln(1+2*x_i^2)
 * 
 * 4. CONSTRUIR TABLA DE RESULTADOS:
 *    Formato de la tabla:
 *    +-----+-------------+------------------+
 *    |  i  |      x_i    |      f(x_i)      |
 *    +-----+-------------+------------------+
 *    |  0  |  0.00000000 |  f(0.00000000)   |
 *    |  1  |  0.16666666 |  f(0.16666666)   |
 *    |  2  |  0.33333333 |  f(0.33333333)   |
 *    |  3  |  0.50000000 |  f(0.50000000)   |
 *    |  4  |  0.66666666 |  f(0.66666666)   |
 *    |  5  |  0.83333333 |  f(0.83333333)   |
 *    |  6  |  1.00000000 |  f(1.00000000)   |
 *    +-----+-------------+------------------+
 * 
 * NOTAS IMPORTANTES:
 * - Los puntos estan EQUIESPACIADOS (distancia constante h entre ellos)
 * - Se incluyen AMBOS EXTREMOS del intervalo (x=0 y x=1)
 * - Total de puntos: 7 (porque 6 subintervalos requieren 7 puntos)
 * - Esta tabla sera la base para los calculos de diferenciacion (parte b)
 *   e integracion (partes c y d)
 * 
 * ===========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_PUNTOS 7        // 6 subintervalos = 7 puntos
#define N_SUBINTERVALOS 6

// Arreglos globales para almacenar los datos
double X[N_PUNTOS];       // Valores de x
double F_X[N_PUNTOS];     // Valores de f(x)
double DF_DX[N_PUNTOS];   // Valores de df/dx (parte b)

// ===========================================================================
// DEFINICION DE LA FUNCION f(x)
// ===========================================================================

/**
 * Funcion f(x) = e^(sqrt(1+x)) * ln(1+2x^2)
 * 
 * Esta es un PRODUCTO de dos terminos:
 *   Termino 1: e^(sqrt(1+x))
 *   Termino 2: ln(1+2x^2)
 * 
 * @param x Punto donde evaluar la funcion
 * @return Valor de f(x)
 */
double f(double x) {
    double raiz = sqrt(1.0 + x);                  // sqrt(1+x)
    double exponencial = exp(raiz);               // e^(sqrt(1+x))
    double argumento_log = 1.0 + 2.0 * x * x;     // (1+2x^2)
    double logaritmo = log(argumento_log);        // ln(1+2x^2)
    double resultado = exponencial * logaritmo;   // e^(sqrt(1+x)) * ln(1+2x^2)
    
    return resultado;
}

// ===========================================================================
// PARTE a) GENERAR TABLA DE VALORES
// ===========================================================================

void generar_tabla() {
    double a = 0.0;           // Extremo izquierdo del intervalo
    double b = 1.0;           // Extremo derecho del intervalo
    double h;                 // Paso (distancia entre puntos consecutivos)
    
    // Calcular el paso h
    h = (b - a) / N_SUBINTERVALOS;
    
    printf("===========================================================================\n");
    printf("PARTE a) GENERACION DE TABLA DE VALORES\n");
    printf("===========================================================================\n\n");
    
    printf("DATOS DEL PROBLEMA:\n");
    printf("  - Funcion: f(x) = e^(sqrt(1+x)) * ln(1+2x^2)\n");
    printf("  - Intervalo: [%.1f, %.1f]\n", a, b);
    printf("  - Numero de subintervalos: %d\n", N_SUBINTERVALOS);
    printf("  - Numero de puntos: %d\n", N_PUNTOS);
    printf("  - Paso h = (%.1f - %.1f) / %d = %.10f\n\n", b, a, N_SUBINTERVALOS, h);
    
    printf("GENERANDO PUNTOS EQUIESPACIADOS:\n");
    printf("  Formula: x_i = a + i*h, donde i = 0, 1, 2, ..., %d\n\n", N_SUBINTERVALOS);
    
    // Generar los puntos y evaluar la funcion
    for (int i = 0; i < N_PUNTOS; i++) {
        X[i] = a + i * h;           // Calcular x_i
        F_X[i] = f(X[i]);           // Evaluar f(x_i)
        
        printf("  x_%d = %.1f + %d * %.10f = %.10f\n", 
               i, a, i, h, X[i]);
    }
    
    printf("\n");
    printf("EVALUANDO LA FUNCION f(x) EN CADA PUNTO:\n\n");
    
    // Mostrar evaluaciones detalladas para algunos puntos
    printf("  Ejemplo para x_0 = %.10f:\n", X[0]);
    printf("    1 + x_0 = %.10f\n", 1.0 + X[0]);
    printf("    sqrt(1 + x_0) = %.10f\n", sqrt(1.0 + X[0]));
    printf("    e^(sqrt(1 + x_0)) = %.10f\n", exp(sqrt(1.0 + X[0])));
    printf("    2*x_0^2 = %.10f\n", 2.0 * X[0] * X[0]);
    printf("    1 + 2*x_0^2 = %.10f\n", 1.0 + 2.0 * X[0] * X[0]);
    printf("    ln(1 + 2*x_0^2) = %.10f\n", log(1.0 + 2.0 * X[0] * X[0]));
    printf("    f(x_0) = %.10f * %.10f = %.10f\n\n", 
           exp(sqrt(1.0 + X[0])), log(1.0 + 2.0 * X[0] * X[0]), F_X[0]);
    
    // Mostrar tabla completa de resultados
    printf("===========================================================================\n");
    printf("TABLA DE RESULTADOS - PARTE a)\n");
    printf("===========================================================================\n\n");
    
    printf("+-----+------------------+------------------+\n");
    printf("|  i  |       x_i        |      f(x_i)      |\n");
    printf("+-----+------------------+------------------+\n");
    
    for (int i = 0; i < N_PUNTOS; i++) {
        printf("| %3d | %16.10f | %16.10f |\n", i, X[i], F_X[i]);
    }
    
    printf("+-----+------------------+------------------+\n\n");
    
    printf("OBSERVACIONES:\n");
    printf("  - Todos los puntos estan EQUIESPACIADOS con h = %.10f\n", h);
    printf("  - Se incluyen ambos extremos: x_0 = %.1f y x_%d = %.1f\n", a, N_SUBINTERVALOS, b);
    printf("  - Total de puntos: %d\n", N_PUNTOS);
    printf("  - Esta tabla sera utilizada en las partes b), c) y d)\n\n");
}

// ===========================================================================
// PARTE b) CALCULAR df/dx USANDO DIFERENCIAS FINITAS
// ===========================================================================
/*
 * OBJETIVO: A partir de la tabla del inciso a), calcule df/dx usando 
 *           aproximaciones por diferencias finitas centradas en los puntos
 *           interiores, y operadores hacia adelante o hacia atras en los
 *           extremos, segun corresponda.
 * 
 * FORMULAS DE DIFERENCIAS FINITAS:
 * 
 * 1. DIFERENCIA HACIA ADELANTE (FORWARD) - En x_0 (extremo izquierdo):
 *    Formula: f'(x_0) = [f(x_1) - f(x_0)] / h
 *    
 *    Justificacion:
 *    - En el extremo izquierdo NO hay punto anterior
 *    - Solo podemos usar el punto actual y el siguiente
 *    - Error de truncamiento: O(h) - precision de primer orden
 * 
 * 2. DIFERENCIA CENTRADA - En x_1, x_2, x_3, x_4, x_5 (puntos interiores):
 *    Formula: f'(x_i) = [f(x_{i+1}) - f(x_{i-1})] / (2h)
 *    
 *    Justificacion:
 *    - En puntos interiores tenemos punto anterior Y siguiente
 *    - Usamos ambos para mejor aproximacion
 *    - Error de truncamiento: O(h^2) - precision de segundo orden
 *    - MAS PRECISA que forward/backward
 * 
 * 3. DIFERENCIA HACIA ATRAS (BACKWARD) - En x_6 (extremo derecho):
 *    Formula: f'(x_6) = [f(x_6) - f(x_5)] / h
 *    
 *    Justificacion:
 *    - En el extremo derecho NO hay punto siguiente
 *    - Solo podemos usar el punto actual y el anterior
 *    - Error de truncamiento: O(h) - precision de primer orden
 * 
 * PROCEDIMIENTO:
 * 1. Calcular h (ya tenemos de la parte a)
 * 2. Para i=0: usar diferencia hacia adelante
 * 3. Para i=1,2,3,4,5: usar diferencia centrada
 * 4. Para i=6: usar diferencia hacia atras
 * 5. Almacenar resultados en arreglo DF_DX[]
 * 6. Agregar columna df/dx a la tabla del inciso a)
 */

void calcular_derivadas() {
    double h = (X[N_PUNTOS-1] - X[0]) / N_SUBINTERVALOS;
    
    printf("===========================================================================\n");
    printf("PARTE b) CALCULO DE df/dx USANDO DIFERENCIAS FINITAS\n");
    printf("===========================================================================\n\n");
    
    printf("FORMULAS A UTILIZAR:\n");
    printf("  1. Forward (x_0):    f'(x_0) = [f(x_1) - f(x_0)] / h\n");
    printf("  2. Centrada (x_i):   f'(x_i) = [f(x_{i+1}) - f(x_{i-1})] / (2h)\n");
    printf("  3. Backward (x_6):   f'(x_6) = [f(x_6) - f(x_5)] / h\n\n");
    
    printf("PARAMETROS:\n");
    printf("  - Paso h = %.10f\n\n", h);
    
    printf("CALCULANDO DERIVADAS:\n\n");
    
    // 1. Punto inicial (x_0): Diferencia hacia adelante
    printf("  i = 0 (extremo izquierdo) - FORWARD:\n");
    printf("    f'(x_0) = [f(x_1) - f(x_0)] / h\n");
    printf("    f'(%.10f) = [%.10f - %.10f] / %.10f\n", 
           X[0], F_X[1], F_X[0], h);
    DF_DX[0] = (F_X[1] - F_X[0]) / h;
    printf("    f'(x_0) = %.10f\n\n", DF_DX[0]);
    
    // 2. Puntos interiores (x_1 a x_5): Diferencia centrada
    for (int i = 1; i < N_PUNTOS - 1; i++) {
        printf("  i = %d (punto interior) - CENTRADA:\n", i);
        printf("    f'(x_%d) = [f(x_%d) - f(x_%d)] / (2h)\n", i, i+1, i-1);
        printf("    f'(%.10f) = [%.10f - %.10f] / (2 * %.10f)\n",
               X[i], F_X[i+1], F_X[i-1], h);
        DF_DX[i] = (F_X[i+1] - F_X[i-1]) / (2.0 * h);
        printf("    f'(x_%d) = %.10f\n\n", i, DF_DX[i]);
    }
    
    // 3. Punto final (x_6): Diferencia hacia atras
    int ultimo = N_PUNTOS - 1;
    printf("  i = %d (extremo derecho) - BACKWARD:\n", ultimo);
    printf("    f'(x_%d) = [f(x_%d) - f(x_%d)] / h\n", ultimo, ultimo, ultimo-1);
    printf("    f'(%.10f) = [%.10f - %.10f] / %.10f\n",
           X[ultimo], F_X[ultimo], F_X[ultimo-1], h);
    DF_DX[ultimo] = (F_X[ultimo] - F_X[ultimo-1]) / h;
    printf("    f'(x_%d) = %.10f\n\n", ultimo, DF_DX[ultimo]);
    
    // Mostrar tabla completa con f(x) y df/dx
    printf("===========================================================================\n");
    printf("TABLA COMPLETA - PARTES a) Y b)\n");
    printf("===========================================================================\n\n");
    
    printf("+-----+------------------+------------------+------------------+----------+\n");
    printf("|  i  |       x_i        |      f(x_i)      |     df/dx(x_i)   |  Metodo  |\n");
    printf("+-----+------------------+------------------+------------------+----------+\n");
    
    for (int i = 0; i < N_PUNTOS; i++) {
        const char* metodo;
        if (i == 0) {
            metodo = "Forward";
        } else if (i == N_PUNTOS - 1) {
            metodo = "Backward";
        } else {
            metodo = "Centrada";
        }
        
        printf("| %3d | %16.10f | %16.10f | %16.10f | %8s |\n", 
               i, X[i], F_X[i], DF_DX[i], metodo);
    }
    
    printf("+-----+------------------+------------------+------------------+----------+\n\n");
    
    printf("OBSERVACIONES:\n");
    printf("  - Diferencia centrada es MAS PRECISA (error O(h^2))\n");
    printf("  - Forward/Backward tienen error O(h), solo en extremos\n");
    printf("  - Los valores de df/dx seran usados en las partes c) y d)\n\n");
}

// ===========================================================================
// PARTE c) CALCULAR INTEGRAL USANDO SIMPSON 1/3
// ===========================================================================
/*
 * OBJETIVO: Calcular integral de 0 a 1 de (df/dx) dx, usando el metodo de
 *           Simpson 1/3 a partir de los valores obtenidos en b).
 *           Escriba el valor obtenido para la integral, y calcule el 
 *           error absoluto exacto.
 * 
 * FORMULA DE SIMPSON 1/3 COMPUESTO:
 * 
 *   I = (h/3) * [f₀ + 4f₁ + 2f₂ + 4f₃ + 2f₄ + 4f₅ + f₆]
 * 
 * Patron de coeficientes: 1, 4, 2, 4, 2, 4, 1
 * 
 * REQUISITO: Numero de subintervalos debe ser PAR
 *   - Tenemos 6 subintervalos (PAR) ✓
 *   - Tenemos 7 puntos
 * 
 * VALOR EXACTO DE LA INTEGRAL:
 * Por el Teorema Fundamental del Calculo:
 *   ∫₀¹ (df/dx) dx = f(1) - f(0)
 * 
 * ERROR ABSOLUTO:
 *   Error = |I_calculado - I_exacto|
 */

void calcular_integral_simpson() {
    double h = (X[N_PUNTOS-1] - X[0]) / N_SUBINTERVALOS;
    double I_simpson;
    double I_exacto;
    double error_absoluto;
    
    printf("===========================================================================\n");
    printf("PARTE c) INTEGRACION CON SIMPSON 1/3\n");
    printf("===========================================================================\n\n");
    
    printf("OBJETIVO: Calcular integral de 0 a 1 de (df/dx) dx\n\n");
    
    printf("FORMULA DE SIMPSON 1/3 COMPUESTO:\n");
    printf("  I = (h/3) * [f₀ + 4f₁ + 2f₂ + 4f₃ + 2f₄ + 4f₅ + f₆]\n");
    printf("  Coeficientes: 1, 4, 2, 4, 2, 4, 1\n\n");
    
    printf("VERIFICACION:\n");
    printf("  - Numero de subintervalos: %d (PAR) ✓\n", N_SUBINTERVALOS);
    printf("  - Numero de puntos: %d\n", N_PUNTOS);
    printf("  - Paso h = %.10f\n\n", h);
    
    printf("CALCULANDO CON LOS VALORES DE df/dx:\n\n");
    
    // Aplicar formula de Simpson
    double suma = DF_DX[0] + DF_DX[N_PUNTOS-1];  // Extremos con coeficiente 1
    
    printf("  Extremos (coef = 1):\n");
    printf("    1 * df/dx(x₀) = 1 * %.10f = %.10f\n", DF_DX[0], DF_DX[0]);
    printf("    1 * df/dx(x₆) = 1 * %.10f = %.10f\n\n", DF_DX[6], DF_DX[6]);
    
    printf("  Puntos interiores:\n");
    double suma_impares = 0.0;
    double suma_pares = 0.0;
    
    for (int i = 1; i < N_PUNTOS - 1; i++) {
        if (i % 2 == 1) {  // Indices impares: coeficiente 4
            suma_impares += DF_DX[i];
            printf("    4 * df/dx(x%d) = 4 * %.10f = %.10f\n", i, DF_DX[i], 4*DF_DX[i]);
        } else {  // Indices pares: coeficiente 2
            suma_pares += DF_DX[i];
            printf("    2 * df/dx(x%d) = 2 * %.10f = %.10f\n", i, DF_DX[i], 2*DF_DX[i]);
        }
    }
    
    suma += 4.0 * suma_impares + 2.0 * suma_pares;
    
    printf("\n  SUMA TOTAL:\n");
    printf("    Extremos:       %.10f + %.10f = %.10f\n", DF_DX[0], DF_DX[6], DF_DX[0] + DF_DX[6]);
    printf("    4*(impares):    4 * %.10f = %.10f\n", suma_impares, 4.0*suma_impares);
    printf("    2*(pares):      2 * %.10f = %.10f\n", suma_pares, 2.0*suma_pares);
    printf("    Total:          %.10f\n\n", suma);
    
    I_simpson = (h / 3.0) * suma;
    
    printf("  RESULTADO:\n");
    printf("    I = (h/3) * suma\n");
    printf("    I = (%.10f / 3) * %.10f\n", h, suma);
    printf("    I = %.10f\n\n", I_simpson);
    
    // Calcular valor exacto usando Teorema Fundamental del Calculo
    I_exacto = F_X[N_PUNTOS-1] - F_X[0];
    
    printf("VALOR EXACTO (Teorema Fundamental del Calculo):\n");
    printf("  ∫₀¹ (df/dx) dx = f(1) - f(0)\n");
    printf("  I_exacto = %.10f - %.10f\n", F_X[N_PUNTOS-1], F_X[0]);
    printf("  I_exacto = %.10f\n\n", I_exacto);
    
    // Calcular error absoluto
    error_absoluto = fabs(I_simpson - I_exacto);
    
    printf("ERROR ABSOLUTO:\n");
    printf("  Error = |I_simpson - I_exacto|\n");
    printf("  Error = |%.10f - %.10f|\n", I_simpson, I_exacto);
    printf("  Error = %.10f\n\n", error_absoluto);
    
    printf("===========================================================================\n");
    printf("RESUMEN PARTE c):\n");
    printf("  - Integral (Simpson 1/3): %.10f\n", I_simpson);
    printf("  - Integral exacta:        %.10f\n", I_exacto);
    printf("  - Error absoluto:         %.10f\n", error_absoluto);
    printf("===========================================================================\n\n");
}

// ===========================================================================
// PARTE d) CALCULAR INTEGRAL USANDO CUADRATURA DE GAUSS (4 PUNTOS)
// ===========================================================================
/*
 * OBJETIVO: Repetir el inciso c) empleando cuadratura de Gauss con 4 puntos.
 * 
 * CUADRATURA DE GAUSS-LEGENDRE CON 4 PUNTOS:
 * 
 * La cuadratura de Gauss transforma la integral de [a,b] a [-1,1]:
 * 
 *   ∫ₐᵇ g(x) dx = ((b-a)/2) * ∫₋₁¹ g(((b-a)*t + (b+a))/2) dt
 * 
 * Para nuestro caso: [0,1] → [-1,1]
 *   x = (t + 1) / 2
 *   dx = 0.5 dt
 * 
 * Formula con 4 puntos:
 *   ∫₋₁¹ g(t) dt ≈ w₁*g(t₁) + w₂*g(t₂) + w₃*g(t₃) + w₄*g(t₄)
 * 
 * PUNTOS Y PESOS DE GAUSS-LEGENDRE (4 PUNTOS):
 *   t₁ = -0.861136311594053    w₁ = 0.347854845137454
 *   t₂ = -0.339981043584856    w₂ = 0.652145154862546
 *   t₃ = +0.339981043584856    w₃ = 0.652145154862546
 *   t₄ = +0.861136311594053    w₄ = 0.347854845137454
 * 
 * PROCEDIMIENTO:
 * 1. Transformar cada t_i a x_i = (t_i + 1) / 2
 * 2. Calcular df/dx en cada x_i transformado
 * 3. Aplicar: I ≈ 0.5 * Σ[w_i * (df/dx)(x_i)]
 */

// Funcion auxiliar para calcular df/dx en cualquier punto x
// usando interpolacion de los valores calculados
double interpolar_derivada(double x) {
    // Buscar intervalo donde esta x
    int k = 0;
    for (k = 0; k < N_PUNTOS - 1; k++) {
        if (x >= X[k] && x <= X[k+1]) {
            break;
        }
    }
    
    // Interpolacion lineal de df/dx
    double pendiente = (DF_DX[k+1] - DF_DX[k]) / (X[k+1] - X[k]);
    return DF_DX[k] + pendiente * (x - X[k]);
}

void calcular_integral_gauss() {
    // Puntos y pesos de Gauss-Legendre para n=4
    const int n_gauss = 4;
    double t_gauss[4] = {-0.861136311594053, -0.339981043584856, 
                          0.339981043584856,  0.861136311594053};
    double w_gauss[4] = { 0.347854845137454,  0.652145154862546,
                          0.652145154862546,  0.347854845137454};
    
    double I_gauss;
    double I_exacto;
    double error_absoluto;
    
    printf("===========================================================================\n");
    printf("PARTE d) INTEGRACION CON CUADRATURA DE GAUSS (4 PUNTOS)\n");
    printf("===========================================================================\n\n");
    
    printf("OBJETIVO: Calcular integral de 0 a 1 de (df/dx) dx\n");
    printf("          usando Cuadratura de Gauss-Legendre con 4 puntos\n\n");
    
    printf("TRANSFORMACION DE INTERVALO:\n");
    printf("  De [0, 1] a [-1, 1]\n");
    printf("  x = (t + 1) / 2\n");
    printf("  dx = 0.5 dt\n\n");
    
    printf("PUNTOS Y PESOS DE GAUSS-LEGENDRE (n=4):\n");
    printf("+---+--------------------+--------------------+\n");
    printf("| i |        t_i         |        w_i         |\n");
    printf("+---+--------------------+--------------------+\n");
    for (int i = 0; i < n_gauss; i++) {
        printf("| %d | %18.15f | %18.15f |\n", i+1, t_gauss[i], w_gauss[i]);
    }
    printf("+---+--------------------+--------------------+\n\n");
    
    printf("EVALUACION EN LOS PUNTOS DE GAUSS:\n\n");
    
    double suma = 0.0;
    
    for (int i = 0; i < n_gauss; i++) {
        // Transformar t a x
        double x_i = (t_gauss[i] + 1.0) / 2.0;
        
        // Calcular df/dx en x_i (interpolando)
        double df_dx_i = interpolar_derivada(x_i);
        
        // Contribucion al resultado
        double contribucion = w_gauss[i] * df_dx_i;
        suma += contribucion;
        
        printf("  Punto %d:\n", i+1);
        printf("    t_%d = %.15f\n", i+1, t_gauss[i]);
        printf("    x_%d = (%.15f + 1) / 2 = %.15f\n", i+1, t_gauss[i], x_i);
        printf("    df/dx(x_%d) = %.15f (interpolado)\n", i+1, df_dx_i);
        printf("    w_%d * df/dx(x_%d) = %.15f * %.15f = %.15f\n\n", 
               i+1, i+1, w_gauss[i], df_dx_i, contribucion);
    }
    
    printf("  SUMA TOTAL: %.15f\n\n", suma);
    
    I_gauss = 0.5 * suma;
    
    printf("  RESULTADO:\n");
    printf("    I = 0.5 * suma\n");
    printf("    I = 0.5 * %.15f\n", suma);
    printf("    I = %.15f\n\n", I_gauss);
    
    // Valor exacto
    I_exacto = F_X[N_PUNTOS-1] - F_X[0];
    
    printf("VALOR EXACTO:\n");
    printf("  I_exacto = f(1) - f(0) = %.15f\n\n", I_exacto);
    
    // Error absoluto
    error_absoluto = fabs(I_gauss - I_exacto);
    
    printf("ERROR ABSOLUTO:\n");
    printf("  Error = |I_gauss - I_exacto|\n");
    printf("  Error = |%.15f - %.15f|\n", I_gauss, I_exacto);
    printf("  Error = %.15f\n\n", error_absoluto);
    
    printf("===========================================================================\n");
    printf("RESUMEN PARTE d):\n");
    printf("  - Integral (Gauss 4 pts): %.15f\n", I_gauss);
    printf("  - Integral exacta:        %.15f\n", I_exacto);
    printf("  - Error absoluto:         %.15f\n", error_absoluto);
    printf("===========================================================================\n\n");
}

// ===========================================================================
// PROGRAMA PRINCIPAL
// ===========================================================================

int main() {
    printf("\n");
    printf("===========================================================================\n");
    printf("       PARCIAL 2 - 2022: EJERCICIO 1                                    \n");
    printf("  Calculo de tabla, diferenciacion e integracion numerica               \n");
    printf("===========================================================================\n\n");
    
    // PARTE a) Generar tabla de valores
    generar_tabla();
    
    // PARTE b) Calcular derivadas con diferencias finitas
    calcular_derivadas();
    
    // PARTE c) Calcular integral con Simpson 1/3
    calcular_integral_simpson();
    
    // PARTE d) Calcular integral con Gauss (4 puntos)
    calcular_integral_gauss();
    
    printf("===========================================================================\n");
    printf("FIN DEL PROGRAMA - EJERCICIO COMPLETO\n");
    printf("===========================================================================\n\n");
    
    return 0;
}
