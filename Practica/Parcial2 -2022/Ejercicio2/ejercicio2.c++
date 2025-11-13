/*
 * ===========================================================================
 * PARCIAL 2 - 2022: EJERCICIO 2
 * Metodo modificado del Trapecio para integracion numerica
 * ===========================================================================
 * 
 * ENUNCIADO: Se propone un metodo para aproximar integrales que consiste 
 *            de una modificacion del metodo del trapecio de manera que:
 * 
 *   I = ∫ f(x)dx ≈ Σ(x_{i+1} - x_i)f(x_i) + (1/2)Σ(x_{i+1} - x_i)(f(x_{i+1}) - f(x_i))
 *                 i=0                        i=0
 * 
 * Donde n es el numero de sub-intervalos usado.
 * 
 * a) Escribir pseudo-codigo e implementarlo
 * b) Calcular ∫₀¹ (3x² + 1)dx con n = 10, 100, 1000
 *    Escribir resultado y error absoluto exacto
 * c) Repetir inciso b) usando metodo del trapecio tradicional
 * 
 * BASADO EN: NewtonCoartsTrapecio.c++ 
 * 
 * ===========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ===========================================================================
// DEFINICION DE LA FUNCION f(x)
// ===========================================================================

/**
 * MODIFICACION 1: Funcion del ejercicio
 * 
 * Funcion f(x) = 3x^2 + 1
 * 
 * Esta es la funcion especifica del Ejercicio 2 del parcial.
 * Integral exacta: ∫₀¹ (3x² + 1)dx = [x³ + x]₀¹ = 2
 * 
 * @param x Punto donde evaluar la funcion
 * @return Valor de f(x)
 */
double f(double x) {
    return 3.0 * x * x + 1.0;
}

// ===========================================================================
// PARTE a) METODO MODIFICADO DEL TRAPECIO
// ===========================================================================
/*
 * FORMULA DEL METODO MODIFICADO:
 * 
 *   I ≈ Σ(x_{i+1} - x_i)f(x_i) + (1/2)Σ(x_{i+1} - x_i)(f(x_{i+1}) - f(x_i))
 *      i=0                          i=0
 * 
 * SIMPLIFICACION (si los subintervalos son iguales, x_{i+1} - x_i = h):
 * 
 *   I ≈ h·Σf(x_i) + (h/2)·Σ(f(x_{i+1}) - f(x_i))
 *      i=0           i=0
 * 
 * Expandiendo la segunda sumatoria:
 *   Σ(f(x_{i+1}) - f(x_i)) = [f(x_1)-f(x_0)] + [f(x_2)-f(x_1)] + ... + [f(x_n)-f(x_{n-1})]
 *                          = f(x_n) - f(x_0)  (suma telescopica)
 * 
 * Por lo tanto:
 *   I ≈ h·Σf(x_i) + (h/2)·[f(x_n) - f(x_0)]
 *      i=0
 * 
 * PSEUDOCODIGO:
 * PSEUDOCODIGO:
 * 
 *   ENTRADA: f(x), a, b, n
 *   
 *   1. Calcular h = (b - a) / n
 *   2. Inicializar suma1 = 0
 *   3. Para i = 0 hasta n-1:
 *        x_i = a + i*h
 *        suma1 = suma1 + f(x_i)
 *   4. suma2 = f(b) - f(a)
 *   5. I = h * suma1 + (h/2) * suma2
 *   
 *   SALIDA: I
 */

double metodo_modificado_trapecio(double a, double b, int n) {
    double h = (b - a) / n;
    double suma1 = 0.0;
    
    // Sumar f(x_i) para i = 0 hasta n-1
    for (int i = 0; i < n; i++) {
        double x_i = a + i * h;
        suma1 += f(x_i);
    }
    
    // Calcular la diferencia telescopica
    double suma2 = f(b) - f(a);
    
    // Formula final
    double I = h * suma1 + (h / 2.0) * suma2;
    
    return I;
}

// ===========================================================================
// PARTE c) METODO DEL TRAPECIO TRADICIONAL
// ===========================================================================
/*
 * FORMULA DEL TRAPECIO COMPUESTO TRADICIONAL:
 * 
 *   I ≈ (h/2) * [f(x₀) + 2f(x₁) + 2f(x₂) + ... + 2f(xₙ₋₁) + f(xₙ)]
 * 
 * Donde:
 *   - h = (b-a)/n es el paso
 *   - x₀ = a, x₁ = a+h, ..., xₙ = b
 */

double trapecio_tradicional(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += 2.0 * f(x);
    }
    
    return (h / 2.0) * sum;
}

// ===========================================================================
// FUNCION PARA IMPRIMIR RESULTADOS
// ===========================================================================

void imprimir_resultados(int n, double I_modificado, double I_tradicional, double I_exacto) {
    double error_mod = fabs(I_modificado - I_exacto);
    double error_trad = fabs(I_tradicional - I_exacto);
    
    printf("  n = %4d:\n", n);
    printf("    Metodo Modificado:  I = %.15f   Error = %.15f\n", I_modificado, error_mod);
    printf("    Trapecio Clasico:   I = %.15f   Error = %.15f\n\n", I_tradicional, error_trad);
}

// ===========================================================================
// FUNCION PRINCIPAL
// ===========================================================================

int main() {
    printf("\n");
    printf("===========================================================================\n");
    printf("       PARCIAL 2 - 2022: EJERCICIO 2                                    \n");
    printf("  Metodo Modificado del Trapecio vs Trapecio Tradicional                \n");
    printf("===========================================================================\n\n");
    
    double a = 0.0;
    double b = 1.0;
    double I_exacto = 2.0;  // ∫₀¹ (3x² + 1)dx = [x³ + x]₀¹ = 2
    
    printf("PROBLEMA:\n");
    printf("  Funcion: f(x) = 3x² + 1\n");
    printf("  Intervalo: [%.1f, %.1f]\n", a, b);
    printf("  Valor exacto de la integral: I = %.15f\n\n", I_exacto);
    
    // ===========================================================================
    // PARTE a) PSEUDO-CODIGO
    // ===========================================================================
    
    printf("===========================================================================\n");
    printf("PARTE a) PSEUDO-CODIGO DEL METODO MODIFICADO\n");
    printf("===========================================================================\n\n");
    
    printf("FORMULA:\n");
    printf("  I = Σ(x_{i+1} - x_i)f(x_i) + (1/2)Σ(x_{i+1} - x_i)(f(x_{i+1}) - f(x_i))\n");
    printf("     i=0                            i=0\n\n");
    
    printf("SIMPLIFICADA (para h constante):\n");
    printf("  I = h·Σf(x_i) + (h/2)·[f(b) - f(a)]\n");
    printf("       i=0\n\n");
    
    printf("PSEUDOCODIGO:\n");
    printf("  1. Calcular h = (b - a) / n\n");
    printf("  2. Inicializar suma1 = 0\n");
    printf("  3. Para i = 0 hasta n-1:\n");
    printf("       x_i = a + i*h\n");
    printf("       suma1 = suma1 + f(x_i)\n");
    printf("  4. suma2 = f(b) - f(a)\n");
    printf("  5. I = h * suma1 + (h/2) * suma2\n\n");
    
    // ===========================================================================
    // PARTE b) CALCULAR CON METODO MODIFICADO
    // ===========================================================================
    
    printf("===========================================================================\n");
    printf("PARTE b) CALCULO CON METODO MODIFICADO\n");
    printf("===========================================================================\n\n");
    
    int valores_n[] = {10, 100, 1000};
    
    for (int k = 0; k < 3; k++) {
        int n = valores_n[k];
        double I_mod = metodo_modificado_trapecio(a, b, n);
        double error = fabs(I_mod - I_exacto);
        
        printf("  n = %4d:  I = %.15f   Error = %.15f\n", n, I_mod, error);
    }
    
    printf("\n");
    
    // ===========================================================================
    // PARTE c) CALCULAR CON TRAPECIO TRADICIONAL
    // ===========================================================================
    
    printf("===========================================================================\n");
    printf("PARTE c) CALCULO CON TRAPECIO TRADICIONAL\n");
    printf("===========================================================================\n\n");
    
    for (int k = 0; k < 3; k++) {
        int n = valores_n[k];
        double I_trad = trapecio_tradicional(a, b, n);
        double error = fabs(I_trad - I_exacto);
        
        printf("  n = %4d:  I = %.15f   Error = %.15f\n", n, I_trad, error);
    }
    
    printf("\n");
    
    // ===========================================================================
    // COMPARACION
    // ===========================================================================
    
    printf("===========================================================================\n");
    printf("COMPARACION DE AMBOS METODOS\n");
    printf("===========================================================================\n\n");
    
    printf("Valor exacto: I = %.15f\n\n", I_exacto);
    
    for (int k = 0; k < 3; k++) {
        int n = valores_n[k];
        double I_mod = metodo_modificado_trapecio(a, b, n);
        double I_trad = trapecio_tradicional(a, b, n);
        imprimir_resultados(n, I_mod, I_trad, I_exacto);
    }
    
    printf("===========================================================================\n");
    printf("FIN DEL EJERCICIO 2\n");
    printf("===========================================================================\n\n");
    
    return 0;
}

/*
 * ===========================================================================
 * RESUMEN DE MODIFICACIONES AL CODIGO BASE (NewtonCoartsTrapecio.c++)
 * ===========================================================================
 * 
 * MODIFICACION 1 - Implementacion del Metodo Modificado:
 *   NUEVO: Funcion metodo_modificado_trapecio()
 *   Formula: I = h·Σf(x_i) + (h/2)·[f(b) - f(a)]
 *   Ubicacion: Lineas 60-78
 * 
 * MODIFICACION 2 - Trapecio Tradicional Simplificado:
 *   ORIGINAL: Codigo en main con opciones interactivas
 *   NUEVO: Funcion trapecio_tradicional() independiente
 *   Ubicacion: Lineas 90-102
 * 
 * MODIFICACION 3 - Estructura del Programa:
 *   ELIMINADO: Menu interactivo, lectura por teclado, opciones multiples
 *   AGREGADO: Ejecucion automatica de partes a), b) y c)
 *   Ubicacion: main() completo (lineas 114-200)
 * 
 * MODIFICACION 4 - Calculos Multiples:
 *   ORIGINAL: Un solo calculo por ejecucion
 *   NUEVO: Calcula con n = 10, 100, 1000 automaticamente
 *   Compara ambos metodos lado a lado
 * 
 * ===========================================================================
 */
