#include <stdio.h>
#include <math.h>

// Función f(x) = x^n - A para encontrar la raíz n-ésima de A
double f(double x, double A, int n) {
    return pow(x, n) - A;
}

// Derivada f'(x) = n * x^(n-1)
double df(double x, int n) {
    if (n == 1) return 1.0;
    return n * pow(x, n-1);
}

// Método de Newton-Raphson generalizado para raíz n-ésima
double newton_raphson_raiz_nesima(double A, int n, double x0, double tolerancia, int max_iter) {
    double x1, error;
    int i = 0;
    
    printf("\n=== METODO DE NEWTON-RAPHSON PARA RAIZ %d-ESIMA ===\n", n);
    printf("Buscando la raiz %d-esima de %.6lf\n", n, A);
    printf("Funcion: f(x) = x^%d - %.6lf\n", n, A);
    printf("Derivada: f'(x) = %d * x^%d\n", n, n-1);
    printf("Formula: x_{k+1} = x_k - (x_k^%d - %.6lf) / (%d * x_k^%d)\n\n", n, A, n, n-1);
    
    printf("Iter\t x_k\t\t f(x_k)\t\t Error\n");
    printf("----\t -------\t --------\t -------\n");
    
    do {
        i++;
        
        // Verificar división por cero
        double derivada = df(x0, n);
        if (fabs(derivada) < 1e-12) {
            printf("Derivada muy pequeña. El metodo puede no converger.\n");
            return x0;
        }
        
        // Fórmula de Newton-Raphson: x1 = x0 - f(x0)/f'(x0)
        x1 = x0 - f(x0, A, n) / derivada;
        
        // Calcular error absoluto
        error = fabs(x1 - x0);
        
        printf("%d\t %.10lf\t %.2e\t %.2e\n", i, x1, f(x1, A, n), error);
        
        // Actualizar para siguiente iteración
        x0 = x1;
        
    } while (error > tolerancia && i < max_iter);
    
    if (i >= max_iter) {
        printf("\nSe alcanzó el número máximo de iteraciones (%d).\n", max_iter);
    }
    
    return x1;
}

// Implementación específica para raíz cuadrada usando la fórmula del problema
double formula_raiz_cuadrada(double A, double p0, double tolerancia, int max_iter) {
    double pk = p0, pk1;
    int k = 0;
    
    printf("\n=== FORMULA ESPECIFICA PARA RAIZ CUADRADA ===\n");
    printf("Formula: p_k = (1/2) * (p_{k-1} + A/p_{k-1})\n");
    printf("Donde A = %.6lf\n\n", A);
    
    printf("Iter\t p_k\t\t Error\n");
    printf("----\t -------\t -------\n");
    
    do {
        k++;
        
        // Verificar división por cero
        if (fabs(pk) < 1e-12) {
            printf("Valor pk muy pequeño. División por cero.\n");
            return pk;
        }
        
        // Fórmula específica: pk = (1/2) * (pk-1 + A/pk-1)
        pk1 = 0.5 * (pk + A / pk);
        
        double error = fabs(pk1 - pk);
        printf("%d\t %.10lf\t %.2e\n", k, pk1, error);
        
        if (error < tolerancia) break;
        pk = pk1;
        
    } while (k < max_iter);
    
    return pk1;
}

// Implementación específica para raíz cúbica usando la fórmula del problema
double formula_raiz_cubica(double A, double p0, double tolerancia, int max_iter) {
    double pk = p0, pk1;
    int k = 0;
    
    printf("\n=== FORMULA ESPECIFICA PARA RAIZ CUBICA ===\n");
    printf("Formula: p_k = (2*p_{k-1} + A/p_{k-1}^2) / 3\n");
    printf("Donde A = %.6lf\n\n", A);
    
    printf("Iter\t p_k\t\t Error\n");
    printf("----\t -------\t -------\n");
    
    do {
        k++;
        
        // Verificar división por cero
        if (fabs(pk) < 1e-12) {
            printf("Valor pk muy pequeño. División por cero.\n");
            return pk;
        }
        
        // Fórmula específica: pk = (2*pk-1 + A/pk-1^2) / 3
        pk1 = (2.0 * pk + A / (pk * pk)) / 3.0;
        
        double error = fabs(pk1 - pk);
        printf("%d\t %.10lf\t %.2e\n", k, pk1, error);
        
        if (error < tolerancia) break;
        pk = pk1;
        
    } while (k < max_iter);
    
    return pk1;
}

void mostrar_menu() {
    printf("\n");
    printf("========================================================\n");
    printf("    PROBLEMA 3: ALGORITMO DE LA RAIZ CUADRADA Y CUBICA\n");
    printf("      Metodo de Newton-Raphson para f(x) = x^n - A\n");
    printf("========================================================\n");
    printf("1. Raiz cuadrada (n=2) - Newton-Raphson general\n");
    printf("2. Raiz cubica (n=3) - Newton-Raphson general\n");
    printf("3. Raiz cuadrada - Formula especifica del problema\n");
    printf("4. Raiz cubica - Formula especifica del problema\n");
    printf("5. Comparar Newton-Raphson vs Formula especifica (raiz cuadrada)\n");
    printf("6. Ejemplo numerico con A = 2 (sqrt(2))\n");
    printf("7. Ejemplo numerico con A = 8 (cbrt(8))\n");
    printf("0. Salir del programa\n");
    printf("========================================================\n");
    printf("Seleccione una opcion: ");
}

void continuar_programa() {
    printf("\n");
    printf("Presione Enter para continuar...");
    getchar(); // Limpiar buffer
    getchar(); // Esperar Enter
}

int main() {
    int opcion, max_iter = 100;
    double A, x0, tolerancia, resultado, resultado1, resultado2;
    
    printf("===== BIENVENIDO AL CALCULADOR DE RAICES =====\n");
    printf("Implementacion del Problema 3 - Newton-Raphson\n");
    
    do {
        mostrar_menu();
        scanf("%d", &opcion);
        
        if (opcion == 0) {
            printf("\n¡Gracias por usar el programa!\n");
            break;
        }
        
        // Solicitar datos según la opción
        if (opcion >= 1 && opcion <= 5) {
            printf("\nIngrese el valor de A: ");
            scanf("%lf", &A);
            printf("Ingrese la aproximacion inicial x0: ");
            scanf("%lf", &x0);
            printf("Ingrese la tolerancia: ");
            scanf("%lf", &tolerancia);
        }
        
        // Ejecutar la opción seleccionada
        switch(opcion) {
            case 1:
                resultado = newton_raphson_raiz_nesima(A, 2, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("Raiz cuadrada de %.6lf = %.10lf\n", A, resultado);
                printf("Verificacion: (%.10lf)^2 = %.10lf\n", resultado, pow(resultado, 2));
                printf("Error de verificacion: %.2e\n", fabs(pow(resultado, 2) - A));
                continuar_programa();
                break;
                
            case 2:
                resultado = newton_raphson_raiz_nesima(A, 3, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("Raiz cubica de %.6lf = %.10lf\n", A, resultado);
                printf("Verificacion: (%.10lf)^3 = %.10lf\n", resultado, pow(resultado, 3));
                printf("Error de verificacion: %.2e\n", fabs(pow(resultado, 3) - A));
                continuar_programa();
                break;
                
            case 3:
                resultado = formula_raiz_cuadrada(A, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("Raiz cuadrada de %.6lf = %.10lf\n", A, resultado);
                printf("Verificacion: (%.10lf)^2 = %.10lf\n", resultado, pow(resultado, 2));
                printf("Error de verificacion: %.2e\n", fabs(pow(resultado, 2) - A));
                continuar_programa();
                break;
                
            case 4:
                resultado = formula_raiz_cubica(A, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("Raiz cubica de %.6lf = %.10lf\n", A, resultado);
                printf("Verificacion: (%.10lf)^3 = %.10lf\n", resultado, pow(resultado, 3));
                printf("Error de verificacion: %.2e\n", fabs(pow(resultado, 3) - A));
                continuar_programa();
                break;
                
            case 5:
                printf("\n=== COMPARACION DE METODOS PARA RAIZ CUADRADA ===\n");
                
                printf("\n--- METODO 1: Newton-Raphson General ---");
                resultado1 = newton_raphson_raiz_nesima(A, 2, x0, tolerancia, max_iter);
                
                printf("\n--- METODO 2: Formula Especifica del Problema ---");
                resultado2 = formula_raiz_cuadrada(A, x0, tolerancia, max_iter);
                
                printf("\n=== COMPARACION FINAL ===\n");
                printf("Newton-Raphson general: %.12lf\n", resultado1);
                printf("Formula especifica:     %.12lf\n", resultado2);
                printf("Diferencia absoluta:    %.2e\n", fabs(resultado1 - resultado2));
                printf("Valor real sqrt(%.6lf): %.12lf\n", A, sqrt(A));
                continuar_programa();
                break;
                
            case 6:
                A = 2.0; x0 = 1.5; tolerancia = 1e-10;
                printf("\n=== EJEMPLO: CALCULANDO sqrt(2) ===\n");
                printf("A = %.1lf, x0 = %.1lf, tolerancia = %.0e\n", A, x0, tolerancia);
                
                resultado = newton_raphson_raiz_nesima(A, 2, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("sqrt(2) ≈ %.12lf\n", resultado);
                printf("Valor real: %.12lf\n", sqrt(2.0));
                printf("Error: %.2e\n", fabs(resultado - sqrt(2.0)));
                continuar_programa();
                break;
                
            case 7:
                A = 8.0; x0 = 2.5; tolerancia = 1e-10;
                printf("\n=== EJEMPLO: CALCULANDO cbrt(8) ===\n");
                printf("A = %.1lf, x0 = %.1lf, tolerancia = %.0e\n", A, x0, tolerancia);
                
                resultado = newton_raphson_raiz_nesima(A, 3, x0, tolerancia, max_iter);
                printf("\n=== RESULTADO ===\n");
                printf("cbrt(8) ≈ %.12lf\n", resultado);
                printf("Valor real: %.12lf\n", cbrt(8.0));
                printf("Error: %.2e\n", fabs(resultado - cbrt(8.0)));
                continuar_programa();
                break;
                
            default:
                if (opcion != 0) {
                    printf("\n❌ Opcion invalida. Por favor seleccione una opcion del 0 al 7.\n");
                    continuar_programa();
                }
                break;
        }
        
    } while (opcion != 0);
    
    printf("\n======================================================\n");
    printf("Nota: Este programa implementa las formulas de recurrencia\n");
    printf("del Problema 3 para encontrar raices usando Newton-Raphson.\n");
    printf("======================================================\n");
    
    return 0;
}
