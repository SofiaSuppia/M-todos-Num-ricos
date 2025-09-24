#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* 
Buscando las raíces de f(x) = 30809.6 * x² - 1026.99 * x³ - 2620861.9
Se puede observar que hay tres raíces. Descartamos la raíz negativa porque necesitamos encontrar una profundidad (es positiva)
Usaremos el método de Newton-Raphson:

Raíz aproximada es: 11.861, con un error de 0.000
Número de iteraciones: 3
Función evaluada en la raíz aproximada es igual a 0.000
✓ El método Newton-Raphson converge correctamente.
  - Convergió en 3 iteraciones (< 10000)
  - f(raíz) = 0.000 está cerca de 0

Raíz aproximada es: 26.314, con un error de 0.000  
Número de iteraciones: 3
Función evaluada en la raíz aproximada es igual a -0.001
✓ El método Newton-Raphson converge correctamente.
  - Convergió en 3 iteraciones (< 10000)
  - f(raíz) = -0.001 está cerca de 0

La primera raíz es la que necesitamos para encontrar la profundidad de la esfera porque 0 < d < 2 * r = 20 cm
Entonces d = 11.861 cm 
*/

// Prototipos de funciones
double f(double x);
double f_prima(double x);
double calcular_error(double x_nuevo, double x_anterior, int tipo_error);
void mostrar_resultados(const char* metodo, double raiz, double error, int iteraciones);
void verificar_convergencia(const char* metodo, int iteraciones, double raiz, int max_iter);

int main() {
    double x0, x1, x2, error, tolerancia;
    int max_iteraciones = 0, tipo_error, metodo;

    // X0 es el valor inicial desde el cual se busca la raíz en el eje x, aproximándose con una línea recta x cada vez más cerca de la raíz
    printf("Ingrese el valor inicial x0 desde donde comenzará a buscar la raíz: ");
    scanf("%lf", &x0);
    printf("¿Necesita el error porcentual o absoluto? (1 para absoluto, 0 para porcentual): ");
    scanf("%d", &tipo_error);
    printf("Ingrese la tolerancia: ");
    scanf("%lf", &tolerancia);
    printf("Ingrese el tipo de método (1 para punto fijo, 2 para Newton-Raphson, 3 para secante): ");
    scanf("%d", &metodo);

    switch(metodo) {
        case 1: 
            // Método de punto fijo
            do {
                max_iteraciones++;

                // Si la pendiente de la curva en el punto donde intersecta la línea x con la función es mayor --> no es posible encontrar la raíz en tiempo razonable
                if(fabs(f_prima(x0)) >= 1) {
                    printf("El método no converge en la iteración %d\n", max_iteraciones);
                    exit(0);
                }
            
                x1 = f(x0);
                error = calcular_error(x1, x0, tipo_error);
                x0 = x1; // Actualizamos x0 para la siguiente iteración
            } while(error > tolerancia);

            mostrar_resultados("punto fijo", x1, error, max_iteraciones);
            break;

        case 2:
            // Método de Newton-Raphson
            do {
                max_iteraciones++;
        
                // Si la derivada es muy pequeña, puede llevar a división por cero o convergencia lenta
                if(fabs(f_prima(x0)) < 1e-6) {
                    printf("La derivada es muy pequeña en la iteración %d\n", max_iteraciones);
                    exit(0);
                }
        
                // Aplicamos la fórmula de Newton-Raphson
                x1 = x0 - (f(x0) / f_prima(x0));
                error = calcular_error(x1, x0, tipo_error);
                x0 = x1; // Actualizamos x0 para la siguiente iteración
            } while(error > tolerancia && max_iteraciones < 10000);
        
            mostrar_resultados("Newton-Raphson", x1, error, max_iteraciones);
            verificar_convergencia("Newton-Raphson", max_iteraciones, x1, 10000);
            break;

        case 3: 
            // Método de la secante - necesita un segundo punto inicial (debe estar cerca de la raíz, como x0)
            printf("Ingrese un segundo valor inicial x1 para el método de la secante: ");
            scanf("%lf", &x1);
            
            do {
                max_iteraciones++;
        
                // Aplicamos la fórmula del método de la secante
                x2 = x1 - ((f(x1) * (x1 - x0)) / (f(x1) - f(x0)));
                error = calcular_error(x2, x1, tipo_error);
        
                x0 = x1; // Actualizamos x0 para la siguiente iteración
                x1 = x2; // Actualizamos x1 para la siguiente iteración
            } while(error > tolerancia && max_iteraciones < 10000);
        
            mostrar_resultados("la secante", x2, error, max_iteraciones);
            verificar_convergencia("la secante", max_iteraciones, x2, 10000);
            break;

        default:
            printf("Opción inválida. Seleccione 1, 2 o 3.\n");
    }
    
    return 0;
}

// Definiciones de funciones

// Función f(x) 
double f(double x) {
    return 30809.6 * pow(x, 2) - 1026.99 * pow(x, 3) - 2620861.9;
}

// Derivada numérica de f(x) usando diferencias finitas
double f_prima(double x) {
    return (f(x + 0.001) - f(x)) / 0.001;
}

// Función para calcular el error
double calcular_error(double x_nuevo, double x_anterior, int tipo_error) {
    if(tipo_error == 1) {
        // Error absoluto
        return fabs(x_nuevo - x_anterior);
    } else {
        // Error porcentual
        return fabs((x_nuevo - x_anterior) / x_nuevo) * 100;
    }
}

// Función para mostrar resultados
void mostrar_resultados(const char* metodo, double raiz, double error, int iteraciones) {
    printf("Raíz aproximada es: %.3lf, con un error de %.3lf\n", raiz, error);
    printf("Número de iteraciones: %d\n", iteraciones);
    if(strcmp(metodo, "punto fijo") != 0) {
        printf("Función evaluada en la raíz aproximada es igual a %.3lf\n", f(raiz));
    }
}

// Función para verificar convergencia
// Esta función verifica si el método converge y si el valor de la función en la raíz está cerca de cero
void verificar_convergencia(const char* metodo, int iteraciones, double raiz, int max_iter) {
    double valor_funcion = f(raiz);
    
    if(iteraciones < max_iter) {
        if(fabs(valor_funcion) < 0.01) {
            printf("✓ El método %s converge correctamente.\n", metodo);
            printf("  - Convergió en %d iteraciones (< %d)\n", iteraciones, max_iter);
            printf("  - f(raíz) = %.3lf está cerca de 0\n", valor_funcion);
        } else {
            printf("⚠ El método converge en iteraciones pero f(raíz) = %.3lf NO está cerca de 0.\n", valor_funcion);
            printf("  Posible problema: la raíz encontrada no es precisa.\n");
        }
    } else {
        printf("✗ El método %s no converge.\n", metodo);
        printf("  - Alcanzó el número máximo de iteraciones (%d).\n", max_iter);
        printf("  - f(raíz) = %.3lf\n", valor_funcion);
    }
}
