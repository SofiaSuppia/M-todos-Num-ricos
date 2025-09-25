#include <iostream>
#include <cmath>

using namespace std;

// Prototipos de funciones
double g(double x);
double gprima(double x);
double calcular_error(double x_nuevo, double x_anterior, int tipo_error);
void mostrar_resultados(const char* metodo, double raiz, double error, int iteraciones);
void verificar_convergencia(const char* metodo, int iteraciones, double raiz, int max_iter);

/* Usando método gráfico:
    raíz f(x) = 0.6627
Usando método de Newton-Raphson:
    X0 = 0.5, tolerancia = 1e-6
    La raíz aproximada es: 0.663, con un error de 0.000
    Número de iteraciones: 4
Usando método de Punto Fijo:
    X0 = 0.5, tolerancia = 1e-6
    No converge porque |g'(x)| >= 1 en la intersección con y = x <--> x = g(x)
     */

int main(int argc, char const *argv[]) {
    double x_cero, x_uno, x_dos, error_actual, tolerancia;
    int max_iteraciones = 0, tipo_error, metodo;

    // x_cero es el valor inicial desde el cual se busca la raíz en el eje x, aproximándose con una línea recta x cada vez más cerca de la raíz
    printf("Ingrese el valor inicial x0 desde donde comenzará a buscar la raíz: ");
    scanf("%lf", &x_cero);
    printf("¿Necesita el error porcentual o absoluto? (1 para absoluto, 0 para porcentual): ");
    scanf("%d", &tipo_error);
    printf("Ingrese la tolerancia: ");
    scanf("%lf", &tolerancia);
    printf("Ingrese el tipo de método (1 para punto fijo, 2 para Newton-Raphson, 3 para secante): ");
    scanf("%d", &metodo);

    switch(metodo) {
        case 1: 
            // Método de Punto Fijo
            do {
                max_iteraciones++;

                // Si la pendiente de la curva en el punto donde se cruza la línea x con la función es mayor, no es posible encontrar la raíz en un tiempo razonable
                if(fabs(gprima(x_cero)) >= 1) {
                    printf("El método no converge en la iteración %d\n", max_iteraciones);
                    exit(0);
                }
            
                x_uno = g(x_cero);
                error_actual = calcular_error(x_uno, x_cero, tipo_error);
                x_cero = x_uno; // Actualizamos x_cero para la siguiente iteración
            } while(error_actual > tolerancia);

            mostrar_resultados("punto fijo", x_uno, error_actual, max_iteraciones);
            break;

        case 2:
            // Método de Newton-Raphson
            do {
                max_iteraciones++;
        
                // Si la derivada es muy pequeña, puede conducir a la división por cero o a una convergencia lenta
                if(fabs(gprima(x_cero)) < 10e-4) {
                    printf("La derivada es muy pequeña en la iteración %d\n", max_iteraciones);
                    exit(0);
                }
        
                // Aplicamos la fórmula de Newton-Raphson
                x_uno = x_cero - (g(x_cero) / gprima(x_cero));
                error_actual = calcular_error(x_uno, x_cero, tipo_error);
                x_cero = x_uno; // Actualizamos x_cero para la siguiente iteración
            } while(error_actual > tolerancia && max_iteraciones < 10000);
        
            mostrar_resultados("Newton-Raphson", x_uno, error_actual, max_iteraciones);
            verificar_convergencia("Newton-Raphson", max_iteraciones, x_uno, 10000);
            break;

        case 3: 
            // Método de la Secante - necesita un segundo punto inicial (debe estar cerca de la raíz, como x_cero)
            printf("Ingrese un segundo valor inicial x1 para el método de la secante: ");
            scanf("%lf", &x_uno);
            
            do {
                max_iteraciones++;
        
                // Aplicamos la fórmula del método de la secante
                x_dos = x_uno - ((g(x_uno) * (x_uno - x_cero)) / (g(x_uno) - g(x_cero)));
                error_actual = calcular_error(x_dos, x_uno, tipo_error);
        
                x_cero = x_uno; // Actualizamos x_cero (que se convierte en x_anterior)
                x_uno = x_dos; // Actualizamos x_uno (que se convierte en x_actual)
            } while(error_actual > tolerancia && max_iteraciones < 10000);
        
            mostrar_resultados("la secante", x_dos, error_actual, max_iteraciones);
            verificar_convergencia("la secante", max_iteraciones, x_dos, 10000);
            break;

        default:
            printf("Opción inválida. Seleccione 1, 2 o 3.\n");
    }
    
    return 0;
}

// Definiciones de funciones

// Función g(x) en el método de punto fijo
// En el método de Netwon-Raphson, se usa f(x); vea la diferencia conceptual en una bibliografía.
double g(double x) {
    return sin(3*x) - log(x)/2;
}

// Derivada numérica de g(x) usando diferencias finitas
// En el método de Netwon-Raphson, se usa f'(x); vea la diferencia conceptual en una bibliografía.
double gprima(double x) {
    double h = 0.01;
    // Fórmula de diferencia finita de tres puntos (aproximación)
    return (3*g(x) - 4*g(x-h) + g(x-2*h))/(2*h); 
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
    printf("La raíz aproximada es: %.6lf, con un error de %.6lf\n", raiz, error);
    printf("Número de iteraciones: %d\n", iteraciones);
    if(metodo != "punto fijo") {
        printf("Función evaluada en la raíz aproximada es igual a %.6lf\n", g(raiz));
    }
}

// Función para verificar la convergencia
// Esta función verifica si el método converge y si el valor de la función en la raíz es cercano a cero
void verificar_convergencia(const char* metodo, int iteraciones, double raiz, int max_iter) {
    double valor_funcion = g(raiz);
    
    if(iteraciones < max_iter) {
        if(fabs(valor_funcion) < 0.01) {
            printf("✓ El método %s converge correctamente.\n", metodo);
            printf("  - Convergido en %d iteraciones (< %d)\n", iteraciones, max_iter);
            printf("  - f(raíz) = %.6lf está cerca de 0\n", valor_funcion);
        } else {
            printf("⚠ El método converge en iteraciones, pero g(raíz) = %.6lf NO está cerca de 0.\n", valor_funcion);
            printf("  Posible problema: la raíz encontrada no es precisa.\n");
        }
    } else {
        printf("✗ El método %s no converge.\n", metodo);
        printf("  - Se alcanzó el número máximo de iteraciones (%d).\n", max_iter);
        printf("  - f(raíz) = %.6lf\n", valor_funcion);
    }
}