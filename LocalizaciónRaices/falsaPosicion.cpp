#include <iostream>
#include <math.h>

// Definir la función f(x)
double f(double x) {
    return 5 * x - 3; // Ejemplo
}

int main() {
    double a, b, tolerancia, c, c_viejo, error;
    int iteraciones = 0;

    printf("Ingrese a: ");
    scanf("%lf", &a);
    printf("Ingrese b: ");
    scanf("%lf", &b);
    printf("Ingrese tolerancia (%%): ");
    scanf("%lf", &tolerancia);

    // Chequear si hay raíz en el intervalo
    if (f(a) * f(b) > 0) {
        printf("No hay raiz en el intervalo [%.6lf, %.6lf]\n", a, b);
        return 0;
    }

    c_viejo = a; // Inicialización para error

    do {
        // Fórmula de falsa posición
        c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        iteraciones++;

        // Chequear signo para acotar intervalo
        if (f(a) * f(c) > 0)
            a = c;
        else if (f(a) * f(c) < 0)
            b = c;
        else {
            // Encontramos raíz exacta
            error = 0;
            break;
        }

        // Calcular error porcentual estimado
        error = fabs((c - c_viejo) / c) * 100;
        c_viejo = c;

    } while (error > tolerancia);

    printf("\nLa raiz es: %.6lf\n", c);
    printf("Con un error de: %.6lf%%\n", error);
    printf("Iteraciones: %d\n", iteraciones);

    return 0;
}
