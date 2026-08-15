#include <cstdio>
#include <math.h>

// Definir la función f(x)
double f(double x) {
    return log(x) + exp(sin(x)) - x; // Ejemplo de la hoja
}

int main() {
    double a, b, tolerancia, c_nuevo, c_viejo, error;
    int iteraciones = 0;

    printf("Ingrese a: ");
    scanf("%lf", &a);
    printf("Ingrese b: ");
    scanf("%lf", &b);
    printf("Ingrese tolerancia (%%): ");
    scanf("%lf", &tolerancia);

    // Verificar si hay raíz en el intervalo
    if (f(a) * f(b) > 0) {
        printf("No hay raiz en el intervalo [%.6f, %.6f]\n", a, b);
        return 0;
    }

    c_viejo = a; // Valor inicial para calcular error

    do {
        // Punto medio
        c_nuevo = (a + b) / 2.0;
        iteraciones++;

        // Determinar en qué subintervalo está la raíz
        if (f(a) * f(c_nuevo) > 0)
            a = c_nuevo;
        else if (f(a) * f(c_nuevo) < 0)
            b = c_nuevo;
        else {
            // Encontramos la raíz exacta
            error = 0;
            break;
        }

        // Calcular error porcentual estimado
        error = fabs((c_nuevo - c_viejo) / c_nuevo) * 100;
        c_viejo = c_nuevo;

    } while (error > tolerancia);

    printf("La raiz es: %f\n", c_nuevo);
    printf("Con un error de: %f%%\n", error);
    printf("Iteraciones: %d\n", iteraciones);

    return 0;
}
