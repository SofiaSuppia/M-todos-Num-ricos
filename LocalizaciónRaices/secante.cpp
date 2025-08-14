#include <iostream>
#include <math.h>

// Definir la función f(x)
double f(double x) {
    return pow(x, 3) - 9 * x + 3; // Ejemplo, cámbialo según el problema
}

int main() {
    double x0, x1, x2, tolerancia, error;
    int i = 0;

    printf("Ingrese x0: ");
    scanf("%lf", &x0);
    printf("Ingrese x1: ");
    scanf("%lf", &x1);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    do {
        i++;

        // Evitar división por cero
        if (fabs(f(x1) - f(x0)) < 1e-12) {
            printf("Diferencia muy pequena entre f(x1) y f(x0). El metodo no puede continuar.\n");
            return 0;
        }

        // Fórmula de la secante
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));

        // Calcular error absoluto
        error = fabs(x2 - x1);

        // Actualizar valores
        x0 = x1;
        x1 = x2;

        printf("Iteracion %d: x = %.6lf, error = %.6lf, f(x) = %.6lf\n",
               i, x2, error, f(x2));

    } while (error > tolerancia && i < 10000);

    printf("\nRaiz aproximada: %.6lf\n", x2);
    printf("Error final: %.6lf\n", error);
    printf("Iteraciones: %d\n", i);

    return 0;
}
