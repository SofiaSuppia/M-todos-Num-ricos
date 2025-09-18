#include <iostream>
#include <math.h>

// Definir la función f(x)
double f(double x) {
    return pow(x, 3) - 9 * x + 3; // Ejemplo, cámbialo según el problema
}

// Derivada aproximada de f(x)
double df(double x) {
    double h = 0.001; // Paso pequeño
    return (f(x + h) - f(x)) / h;
}

int main() {
    double x0, x1, tolerancia, error;
    int i = 0;

    printf("Ingrese x0: ");
    scanf("%lf", &x0);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    do {
        i++;

        // Verificar si la derivada es muy pequeña
        if (fabs(df(x0)) < 1e-4) {
            printf("Derivada muy pequena. El metodo puede no converger.\n");
            return 0;
        }

        // Fórmula de Newton-Raphson
        x1 = x0 - f(x0) / df(x0);

        // Calcular error absoluto
        error = fabs(x1 - x0);

        // Actualizar x0 para la siguiente iteración
        x0 = x1;

        printf("Iteracion %d: x = %.6lf, error = %.6lf, f(x) = %.6lf\n",
               i, x1, error, f(x1));

    } while (error > tolerancia && i < 10000);

    printf("\nRaiz aproximada: %.6lf\n", x1);
    printf("Error final: %.6lf\n", error);
    printf("Iteraciones: %d\n", i);

    return 0;
}
