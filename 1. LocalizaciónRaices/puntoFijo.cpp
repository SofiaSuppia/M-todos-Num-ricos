#include <iostream>
#include <math.h>

// Definir g(x)
double g(double x) {
    return sqrt(sin(sqrt(x))); // Ejemplo: cambiar según problema
}

// Derivada aproximada g'(x)
double gp(double x) {
    double h = 0.001; // Paso pequeño
    return (g(x + h) - g(x)) / h;
}

int main() {
    double x0, x1, tolerancia, error;
    int i = 0;

    printf("Ingrese x0: ");
    scanf("%lf", &x0);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    // Verificar criterio de convergencia
    if (fabs(gp(x0)) >= 1) {
        printf("El metodo no converge porque |g'(x0)| >= 1\n");
        return 0;
    }

    do {
        i++;
        x1 = g(x0);
        error = fabs((x1 - x0) / x1) * 100; // Error porcentual estimado
        printf("Iteracion %d: x = %.6lf, error = %.6lf%%\n", i, x1, error);
        x0 = x1;
    } while (error > tolerancia);

    printf("\nRaiz aproximada: %.6lf\n", x1);
    printf("Error final: %.6lf%%\n", error);
    printf("Iteraciones: %d\n", i);

    return 0;
}
