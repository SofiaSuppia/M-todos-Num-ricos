#include <iostream>
#include <math.h>

// Definir g(x)
double g(double x) {
    return exp(-x) ; // g(x) = e^{-e^{\sin(x)}}
}

// Derivada aproximada g'(x)
double gp(double x) {
    double h = 0.001; // Paso pequeño
    return (g(x + h) - g(x)) / h;
}

int main() {
    double xviejo, xnuevo, tolerancia, error;
    int i = 0;

    printf("Ingrese xviejo: ");
    scanf("%lf", &xviejo);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    // Verificar criterio de convergencia
    if (fabs(gp(xviejo)) >= 1) {
        printf("El metodo no converge porque |g'(xviejo)| >= 1\n");
        return 0;
    }

    do {
        i++;
        xnuevo = g(xviejo);
        error = fabs((xnuevo - xviejo)); // Error absoluto estimado   error = %.6lf%%
        printf("Iteracion %d: x = %.6lf, error = %.6lf\n", i, xnuevo, error);
        xviejo = xnuevo;
    } while (tolerancia < error);

    printf("\nRaiz aproximada: %.6lf\n", xnuevo);
    printf("Error final: %.6lf%%\n", error);
    printf("Iteraciones: %d\n", i);

    return 0;
}
