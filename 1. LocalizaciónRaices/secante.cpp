#include <iostream>
#include <math.h>

// ============================================================
// METODO DE LA SECANTE (mejorado)
// x_{n+1} = x_n  -  f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
//
// Es como Newton-Raphson pero reemplaza la derivada f'(x) por el
// cociente incremental entre las dos ultimas aproximaciones.
// Trabaja con DOS puntos iniciales (x0 y x1) y NO usa derivada.
//
// Muestra en CADA iteracion: x0, x1, x2, f(x2) y el error
// absoluto |x2 - x1|.
// ============================================================

// Definir la función f(x)
double f(double x) {
    return log(x * x + 1) - sin(x); // Parcial 1 - 2024
}

int main() {
    double x0, x1, x2, tolerancia, error;
    int i = 0, max_iter = 1000;

    printf("Ingrese x0: ");
    scanf("%lf", &x0);
    printf("Ingrese x1: ");
    scanf("%lf", &x1);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    printf("\n iter |      x0      |      x1      |      x2      |    f(x2)     |     error     \n");
    printf("------+-------------+-------------+-------------+--------------+--------------+\n");

    do {
        i++;

        // Si f(x1) == f(x0) la secante es horizontal y no corta al eje x:
        // la formula divide por cero. Avisamos pero seguimos iterando.
        if (fabs(f(x1) - f(x0)) < 1e-12) {
            printf("Ojo: f(x1) y f(x0) muy parecidos (secante horizontal).\n");
        }

        // Fórmula de la secante: interseccion de la recta que pasa por
        // (x0, f(x0)) y (x1, f(x1)) con el eje x.
        x2 = x1 - (f(x1) * (x1 - x0) / (f(x1) - f(x0)));

        // Error absoluto estimado: que tan lejos se movio la aprox.
        error = fabs(x2 - x1);

        // Mostrar los valores de la iteracion actual
        printf("%5d | %11.8f | %11.8f | %11.8f | %12.8f | %12.8f\n",
               i, x0, x1, x2, f(x2), error);

        // Si f(x2) = 0 la raiz es exacta: terminar
        if (f(x2) == 0)
            break;

        // Actualizar valores: x2 pasa a ser el "nuevo", x1 el "viejo"
        x0 = x1;
        x1 = x2;

    } while (error > tolerancia && i < max_iter);

    printf("------+-------------+-------------+-------------+--------------+--------------+\n");
    printf("\nRaiz aproximada: %.10f\n", x2);
    printf("Error final: %.10f\n", error);
    printf("Iteraciones: %d\n", i);

    // Avisar si se corto por el limite de iteraciones (no convergio)
    if (i == max_iter && error > tolerancia)
        printf("Ojo: se alcanzo el maximo de iteraciones (%d).\n", max_iter);

    return 0;
}