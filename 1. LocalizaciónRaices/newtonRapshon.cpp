#include <iostream>
#include <math.h>

// ============================================================
// METODO DE NEWTON-RAPHSON (mejorado)
// x_{n+1} = x_n - f(x_n) / f'(x_n)
//
// - La derivada es numerica (sirve para CUALQUIER f).
// - Muestra en CADA iteracion: x, f(x), f'(x) y el error
//   absoluto |x_nuevo - x_viejo|.
// - Tope de 1000 iteraciones para evitar ciclos infinitos.
// ============================================================

// Definir la función f(x)
double f(double x) {
    return log(x * x + 1) - sin(x); // Parcial 1 - 2024: f(x)=ln(x^2+1)-sin(x)
}

// Derivada aproximada (numerica) de f(x)
double df(double x) {
    double h = 0.001; // Paso pequeño
    return (f(x + h) - f(x)) / h;
}

int main() {
    double x0, x1, tolerancia, error;
    int i = 0, max_iter = 1000;

    printf("Ingrese x0: ");
    scanf("%lf", &x0);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    printf("\n iter |      x          |      f(x)      |     f'(x)      |     error     \n");
    printf("------+-----------------+----------------+----------------+---------------+\n");

    do {
        i++;

        // Si la derivada es (casi) cero la formula divide por cero:
        // la tangente es horizontal y x_nuevo se va "al infinito".
        // Avisamos pero seguimos (puede recuperarse en la prox. iteracion).
        if (fabs(df(x0)) < 1e-4) {
            printf("Ojo: |f'(x)| = %.10f muy pequena en x = %.10f\n",
                   df(x0), x0);
        }

        // Fórmula de Newton-Raphson
        x1 = x0 - f(x0) / df(x0);

        // Error absoluto estimado: que tan lejos se movio la aprox.
        error = fabs(x1 - x0);

        // Mostrar los valores de la iteracion actual
        printf("%5d | %17.10f | %16.10f | %16.10f | %15.10f\n",
               i, x1, f(x1), df(x0), error);

        // Si f(x1) = 0 la raiz es exacta: terminar
        if (f(x1) == 0)
            break;

        // Actualizar x0 para la siguiente iteración
        x0 = x1;

    } while (error > tolerancia && i < max_iter);

    printf("------+-----------------+----------------+----------------+---------------+\n");
    printf("\nRaiz aproximada: %.10f\n", x1);
    printf("Error final: %.10f\n", error);
    printf("Iteraciones: %d\n", i);

    // Avisar si se corto por el limite de iteraciones (no convergio)
    if (i == max_iter && error > tolerancia)
        printf("Ojo: se alcanzo el maximo de iteraciones (%d).\n", max_iter);

    return 0;
}