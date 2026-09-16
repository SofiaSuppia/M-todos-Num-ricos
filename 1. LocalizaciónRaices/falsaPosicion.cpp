#include <iostream>
#include <math.h>

// ============================================================
// METODO DE REGULA FALSI (mejorado)
// Muestra en CADA iteracion: a, b, f(a), c, f(c), el error
// absoluto estimado y el error relativo porcentual.
// ============================================================

// Definir la función f(x)
double f(double x) {
    // f(x) = ln(x^2 + 1) - sen(x)
    return log(x * x + 1) - sin(x);
}

int main() {
    double a, b, tolerancia, c, c_viejo;
    double error_abs, error_rel;
    int iteraciones = 0, max_iter = 1000;

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

    c_viejo = a; // Inicialización para calcular el error relativo

    printf("\n iter |      a      |      b      |    f(a)     |      c      |    f(c)     |  error abs  |  error rel %%\n");
    printf("------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+\n");

    do {
        // Fórmula de falsa posición (interseccion de la
        // recta por (a,f(a)) y (b,f(b)) con el eje x)
        c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        iteraciones++;

        // Error absoluto estimado: cambio entre la aprox.
        // actual y la anterior
        error_abs = fabs(c - c_viejo);

        // Error relativo porcentual entre la aprox. actual y la anterior
        error_rel = fabs((c - c_viejo) / c) * 100;

        // Mostrar los valores de la iteracion actual
        printf("%5d | %11.8lf | %11.8lf | %11.8lf | %11.8lf | %11.8lf | %11.8lf | %11.8lf\n",
               iteraciones, a, b, f(a), c, f(c), error_abs, error_rel);

        // Si f(c) = 0, la raiz es exacta: terminar
        if (f(c) == 0) {
            error_rel = 0;
            break;
        }

        // Chequear signo para acotar intervalo
        if (f(a) * f(c) > 0)
            a = c; // la raiz esta en [c, b]
        else
            b = c; // la raiz esta en [a, c]

        c_viejo = c; // guardar aprox. anterior para el error

    } while (error_rel > tolerancia && iteraciones < max_iter);

    printf("------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+\n");
    printf("\nLa raiz es: %.8lf\n", c);
    printf("Con un error relativo de: %.8lf%%\n", error_rel);
    printf("Error absoluto estimado: %.8lf\n", error_abs);
    printf("Iteraciones: %d\n", iteraciones);

    // Avisar si se corto por el limite de iteraciones (no convergio)
    if (iteraciones == max_iter && error_rel > tolerancia)
        printf("Ojo: se alcanzo el maximo de iteraciones (%d).\n", max_iter);

    return 0;
}