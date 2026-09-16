#include <cstdio>
#include <math.h>

// ============================================================
// METODO DE BISECCION (mejorado)
// Muestra en CADA iteracion: a, b, f(a), c, f(c), el error
// absoluto estimado y el error relativo porcentual.
// ============================================================

// Definir la función f(x)
double f(double x) {
    // f(x) = ln(x^2 + 1) - sen(x)
    return log(x * x + 1) - sin(x);
}

int main() {
    double a, b, tolerancia, c_nuevo, c_viejo;
    double error_abs, error_rel;
    int iteraciones = 0, max_iter = 1000;

    printf("Ingrese a: ");
    scanf("%lf", &a);
    printf("Ingrese b: ");
    scanf("%lf", &b);
    printf("Ingrese tolerancia (%%): ");
    scanf("%lf", &tolerancia);

    // Verificar si hay raiz en el intervalo
    if (f(a) * f(b) > 0) {
        printf("No hay raiz en el intervalo [%.6f, %.6f]\n", a, b);
        return 0;
    }

    c_viejo = a; // Valor inicial para calcular el error relativo

    printf("\n iter |      a      |      b      |    f(a)     |      c      |    f(c)     |  error abs  |  error rel %%\n");
    printf("------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+\n");

    do {
        // Punto medio del intervalo (aproximacion de la raiz)
        c_nuevo = (a + b) / 2.0;
        iteraciones++;

        // Error absoluto estimado: la raiz esta a menos de la
        // mitad de la longitud del intervalo respecto de c
        error_abs = (b - a) / 2.0;

        // Error relativo porcentual entre la aprox. actual y la anterior
        error_rel = fabs((c_nuevo - c_viejo) / c_nuevo) * 100;

        // Mostrar los valores de la iteracion actual
        printf("%5d | %11.8f | %11.8f | %11.8f | %11.8f | %11.8f | %11.8f | %11.8f\n",
               iteraciones, a, b, f(a), c_nuevo, f(c_nuevo), error_abs, error_rel);

        // Si f(c) = 0, la raiz es exacta: terminar
        if (f(c_nuevo) == 0) {
            error_rel = 0;
            break;
        }

        // Determinar en qué subintervalo está la raiz
        if (f(a) * f(c_nuevo) > 0)
            a = c_nuevo; // la raiz esta en [c_nuevo, b]
        else
            b = c_nuevo; // la raiz esta en [a, c_nuevo]

        c_viejo = c_nuevo; // guardar aprox. anterior para el error relativo

    } while (error_rel > tolerancia && iteraciones < max_iter);

    printf("------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+\n");
    printf("\nLa raiz es: %.8f\n", c_nuevo);
    printf("Con un error relativo de: %.8f%%\n", error_rel);
    printf("Error absoluto estimado: %.8f\n", error_abs);
    printf("Iteraciones: %d\n", iteraciones);

    // Avisar si se corto por el limite de iteraciones (no convergio)
    if (iteraciones == max_iter && error_rel > tolerancia)
        printf("Ojo: se alcanzo el maximo de iteraciones (%d).\n", max_iter);

    return 0;
}