#include <iostream>
#include <math.h>

// ============================================================
// METODO DE PUNTO FIJO (mejorado)
// Resuelve f(x) = ln(x^2 + 1) - sin(x) = 0
//
// Despeje:  f(x)=0  ->  x^2 + 1 = e^{sin(x)}  ->  x = sqrt(e^{sin(x)} - 1)
// Entonces g(x) = sqrt(e^{sin(x)} - 1).
//
// Criterio de convergencia: si |g'(x)| < 1 en un entorno de la
// raiz, la iteracion  x_{k+1} = g(x_k)  converge.
// Muestra en CADA iteracion: xviejo, g(xviejo), |g'(xviejo)| y error.
// ============================================================

// g(x) definida para el problema (raiz en [1, 2])
double g(double x) {
    return sqrt(exp(sin(x)) - 1);
}

// Derivada aproximada (numerica), sirve para CUALQUIER g(x)
double gp(double x) {
    double h = 0.001; // Paso pequeño
    return (g(x + h) - g(x)) / h;
}

int main() {
    double xviejo, xnuevo, tolerancia, error;
    int i = 0, max_iter = 1000;

    printf("Ingrese xviejo: ");
    scanf("%lf", &xviejo);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);

    // Criterio de convergencia: |g'(x)| < 1 alrededor de la raiz
    // es condicion suficiente. Si en xviejo da >= 1, avisamos pero
    // seguimos iterando: puede converger igual (no es condicion necesaria).
    if (fabs(gp(xviejo)) >= 1) {
        printf("Ojo: |g'(xviejo)| = %.6lf >= 1 (no garantiza convergencia).\n",
               fabs(gp(xviejo)));
        printf("Igualmente itero hasta el maximo permitido...\n");
    } else {
        printf("|g'(xviejo)| = %.6lf < 1  ->  convergencia garantizada\n",
               fabs(gp(xviejo)));
    }

    printf("\n iter |    x viejo   |  g(x) = x nuevo |  |g'(xviejo)|  |   error abs \n");
    printf("------+--------------+-----------------+---------------+--------------+\n");

    do {
        i++;
        // Aproximacion nueva: x_{k+1} = g(x_k)
        xnuevo = g(xviejo);

        // Error absoluto estimado: cambio entre dos aproximaciones
        error = fabs(xnuevo - xviejo);

        // Mostrar los valores de la iteracion actual
        printf("%5d | %12.8f | %15.8f | %13.8f | %12.8f\n",
               i, xviejo, xnuevo, fabs(gp(xviejo)), error);

        // Actualizar para la siguiente iteracion
        xviejo = xnuevo;

    } while (error > tolerancia && i < max_iter);

    printf("------+--------------+-----------------+---------------+--------------+\n");
    printf("\nRaiz aproximada: %.8f\n", xnuevo);
    printf("Error final: %.8f\n", error);
    printf("Iteraciones: %d\n", i);

    // Avisar si se corto por el limite de iteraciones (no convergio)
    if (i == max_iter && error > tolerancia)
        printf("Ojo: se alcanzo el maximo de iteraciones (%d).\n", max_iter);

    return 0;
}