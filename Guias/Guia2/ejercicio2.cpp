/*
Realice un programa adecuado para evaluar el punto fijo de funciones arbitrarias.
Utilice este programa para aproximar los puntos fijos (si es que hay alguno) de:
1. g(x) = x^5 − 3x^3 − 2x^2 + 2
2. g(a) = cos(sin(a))
3. g(n) = n^n − cos(n)

x0 sugerido: 1
tolerancia: 0.0001
*/

#include <stdio.h>
#include <math.h>

// Definir las funciones g(x)
double g1(double x) { return pow(x,5) - 3*pow(x,3) - 2*pow(x,2) + 2; }
double g2(double x) { return cos(sin(x)); }
double g3(double x) { return pow(x, x) - cos(x); }

// Derivada aproximada g'(x)
double gp(double (*g)(double), double x) {
	double h = 0.0001;
	return (g(x + h) - g(x)) / h;
}

int main() {
	int opcion;
	double x0, x1, tolerancia, error;
	int i = 0;
	double (*g)(double) = NULL;

	printf("Seleccione la función para punto fijo:\n");
	printf("1. g(x) = x^5 - 3x^3 - 2x^2 + 2\n");
	printf("2. g(x) = cos(sin(x))\n");
	printf("3. g(x) = x^x - cos(x)\n");
	printf("Opción: ");
	scanf("%d", &opcion);

	switch(opcion) {
		case 1: g = g1; break;
		case 2: g = g2; break;
		case 3: g = g3; break;
		default: printf("Opción inválida\n"); return 1;
	}

	printf("Ingrese x0: ");
	scanf("%lf", &x0);
	printf("Ingrese tolerancia: ");
	scanf("%lf", &tolerancia);

	// Verificar criterio de convergencia
	if (fabs(gp(g, x0)) >= 1) {
		printf("El metodo no converge porque |g'(x0)| >= 1 (g'(x0) = %.4lf)\n", gp(g, x0));
		return 0;
	}

	do {
		i++;
		x1 = g(x0);
		error = fabs((x1 - x0) / x1) * 100;
		printf("Iteracion %d: x = %.8lf, error = %.8lf%%\n", i, x1, error);
		x0 = x1;
	} while (error > tolerancia);

	printf("\nPunto fijo aproximado: %.8lf\n", x1);
	printf("Error final: %.8lf%%\n", error);
	printf("Iteraciones: %d\n", i);
	return 0;
}


