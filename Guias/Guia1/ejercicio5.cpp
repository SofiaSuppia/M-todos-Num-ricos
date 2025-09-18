#include <stdio.h>
#include <math.h>

//Método bisección
// Método de régula falsi

// Definir la función eficiencia
double eficiencia(double x) {
	// x = T2/T1
	double gamma = 5.0/3.0;
	double num = log(x) - (1.0 - 1.0/x);
	double den = log(x) + (1.0 + 1.0/(gamma * (x - 1.0)));
	return num / den;
}

// Función para la raíz: f(x) = eficiencia(x) - 0.3
double f(double x) {
	return eficiencia(x) - 0.3;
}

int main() {
	double a = 1.01, b = 10.0, tol = 1e-6, c, error;
	int max_iter = 100, iter = 0;

	printf("Método de bisección para encontrar T2/T1 con eficiencia 0.3\n");
	printf("Intervalo inicial sugerido: a = 1.01, b = 10\n");
	printf("Ingrese tolerancia: "); //0.000001
	scanf("%lf", &tol);
	printf("Ingrese máximo de iteraciones: "); //100
	scanf("%d", &max_iter);

	if (f(a) * f(b) > 0) {
		printf("No hay cambio de signo en el intervalo inicial.\n");
		return 1;
	}

	do {
		c = (a + b) / 2.0;
		error = fabs(f(c));
		iter++;
		if (f(a) * f(c) < 0)
			b = c;
		else
			a = c;
	} while (error > tol && iter < max_iter);

	printf("\nEl valor de T2/T1 para eficiencia 0.3 es: %.8lf\n", c);
	printf("Error final: %.8lf\n", error);
	printf("Iteraciones: %d\n", iter);
	return 0;
}
