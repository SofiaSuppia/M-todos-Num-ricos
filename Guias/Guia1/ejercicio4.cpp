#include <stdio.h>
#include <math.h>

//Método bisección
// Método de régula falsi

// Definir la función f(x)
double f(double x) {
	return 5 * x - 3; // Ejemplo de la hoja
}

int main() {
	double a, b, delta, c, c_viejo, error;
	int iteraciones = 0;

	printf("Ingrese a: "); //0
	scanf("%lf", &a);
	printf("Ingrese b: ");//2
	scanf("%lf", &b);
	printf("Ingrese tolerancia absoluta (delta): "); //0.0001
	scanf("%lf", &delta);

	// Calcular el número teórico de iteraciones necesarias
	int n_teorico = (int)ceil(log2((b - a) / delta));
	printf("\nNúmero mínimo de iteraciones teóricas para delta = %.10lf: %d\n\n", delta, n_teorico);

	// Verificar si hay raíz en el intervalo
	if (f(a) * f(b) > 0) {
		printf("No hay raiz en el intervalo [%.6lf, %.6lf]\n", a, b);
		return 0;
	}

	c_viejo = a; // Valor inicial para calcular error

	do {
		// Punto medio
		c = (a + b) / 2.0;
		iteraciones++;

		// Determinar en qué subintervalo está la raíz
		if (f(a) * f(c) > 0)
			a = c;
		else if (f(a) * f(c) < 0)
			b = c;
		else {
			// Encontramos la raíz exacta
			error = 0;
			break;
		}

		// Calcular error absoluto
		error = fabs((b - a) / 2.0);
		c_viejo = c;

	} while (error > delta);

	printf("La raiz es: %.10lf\n", c);
	printf("Con un error absoluto de: %.10lf\n", error);
	printf("Iteraciones realizadas: %d\n", iteraciones);

	return 0;
}
