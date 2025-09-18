/*
Sea g(x) = x^2 + x − 4.
¿Podemos utilizar iteración de punto fijo para hallar las soluciones de la ecuación x = g(x)?

Respuesta:
No, no es recomendable usar iteración de punto fijo con g(x) = x^2 + x − 4, porque 
para que el método converja debe cumplirse que |g'(x)| < 1 en el entorno del punto fijo. 
En este caso, g'(x) = 2x + 1, que para valores reales puede ser mayor que 1 en valor 
absoluto, por lo que el método diverge.

Por ejemplo, si intentamos iterar desde x0 = 1:
g'(1) = 2*1 + 1 = 3 > 1
Esto indica que la iteración diverge.

Ventaja de tener g'(P) ≈ 0:
Si en un punto fijo P se cumple que g'(P) ≈ 0, la convergencia del método de punto 
fijo es muy rápida (superlineal), ya que el error en cada iteración se reduce mucho 
más que si g'(P) está cerca de 1.

Ejemplo de código para probar la convergencia:
*/

#include <stdio.h>
#include <math.h>

// Definir g(x)
double g(double x) {
	return pow(x, 2) + x - 4;
}

// Derivada aproximada g'(x)
double gp(double x) {
	double h = 0.0001;
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
		printf("El metodo NO converge porque |g'(x0)| >= 1 (g'(x0) = %.4lf)\n", gp(x0));
		return 0;
	}

	do {
		i++;
		x1 = g(x0);
		error = fabs((x1 - x0) / x1) * 100;
		printf("Iteracion %d: x = %.6lf, error = %.6lf%%\n", i, x1, error);
		x0 = x1;
	} while (error > tolerancia);

	printf("\nRaiz aproximada: %.6lf\n", x1);
	printf("Error final: %.6lf%%\n", error);
	printf("Iteraciones: %d\n", i);
	return 0;
}


