#include <stdio.h>
#include <math.h>

// Recurrencia para raíz cuadrada (n=2)
double raiz_cuadrada(double A, double p0, double tol, int max_iter) {
	double pk = p0, pk1;
	int k = 0;
	do {
		pk1 = 0.5 * (pk + A / pk);
		printf("Iter %d: p = %.10lf\n", k+1, pk1);
		if (fabs(pk1 - pk) < tol) break;
		pk = pk1;
		k++;
	} while (k < max_iter);
	return pk1;
}

// Recurrencia para raíz cúbica (n=3)
double raiz_cubica(double A, double p0, double tol, int max_iter) {
	double pk = p0, pk1;
	int k = 0;
	do {
		pk1 = (2.0*pk + A/(pk*pk)) / 3.0;
		printf("Iter %d: p = %.10lf\n", k+1, pk1);
		if (fabs(pk1 - pk) < tol) break;
		pk = pk1;
		k++;
	} while (k < max_iter);
	return pk1;
}

int main() {
	int opcion, max_iter = 100;
	double A, p0, tol;

	printf("Calculo de raiz cuadrada y cubica usando Newton-Raphson\n");
	printf("1. Raiz cuadrada\n2. Raiz cubica\nOpcion: ");
	scanf("%d", &opcion);
	printf("Ingrese el valor de A: ");
	scanf("%lf", &A);
	printf("Ingrese la aproximacion inicial p0: ");
	scanf("%lf", &p0);
	printf("Ingrese la tolerancia: ");
	scanf("%lf", &tol);

	double resultado = 0;
	if (opcion == 1) {
		resultado = raiz_cuadrada(A, p0, tol, max_iter);
		printf("\nAproximacion de sqrt(%.4lf) = %.10lf\n", A, resultado);
	} else if (opcion == 2) {
		resultado = raiz_cubica(A, p0, tol, max_iter);
		printf("\nAproximacion de cbrt(%.4lf) = %.10lf\n", A, resultado);
	} else {
		printf("Opcion invalida\n");
	}
	return 0;
}
