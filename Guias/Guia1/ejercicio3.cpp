#include <stdio.h>
#include <math.h>

//Método bisección
// Método de régula falsi

// Constantes dadas
const double g = 9.81; // m/s^2
const double c = 14.0; // kg/s
const double t = 7.0;  // s
const double v_objetivo = 35.0; // m/s

// Función f(m):
double f(double m) {
	return (g * m / c) * (1 - exp(-c * t / m)) - v_objetivo;
}

// Método de bisección
double biseccion(double a, double b, double tol, int max_iter) {
	double fa = f(a);
	double fb = f(b);
	if (fa * fb > 0) {
		printf("No hay cambio de signo en el intervalo inicial.\n");
		return NAN;
	}
	double c, fc;
	for (int i = 0; i < max_iter; ++i) {
		c = (a + b) / 2.0;
		fc = f(c);
		if (fabs(fc) < tol || (b - a) / 2.0 < tol)
			return c;
		if (fa * fc < 0) {
			b = c;
			fb = fc;
		} else {
			a = c;
			fa = fc;
		}
	}
	return c;
}

// Método de regula falsi (falsa posición)
double regula_falsi(double a, double b, double tol, int max_iter) {
	double fa = f(a);
	double fb = f(b);
	if (fa * fb > 0) {
		printf("No hay cambio de signo en el intervalo inicial.\n");
		return NAN;
	}
	double c, fc;
	for (int i = 0; i < max_iter; ++i) {
		c = b - fb * (b - a) / (fb - fa);
		fc = f(c);
		if (fabs(fc) < tol)
			return c;
		if (fa * fc < 0) {
			b = c;
			fb = fc;
		} else {
			a = c;
			fa = fc;
		}
	}
	return c;
}

int main() {
	int opcion;
	double a = 40.0, b = 90.0, tol = 1e-6;
	int max_iter = 100;

	printf("\n--- Cálculo de masa de paracaidista ---\n");
	printf("Sugerencia: masa mínima = 40 kg, masa máxima = 90 kg\n");
	printf("1. Método de Bisección\n");
	printf("2. Método de Regula Falsi\n");
	printf("3. Ambos métodos\n");
	printf("Seleccione una opción: ");
	scanf("%d", &opcion);

	printf("Ingrese masa mínima a (kg) [sugerido 40]: ");
	scanf("%lf", &a);
	printf("Ingrese masa máxima b (kg) [sugerido 90]: ");
	scanf("%lf", &b);
	printf("Ingrese tolerancia [sugerido 0.000001]: ");
	scanf("%lf", &tol);
	printf("Ingrese máximo de iteraciones [sugerido 100]: ");
	scanf("%d", &max_iter);

	if (opcion == 1) {
		double m_bis = biseccion(a, b, tol, max_iter);
		printf("Masa encontrada por bisección: %.8lf kg\n", m_bis);
	} else if (opcion == 2) {
		double m_rf = regula_falsi(a, b, tol, max_iter);
		printf("Masa encontrada por regula falsi: %.8lf kg\n", m_rf);
	} else if (opcion == 3) {
		double m_bis = biseccion(a, b, tol, max_iter);
		double m_rf = regula_falsi(a, b, tol, max_iter);
		printf("Masa encontrada por bisección: %.8lf kg\n", m_bis);
		printf("Masa encontrada por regula falsi: %.8lf kg\n", m_rf);
	} else {
		printf("Opción no válida.\n");
	}
	return 0;
}
