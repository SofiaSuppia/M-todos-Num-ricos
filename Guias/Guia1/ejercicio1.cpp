/*Determine las raíces reales de:
f (x) = −2 + 7x − 5x^2 + 6x^3

1. Un método de utilidad y que muchas veces sirve
como guía para aproximar el valor de las raíces de una ecuación determinada, es graficar la función para visualizar las raíces. Al graficar
la función se logra acotar el dominio de búsqueda. Utilizando algún graficador de su agrado,
obtenga una estimación de las raíces de f (x) de manera gráfica.
2. Utilizando el método de bisección para encontrar la raíz más pequeña de f (x). Adopte
como valores iniciales xl = 0 y xu = 1, realice el proceso iterativo hasta que el error de la
aproximación se encuentre por debajo de 1 × 10−4.
*/
#include <iostream>
#include <math.h>
using namespace std;

// Definir la función f(x)
double f(double x) {
    return -2.0 + 7.0 * x - 5.0 * pow(x, 2) + 6.0 * pow(x, 3);
}

int main() {
    double a, b, tolerancia, c, c_viejo, error;
    int iteraciones = 0;

    cout << "Ingrese a: ";
    cin >> a;
    cout << "Ingrese b: ";
    cin >> b;
    cout << "Ingrese tolerancia (%): ";
    cin >> tolerancia;

    // Verificar si hay raíz en el intervalo
    if (f(a) * f(b) > 0) {
        cout << "No hay raiz en el intervalo [" << a << ", " << b << "]\n";
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

        // Calcular error porcentual estimado
        error = fabs((c - c_viejo) / c) * 100;
        c_viejo = c;

    } while (error > tolerancia);

    cout << "La raiz es: " << c << endl;
    cout << "Con un error de: " << error << "%" << endl;
    cout << "Iteraciones: " << iteraciones << endl;

    return 0;
}

/*Resultado
Ingrese a: 0
Ingrese b: 1
Ingrese tolerancia (%): 0.001
La raiz es: 0.333334
Con un error de: 0.000572203%
Iteraciones: 19*/