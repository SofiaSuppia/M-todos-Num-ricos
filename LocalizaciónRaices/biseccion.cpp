#include <iostream>
#include <math.h>
using namespace std;

// Definir la función f(x)
double f(double x) {
    return 5 * x - 3; // Ejemplo de la hoja
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
