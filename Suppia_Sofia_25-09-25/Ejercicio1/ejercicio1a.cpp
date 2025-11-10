#include <iostream>
#include <math.h>
using namespace std;

// Definir la función f(x)
double f(double x) {
    return ((x+1)/(x+4))-0.25*x; 
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
        //error = fabs((c - c_viejo) / c) * 100;
        error = fabs((c - c_viejo) / (1/2 * (c + c_viejo))) * 100;
        c_viejo = c;

    } while (error > tolerancia);

    printf("La raiz es: %.15lf\n", c);
    printf("Con un error de: %.15lf%%\n", error);
    printf("Error en notacion cientifica: %.15e%%\n", error);
    printf("Iteraciones: %d\n", iteraciones);

    return 0;
}
