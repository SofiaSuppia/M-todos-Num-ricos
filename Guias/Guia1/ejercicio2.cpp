/*Obtenga las raíces de la siguiente función utilizando los métodos de bisección y de regla falsi:
g (a) = a^10 − 1 
1. Verifique la convergencia de ambos métodos para un error ε = 1 × 10−5. ¿Cuantas iteraciones
le tomó a cada uno de los métodos en encontrar las raíces con la precisión deseada?. Obtenga
una gráfica comparativa del error de aproximación de cada método en función del número de
iteraciones.
2. Cuantas iteraciones le tomó a cada uno de los métodos en encontrar las raíces con la precisión deseada?. Obtenga
una gráfica comparativa del error de aproximación de cada método en función del número de
iteraciones.*/

//Método bisección
// Método de régula falsi

#include <iostream>
#include <math.h>
#include <limits> // Necesario para limpiar el buffer de entrada

using namespace std;

// --- Definición de Funciones ---
double f(double x) {
    return pow(x, 10) - 1;
}

// Derivada aproximada de f(x)
double df(double x) {
    double h = 0.0001;
    return (f(x + h) - f(x)) / h;
}

// Para el método de Punto Fijo, necesitamos una función g(x) tal que x = g(x).
// Usaremos la del ejemplo: x = sqrt(sin(sqrt(x))) para la ecuación x^2 - sin(sqrt(x)) = 0
double g(double x) {
    /*if (x < 0 || sin(sqrt(x)) < 0) {
        return NAN; // Not a Number, para indicar un dominio inválido
    }*/
    return pow(x, 10) - 1;
}

// Derivada aproximada de g(x)
double gp(double x) {
    double h = 0.0001;
    return (g(x + h) - g(x)) / h;
}

void metodo_biseccion();
void metodo_falsa_posicion();
void metodo_punto_fijo();
void metodo_newton_raphson();
void metodo_secante();


// --- Menú Principal ---
int main() {
    int opcion;
    do {
        cout << "\n\n===== MENU DE METODOS NUMERICOS =====" << endl;
        cout << "  g (a) = a^10 − 1" << endl;
        cout << "=====================================" << endl;
        cout << "\n--- METODOS CERRADOS (requieren un intervalo) ---" << endl;
        cout << "1. Biseccion" << endl;
        cout << "2. Falsa Posicion" << endl;
        cout << "\n--- METODOS ABIERTOS (requieren puntos iniciales) ---" << endl;
        cout << "3. Punto Fijo (usa una funcion g(x) diferente)" << endl;
        cout << "4. Newton-Raphson" << endl;
        cout << "5. Secante" << endl;
        cout << "\n0. Salir" << endl;
        cout << "=====================================" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        // Limpiar buffer en caso de entrada no numérica
        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            opcion = -1; // Asignar opción inválida
        }

        switch (opcion) {
            case 1: metodo_biseccion(); break;
            case 2: metodo_falsa_posicion(); break;
            case 3: metodo_punto_fijo(); break;
            case 4: metodo_newton_raphson(); break;
            case 5: metodo_secante(); break;
            case 0: cout << "Saliendo del programa.\n"; break;
            default: cout << "Opcion no valida. Por favor, intente de nuevo.\n"; break;
        }

    } while (opcion != 0);
    return 0;
}

// --- Implementación de los Métodos ---

void metodo_biseccion() {
    cout << "\n--- Metodo de Biseccion ---\n";
    double a, b, tolerancia, c = 0, c_viejo, error = 100;
    int iteraciones = 0;

    cout << "Ingrese el inicio del intervalo (a): "; cin >> a;
    cout << "Ingrese el final del intervalo (b): "; cin >> b;
    cout << "Ingrese la tolerancia (%): "; cin >> tolerancia;

    if (f(a) * f(b) > 0) {
        cout << "No se puede asegurar una raiz en el intervalo dado.\n";
        return;
    }

    c_viejo = a;
    do {
        c = (a + b) / 2.0;
        iteraciones++;

        if (f(a) * f(c) > 0) a = c;
        else if (f(a) * f(c) < 0) b = c;
        else break; // Raíz exacta encontrada

        error = fabs((c - c_viejo) / c) * 100;
        c_viejo = c;

    } while (error > tolerancia && iteraciones < 100);

    cout << "\nRaiz aproximada: " << c << endl;
    cout << "Error final: " << error << "%" << endl;
    cout << "Iteraciones: " << iteraciones << endl;
}

void metodo_falsa_posicion() {
    cout << "\n--- Metodo de Falsa Posicion ---\n";
    double a, b, tolerancia, c = 0, c_viejo, error = 100;
    int iteraciones = 0;

    cout << "Ingrese el inicio del intervalo (a): "; cin >> a;
    cout << "Ingrese el final del intervalo (b): "; cin >> b;
    cout << "Ingrese la tolerancia (%): "; cin >> tolerancia;

    if (f(a) * f(b) > 0) {
        cout << "No se puede asegurar una raiz en el intervalo dado.\n";
        return;
    }

    c_viejo = a;
    do {
        c = (a * f(b) - b * f(a)) / (f(b) - f(a));
        iteraciones++;

        if (f(a) * f(c) > 0) a = c;
        else if (f(a) * f(c) < 0) b = c;
        else break; // Raíz exacta encontrada

        error = fabs((c - c_viejo) / c) * 100;
        c_viejo = c;

    } while (error > tolerancia && iteraciones < 100);

    cout << "\nRaiz aproximada: " << c << endl;
    cout << "Error final: " << error << "%" << endl;
    cout << "Iteraciones: " << iteraciones << endl;
}

void metodo_punto_fijo() {
    cout << "\n--- Metodo de Punto Fijo ---\n";
    cout << "Este metodo usa g(x) = sqrt(sin(sqrt(x)))\n";
    double x0, x1, tolerancia, error = 100;
    int i = 0;

    cout << "Ingrese el valor inicial (x0): "; cin >> x0;
    cout << "Ingrese la tolerancia (%): "; cin >> tolerancia;

    if (fabs(gp(x0)) >= 1) {
        cout << "El metodo podria no converger porque |g'(x0)| >= 1\n";
    }

    do {
        i++;
        x1 = g(x0);
        if (isnan(x1)) {
            cout << "Error: Se salio del dominio de g(x).\n";
            return;
        }
        error = fabs((x1 - x0) / x1) * 100;
        cout << "Iteracion " << i << ": x = " << x1 << ", error = " << error << "%\n";
        x0 = x1;
    } while (error > tolerancia && i < 100);

    cout << "\nRaiz aproximada: " << x1 << endl;
}

void metodo_newton_raphson() {
    cout << "\n--- Metodo de Newton-Raphson ---\n";
    double x0, x1, tolerancia, error = 100;
    int i = 0;

    cout << "Ingrese el valor inicial (x0): "; cin >> x0;
    cout << "Ingrese la tolerancia (error absoluto): "; cin >> tolerancia;

    do {
        i++;
        if (fabs(df(x0)) < 1e-6) {
            cout << "Derivada muy cercana a cero. El metodo falla.\n";
            return;
        }
        x1 = x0 - f(x0) / df(x0);
        error = fabs(x1 - x0);
        x0 = x1;
        cout << "Iteracion " << i << ": x = " << x1 << ", error = " << error << "\n";
    } while (error > tolerancia && i < 100);

    cout << "\nRaiz aproximada: " << x1 << endl;
}

void metodo_secante() {
    cout << "\n--- Metodo de la Secante ---\n";
    double x0, x1, x2, tolerancia, error = 100;
    int i = 0;

    cout << "Ingrese el primer valor (x0): "; cin >> x0;
    cout << "Ingrese el segundo valor (x1): "; cin >> x1;
    cout << "Ingrese la tolerancia (error absoluto): "; cin >> tolerancia;

    do {
        i++;
        if (fabs(f(x1) - f(x0)) < 1e-6) {
            cout << "Diferencia f(x1)-f(x0) muy pequena. El metodo falla.\n";
            return;
        }
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        error = fabs(x2 - x1);
        x0 = x1;
        x1 = x2;
        cout << "Iteracion " << i << ": x = " << x2 << ", error = " << error << "\n";
    } while (error > tolerancia && i < 100);

    cout << "\nRaiz aproximada: " << x2 << endl;
}

/*RESULTADOS
--- Metodo de Biseccion ---
Ingrese el inicio del intervalo (a): 0.5
Ingrese el final del intervalo (b): 1.5
Ingrese la tolerancia (%): 0.00001

Raiz aproximada: 1
Error final: 100%
Iteraciones: 1

--- Metodo de Falsa Posicion ---
Ingrese el inicio del intervalo (a): 0.5
Ingrese el final del intervalo (b): 1.5
Ingrese la tolerancia (%): 0.00001

Raiz aproximada: 0.999745
Error final: 0.00246276%
Iteraciones: 100
*/

