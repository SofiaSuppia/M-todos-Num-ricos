#include <iostream>
#include <cmath>

using namespace std;

/* Both methods are closed methods, meaning they require an interval [a, b] where the function changes sign.
Also, both methods require an initial guess for the root.
They use root location */

// Function prototypes
double math_function(double x);
int get_error_type();
double calculate_error(double c, double old_c, int error_type);
void update_interval(double *a, double *b, double c);
void print_results(double c, double error, int max_iterations);

/* f(x) = exp(-x) - x

Biseccion:
    Raiz = 0.567143 +- 0.000001
    Numero iteraciones = 23
Regula Falsi:
    Raiz = 0.567235 +- 0.000001
    Numero iteraciones = 785

Se puede notar que el metodo de biseccion converge mucho mas rapido que el de Regula Falsi ya que, en Regula Falsi, al trazar rectas en el intervalo [a, b], como la funcion 
es muy creciente en el dominio, esto hace que las rectas que van intersectando el dominio, sean muchas hasta llegar a la raiz, demostrando 
que este metodo converge mucho mas lento.
*/

int main(int argc, char const *argv[]) {
    // Variable definitions
    double a, b, tolerance, old_c, c, error;
    int max_iterations = 0;
    int method, error_type;
    
    // Request user to input values
    printf("Ingrese el valor de a: ");
    scanf("%lf", &a);
    
    printf("Ingrese el valor de b: ");
    scanf("%lf", &b);
    
    printf("Ingrese la tolerancia del error: ");
    scanf("%lf", &tolerance);
    
    // Check if the function has a root in the interval [a, b]
    if(math_function(a) * math_function(b) > 0) {
        printf("No se puede garantizar la existencia de una raíz en el intervalo [%lf, %lf]\n", a, b);
        exit(0);
    }
    old_c = a; // old_c will be equal to a at minus one iteration

    printf("Ingrese el tipo de metodo (1 para bisección, 2 para regula falsi): ");
    scanf("%d", &method);
    
    switch(method) {
        case 1:
            printf("Método de Bisección seleccionado\n");
            error_type = get_error_type();
            do {
                c = (a + b) / 2; // We divide the interval in half
                max_iterations++; // We count the number of iterations
                
                update_interval(&a, &b, c);
                error = calculate_error(c, old_c, error_type);
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            break;
        case 2:
            printf("Método de Regula Falsi seleccionado\n");
            error_type = get_error_type();
            do {
                // This formula is derived from the equation of a line
                c = (a * math_function(b) - b * math_function(a)) / (math_function(b) - math_function(a)); // We calculate the point c using the Regula Falsi formula
                max_iterations++; // We count the number of iterations
                
                update_interval(&a, &b, c);
                error = calculate_error(c, old_c, error_type);
                old_c = c; // We update old_c for the next iteration
            } while(error > tolerance);

            print_results(c, error, max_iterations);
            break;
        default:
            printf("Opción inválida. Seleccione 1 o 2.\n");
    }
    
    return 0;
}

// Mathematical function to find its root
double math_function(double x) {
    return exp(-x) - x;
}

// Function to get error type from user
int get_error_type() {
    int error_type;
    printf("¿Solicita error absoluto o porcentual? (1 para absoluto, 2 para porcentual): ");
    scanf("%d", &error_type);
    return error_type;
}

// Function to calculate error (absolute or percentage)
double calculate_error(double c, double old_c, int error_type) {
    if(error_type == 2) {
        return fabs((c - old_c) / c) * 100; // Percentage error
    } else if(error_type == 1) {
        return fabs(c - old_c); // Absolute error
    } else {
        printf("Tipo de error inválido. Seleccione 1 o 2.\n");
        exit(0);
    }
}

// Function to update interval based on sign analysis
void update_interval(double *a, double *b, double c) {
    if(math_function(*a) * math_function(c) > 0) {
        *a = c; // The root is in the interval [c, b]
    } else if(math_function(*a) * math_function(c) < 0) {
        *b = c; // The root is in the interval [a, c]
    } else {
        printf("La raíz exacta es (sin error): %lf\n", c);
        exit(1);
    }
}

// Function to print final results
void print_results(double c, double error, int max_iterations) {
    printf("La raíz aproximada es: %lf, con un error de %lf\n", c, error);
    printf("Número de iteraciones: %d\n", max_iterations);
}


