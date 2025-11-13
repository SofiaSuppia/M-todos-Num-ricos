/*
 * PROBLEMA 1: Regresión No Lineal por Mínimos Cuadrados
 * 
 * Ajustar los datos a la función: f(x) = a + bx + ce^x
 * 
 * Datos:
 * x: 0.2, 0.6, 0.8, 1.1, 1.5, 1.6, 2.1, 2.3, 2.5, 3.0
 * y: 2.71, 3.92, 4.79, 6.22, 8.48, 9.15, 13.40, 15.58, 18.17, 27.08
 * 
 * El método de mínimos cuadrados minimiza: E = Σ(y_k - f(x_k))²
 * 
 * a) Obtener el sistema de ecuaciones lineales para a, b, c
 * b) Resolver el sistema usando método iterativo
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_DATOS 10
#define MAX_ITER 1000
#define TOLERANCIA 1e-6

// Datos del problema
double x_datos[N_DATOS] = {0.2, 0.6, 0.8, 1.1, 1.5, 1.6, 2.1, 2.3, 2.5, 3.0};
double y_datos[N_DATOS] = {2.71, 3.92, 4.79, 6.22, 8.48, 9.15, 13.40, 15.58, 18.17, 27.08};

// Función f(x) = a + bx + ce^x
double funcion_ajuste(double x, double a, double b, double c) {
    return a + b*x + c*exp(x);
}

// Calcular el error cuadrático medio
double calcular_error(double a, double b, double c) {
    double suma = 0.0;
    for (int k = 0; k < N_DATOS; k++) {
        double diff = y_datos[k] - funcion_ajuste(x_datos[k], a, b, c);
        suma += diff * diff;
    }
    return suma;
}

/*
 * PARTE A: DERIVACIÓN DEL SISTEMA DE ECUACIONES
 * 
 * Para minimizar E = Σ(y_k - (a + bx_k + ce^(x_k)))²
 * 
 * Derivamos respecto a cada parámetro e igualamos a cero:
 * 
 * ∂E/∂a = -2Σ(y_k - a - bx_k - ce^(x_k)) = 0
 * ∂E/∂b = -2Σ(y_k - a - bx_k - ce^(x_k))·x_k = 0
 * ∂E/∂c = -2Σ(y_k - a - bx_k - ce^(x_k))·e^(x_k) = 0
 * 
 * Reorganizando, obtenemos el sistema lineal 3x3:
 * 
 * [Σ1      Σx_k      Σe^(x_k)    ] [a]   [Σy_k        ]
 * [Σx_k    Σx_k²     Σx_k·e^(x_k)] [b] = [Σx_k·y_k    ]
 * [Σe^(x_k) Σx_k·e^(x_k) Σe^(2x_k)] [c]   [Σe^(x_k)·y_k]
 */

void construir_sistema(double A[3][3], double b[3]) {
    // Inicializar a cero
    for (int i = 0; i < 3; i++) {
        b[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            A[i][j] = 0.0;
        }
    }
    
    // Calcular las sumatorias
    for (int k = 0; k < N_DATOS; k++) {
        double x = x_datos[k];
        double y = y_datos[k];
        double ex = exp(x);
        
        // Fila 1: derivada respecto a 'a'
        A[0][0] += 1.0;           // Σ1
        A[0][1] += x;             // Σx_k
        A[0][2] += ex;            // Σe^(x_k)
        b[0] += y;                // Σy_k
        
        // Fila 2: derivada respecto a 'b'
        A[1][0] += x;             // Σx_k
        A[1][1] += x * x;         // Σx_k²
        A[1][2] += x * ex;        // Σx_k·e^(x_k)
        b[1] += x * y;            // Σx_k·y_k
        
        // Fila 3: derivada respecto a 'c'
        A[2][0] += ex;            // Σe^(x_k)
        A[2][1] += x * ex;        // Σx_k·e^(x_k)
        A[2][2] += ex * ex;       // Σe^(2x_k)
        b[2] += ex * y;           // Σe^(x_k)·y_k
    }
}

void imprimir_sistema(double A[3][3], double b[3]) {
    printf("\n=== SISTEMA DE ECUACIONES LINEALES 3x3 ===\n\n");
    printf("Matriz de coeficientes A y vector b:\n\n");
    
    char vars[] = {'a', 'b', 'c'};
    
    for (int i = 0; i < 3; i++) {
        printf("  ");
        for (int j = 0; j < 3; j++) {
            printf("%12.6f ", A[i][j]);
        }
        printf("| %12.6f\n", b[i]);
    }
    
    printf("\nEn forma de ecuaciones:\n");
    for (int i = 0; i < 3; i++) {
        printf("  %.6f·a + %.6f·b + %.6f·c = %.6f\n", 
               A[i][0], A[i][1], A[i][2], b[i]);
    }
}

/*
 * PARTE B: RESOLUCIÓN DEL SISTEMA USANDO GAUSS-SEIDEL
 * 
 * Método iterativo que actualiza cada variable usando los valores más recientes:
 * 
 * a^(k+1) = (b[0] - A[0][1]*b^(k) - A[0][2]*c^(k)) / A[0][0]
 * b^(k+1) = (b[1] - A[1][0]*a^(k+1) - A[1][2]*c^(k)) / A[1][1]
 * c^(k+1) = (b[2] - A[2][0]*a^(k+1) - A[2][1]*b^(k+1)) / A[2][2]
 */

int gauss_seidel(double A[3][3], double b[3], double x[3], double tol, int max_iter) {
    double x_old[3];
    int iter = 0;
    double error;
    
    // Valores iniciales (aproximación inicial)
    x[0] = 0.0;  // a inicial
    x[1] = 0.0;  // b inicial
    x[2] = 0.0;  // c inicial
    
    printf("\n=== MÉTODO DE GAUSS-SEIDEL ===\n");
    printf("\nValores iniciales: a = %.6f, b = %.6f, c = %.6f\n", x[0], x[1], x[2]);
    printf("Tolerancia: %e\n", tol);
    printf("Máximo de iteraciones: %d\n\n", max_iter);
    
    printf("%-5s %-15s %-15s %-15s %-15s\n", "Iter", "a", "b", "c", "Error");
    printf("--------------------------------------------------------------------------------\n");
    
    do {
        // Guardar valores anteriores
        for (int i = 0; i < 3; i++) {
            x_old[i] = x[i];
        }
        
        // Actualizar a (usando b y c anteriores)
        x[0] = (b[0] - A[0][1]*x[1] - A[0][2]*x[2]) / A[0][0];
        
        // Actualizar b (usando a nuevo y c anterior)
        x[1] = (b[1] - A[1][0]*x[0] - A[1][2]*x[2]) / A[1][1];
        
        // Actualizar c (usando a y b nuevos)
        x[2] = (b[2] - A[2][0]*x[0] - A[2][1]*x[1]) / A[2][2];
        
        // Calcular error relativo máximo
        error = 0.0;
        for (int i = 0; i < 3; i++) {
            double err_i = fabs((x[i] - x_old[i]) / x[i]);
            if (err_i > error) {
                error = err_i;
            }
        }
        
        iter++;
        
        if (iter <= 10 || iter % 10 == 0 || error < tol) {
            printf("%-5d %-15.8f %-15.8f %-15.8f %-15.8e\n", 
                   iter, x[0], x[1], x[2], error);
        }
        
    } while (error > tol && iter < max_iter);
    
    if (error <= tol) {
        printf("\nConvergencia alcanzada en %d iteraciones.\n", iter);
        return iter;
    } else {
        printf("\nNo se alcanzó la convergencia en %d iteraciones.\n", max_iter);
        return -1;
    }
}

/*
 * Método de Jacobi (alternativo)
 * Actualiza todas las variables simultáneamente usando solo valores de la iteración anterior
 */
int jacobi(double A[3][3], double b[3], double x[3], double tol, int max_iter) {
    double x_old[3];
    int iter = 0;
    double error;
    
    // Valores iniciales
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    
    printf("\n=== MÉTODO DE JACOBI ===\n");
    printf("\nValores iniciales: a = %.6f, b = %.6f, c = %.6f\n", x[0], x[1], x[2]);
    printf("Tolerancia: %e\n", tol);
    printf("Máximo de iteraciones: %d\n\n", max_iter);
    
    printf("%-5s %-15s %-15s %-15s %-15s\n", "Iter", "a", "b", "c", "Error");
    printf("--------------------------------------------------------------------------------\n");
    
    do {
        // Guardar valores anteriores
        for (int i = 0; i < 3; i++) {
            x_old[i] = x[i];
        }
        
        // Actualizar todas las variables simultáneamente
        double a_new = (b[0] - A[0][1]*x_old[1] - A[0][2]*x_old[2]) / A[0][0];
        double b_new = (b[1] - A[1][0]*x_old[0] - A[1][2]*x_old[2]) / A[1][1];
        double c_new = (b[2] - A[2][0]*x_old[0] - A[2][1]*x_old[1]) / A[2][2];
        
        x[0] = a_new;
        x[1] = b_new;
        x[2] = c_new;
        
        // Calcular error
        error = 0.0;
        for (int i = 0; i < 3; i++) {
            double err_i = fabs((x[i] - x_old[i]) / x[i]);
            if (err_i > error) {
                error = err_i;
            }
        }
        
        iter++;
        
        if (iter <= 10 || iter % 10 == 0 || error < tol) {
            printf("%-5d %-15.8f %-15.8f %-15.8f %-15.8e\n", 
                   iter, x[0], x[1], x[2], error);
        }
        
    } while (error > tol && iter < max_iter);
    
    if (error <= tol) {
        printf("\nConvergencia alcanzada en %d iteraciones.\n", iter);
        return iter;
    } else {
        printf("\nNo se alcanzó la convergencia en %d iteraciones.\n", max_iter);
        return -1;
    }
}

void mostrar_resultados(double a, double b, double c) {
    printf("\n========================================\n");
    printf("RESULTADOS FINALES\n");
    printf("========================================\n\n");
    
    printf("Función ajustada: f(x) = %.8f + %.8f·x + %.8f·e^x\n\n", a, b, c);
    
    // Calcular error cuadrático
    double E = calcular_error(a, b, c);
    printf("Error cuadrático total E = %.10f\n", E);
    printf("Error cuadrático medio = %.10f\n\n", E / N_DATOS);
    
    // Mostrar tabla de comparación
    printf("Comparación entre datos originales y valores ajustados:\n\n");
    printf("%-5s %-10s %-15s %-15s %-15s\n", "k", "x_k", "y_k (dato)", "f(x_k) (ajuste)", "Residuo");
    printf("--------------------------------------------------------------------------------\n");
    
    for (int k = 0; k < N_DATOS; k++) {
        double x = x_datos[k];
        double y_real = y_datos[k];
        double y_ajuste = funcion_ajuste(x, a, b, c);
        double residuo = y_real - y_ajuste;
        
        printf("%-5d %-10.2f %-15.6f %-15.6f %-15.6f\n", 
               k+1, x, y_real, y_ajuste, residuo);
    }
}

int main() {
    double A[3][3];
    double b_vec[3];
    double solucion[3];
    
    printf("╔════════════════════════════════════════════════════════════════╗\n");
    printf("║   PARCIAL 2 - 2017: EJERCICIO 1                               ║\n");
    printf("║   Regresión No Lineal por Mínimos Cuadrados                   ║\n");
    printf("║   f(x) = a + bx + ce^x                                        ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    
    // PARTE A: Construir el sistema de ecuaciones lineales
    printf("\n--- PARTE A: SISTEMA DE ECUACIONES LINEALES ---\n");
    construir_sistema(A, b_vec);
    imprimir_sistema(A, b_vec);
    
    // PARTE B: Resolver el sistema
    printf("\n\n--- PARTE B: RESOLUCIÓN DEL SISTEMA ---\n");
    printf("\nSeleccione el método a usar:\n");
    printf("1. Gauss-Seidel (Recomendado)\n");
    printf("2. Jacobi\n");
    printf("Opción: ");
    
    int opcion;
    if (scanf("%d", &opcion) != 1) {
        printf("Error al leer la opción.\n");
        return 1;
    }
    
    int iteraciones;
    if (opcion == 1) {
        iteraciones = gauss_seidel(A, b_vec, solucion, TOLERANCIA, MAX_ITER);
    } else if (opcion == 2) {
        iteraciones = jacobi(A, b_vec, solucion, TOLERANCIA, MAX_ITER);
    } else {
        printf("Opción no válida.\n");
        return 1;
    }
    
    if (iteraciones > 0) {
        // Mostrar resultados finales
        mostrar_resultados(solucion[0], solucion[1], solucion[2]);
        
        // Justificación de la elección del método
        printf("\n========================================\n");
        printf("JUSTIFICACIÓN DEL MÉTODO\n");
        printf("========================================\n\n");
        
        if (opcion == 1) {
            printf("Método elegido: GAUSS-SEIDEL\n\n");
            printf("Ventajas:\n");
            printf("  • Converge más rápido que Jacobi en la mayoría de los casos\n");
            printf("  • Usa los valores más actualizados en cada iteración\n");
            printf("  • Requiere menos iteraciones para alcanzar la tolerancia\n");
            printf("  • Mejor eficiencia computacional\n\n");
            printf("Este método es apropiado porque la matriz del sistema es\n");
            printf("diagonalmente dominante o simétrica positiva definida.\n");
        } else {
            printf("Método elegido: JACOBI\n\n");
            printf("Ventajas:\n");
            printf("  • Más simple de implementar\n");
            printf("  • Fácilmente paralelizable\n");
            printf("  • Garantiza convergencia si la matriz es diagonalmente dominante\n\n");
            printf("Desventajas:\n");
            printf("  • Generalmente más lento que Gauss-Seidel\n");
            printf("  • Requiere más iteraciones para converger\n");
        }
        
        printf("\nNúmero de iteraciones realizadas: %d\n", iteraciones);
        printf("Error deseado alcanzado: %e\n", TOLERANCIA);
    }
    
    return 0;
}
