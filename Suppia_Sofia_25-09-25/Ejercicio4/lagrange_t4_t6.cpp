#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PUNTOS 50

// Función para leer datos del archivo
int leer_puntos_datos(const char* nombre_archivo, double X[], double Y[], int* n) {
    FILE *p_archivo;
    
    p_archivo = fopen(nombre_archivo, "r");
    if (p_archivo == NULL) {
        printf("Error: No se puede abrir el archivo '%s'\n", nombre_archivo);
        return 0;
    }
    
    printf("Archivo '%s' abierto exitosamente\n", nombre_archivo);
    
    // Leer número de puntos de datos
    if (fscanf(p_archivo, "%d", n) != 1) {
        printf("Error: No se pudo leer el numero de puntos de datos\n");
        fclose(p_archivo);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("Error: Numero de puntos invalido (%d). Maximo permitido: %d\n", *n, MAX_PUNTOS);
        fclose(p_archivo);
        return 0;
    }
    
    // Leer puntos de datos
    for (int i = 0; i < *n; i++) {
        if (fscanf(p_archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: No se pudo leer el punto de dato %d\n", i + 1);
            fclose(p_archivo);
            return 0;
        }
    }
    
    fclose(p_archivo);
    printf("Se leyeron %d puntos de datos exitosamente\n\n", *n);
    return 1;
}

// Función para mostrar los puntos de datos
void imprimir_puntos_datos(double X[], double Y[], int n) {
    printf("Puntos de Datos:\n");
    printf("=================\n");
    printf("     i  |    t (tiempo)    |   T (temperatura)   \n");
    printf("--------|------------------|--------------------\n");
    for (int i = 0; i < n; i++) {
        printf("%6d  | %14.6f | %18.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

// Función para mostrar el polinomio de Lagrange en forma simbólica
void mostrar_polinomio_lagrange(double X[], double Y[], int n) {
    printf("\n======================================================\n");
    printf("          POLINOMIO DE LAGRANGE SIMBÓLICO\n");
    printf("======================================================\n");
    printf("P(t) = Σ T_i * L_i(t)  donde i = 0, 1, 2, ..., %d\n\n", n-1);
    
    // Mostrar cada término del polinomio
    printf("Desarrollo completo:\n");
    printf("P(t) = ");
    
    for(int k = 0; k < n; k++) {
        if(k > 0) {
            if(Y[k] >= 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
        } else if(Y[k] < 0) {
            printf("-");
        }
        
        printf("%.1f * L_%d(t)", fabs(Y[k]), k);
    }
    printf("\n\n");
    
    // Mostrar cada polinomio de Lagrange Li(t)
    printf("Donde los polinomios de Lagrange son:\n\n");
    
    for(int k = 0; k < n; k++) {
        printf("L_%d(t) = ", k);
        
        // Construir el producto para Lk(t)
        int primer_termino = 1;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                if(!primer_termino) printf(" * ");
                
                printf("(t - %.0f)", X[i]);
                primer_termino = 0;
            }
        }
        
        printf(" / ");
        
        // Construir el denominador
        primer_termino = 1;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                if(!primer_termino) printf(" * ");
                
                double denominador = X[k] - X[i];
                printf("(%.0f - %.0f)", X[k], X[i]);
                primer_termino = 0;
            }
        }
        
        // Calcular y mostrar el denominador numérico
        double denominador_numerico = 1.0;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                denominador_numerico *= (X[k] - X[i]);
            }
        }
        printf(" = ... / %.0f\n", denominador_numerico);
    }
    
    printf("\n======================================================\n");
    printf("          POLINOMIO EXPANDIDO (forma desarrollada)\n");
    printf("======================================================\n");
    
    // Para 5 puntos, mostrar algunos términos expandidos como ejemplo
    printf("El polinomio completo expandido sería de grado %d:\n", n-1);
    printf("P(t) = a₀ + a₁*t + a₂*t² + a₃*t³ + a₄*t⁴\n\n");
    
    printf("Sustituyendo los valores:\n");
    for(int k = 0; k < n; k++) {
        if(k > 0) {
            if(Y[k] >= 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
        } else if(Y[k] < 0) {
            printf("-");
        }
        
        printf("%.1f * [", fabs(Y[k]));
        
        // Mostrar el producto expandido parcialmente
        int contador = 0;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                if(contador > 0) printf(" * ");
                printf("(t-%.0f)", X[i]);
                contador++;
            }
        }
        
        // Calcular denominador
        double denominador = 1.0;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                denominador *= (X[k] - X[i]);
            }
        }
        printf(" / %.0f]", denominador);
    }
    
    printf("\n\n");
    printf("Nota: El polinomio completo expandido requiere algebra simbólica.\n");
    printf("      Para evaluaciones numéricas, use la forma de Lagrange mostrada arriba.\n");
}
double interpolacion_lagrange(double x, double X[], double Y[], int n) {
    double suma = 0.0;
    
    printf("Coeficientes de Lagrange para t = %.1f:\n", x);
    printf("---------------------------------------\n");
    
    for(int k = 0; k < n; k++) {        
        double producto = 1.0;
        for(int i = 0; i < n; i++) {    
            if(i != k) {
                // Lk(x) = ∏ (x - xi) / (xk - xi)
                producto *= (x - X[i]) / (X[k] - X[i]); 
            }
        }
        // Mostrar cada coeficiente de Lagrange
        printf("L%d(%.1f) = %.8f\n", k, x, producto);
        
        // P(x) = Σ Yi * Li(x)
        suma += Y[k] * producto;
    }
    
    return suma;
}

int main() {
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];
    int n;

    // Leer puntos de datos desde el archivo
    if (!leer_puntos_datos("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Suppia_Sofia_25-09-25\\Ejercicio4\\data.dat", X, Y, &n)) {
        printf("Error al leer los datos del archivo. Saliendo.\n");
        return 1;
    }
    
    // Mostrar los puntos de datos leídos
    imprimir_puntos_datos(X, Y, n);

    printf("======================================================\n");
    printf("     INTERPOLACIÓN DE LAGRANGE PARA t=4 y t=6       \n");
    printf("======================================================\n");
    printf("Datos: tiempo (t) vs temperatura (T)\n\n");
    
    // Puntos solicitados: t=4 y t=6
    double puntos[] = {4.0, 6.0};
    int num_puntos = 2;
    
    for (int p = 0; p < num_puntos; p++) {
        double t = puntos[p];
        
        printf("************************************************\n");
        printf("           EVALUANDO EN t = %.1f\n", t);
        printf("************************************************\n");
        
        // Realizar interpolación de Lagrange
        double temperatura = interpolacion_lagrange(t, X, Y, n);
        
        printf("\n*** RESULTADO ***\n");
        printf("Para t = %.1f\n", t);
        printf("Temperatura interpolada: T = %.6f\n", temperatura);
        
        // Verificar si es interpolación o extrapolación
        double t_min = X[0], t_max = X[0];
        for (int i = 1; i < n; i++) {
            if (X[i] < t_min) t_min = X[i];
            if (X[i] > t_max) t_max = X[i];
        }
        
        if (t >= t_min && t <= t_max) {
            printf("Nota: Este es un caso de INTERPOLACIÓN (t está dentro del rango de datos)\n");
        } else {
            printf("Nota: Este es un caso de EXTRAPOLACIÓN (t está fuera del rango de datos)\n");
        }
        
        printf("\nNota: Para datos experimentales no existe función teórica exacta,\n");
        printf("      por tanto no se puede calcular error teórico.\n");
        printf("\n\n");
    }
    
    printf("========================================\n");
    printf("    INTERPOLACIÓN COMPLETADA\n");
    printf("========================================\n");
    printf("\nResultados finales:\n");
    printf("T(4.0) = %.6f\n", interpolacion_lagrange(4.0, X, Y, n));
    printf("T(6.0) = %.6f\n", interpolacion_lagrange(6.0, X, Y, n));

    return 0;
}
