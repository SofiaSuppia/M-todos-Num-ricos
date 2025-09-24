#include <stdio.h>
#include <math.h>

/*
PROBLEMA DEL CONTENEDOR SIN TAPA
Hoja metálica: 10 × 16 cm
Volumen deseado: 100 cm³
Precisión: 1 × 10^-9 cm

Si recortamos cuadrados de lado x en cada esquina:
- Largo: 16 - 2x
- Ancho: 10 - 2x  
- Alto: x
- Volumen: V(x) = x(10-2x)(16-2x) = 100

Ecuación a resolver: f(x) = 4x³ - 52x² + 160x - 100 = 0
Derivada: f'(x) = 12x² - 104x + 160
*/

// Función del volumen menos 100
double f(double x) {
    return 4*x*x*x - 52*x*x + 160*x - 100;
}

// Derivada de f(x)
double df(double x) {
    return 12*x*x - 104*x + 160;
}

// Función para calcular las dimensiones del contenedor
void calcular_dimensiones(double x) {
    double largo = 16 - 2*x;
    double ancho = 10 - 2*x;
    double alto = x;
    double volumen = largo * ancho * alto;
    
    printf("\n=== DIMENSIONES DEL CONTENEDOR ===\n");
    printf("Lado del cuadrado recortado: %.10lf cm\n", x);
    printf("Largo del contenedor: %.6lf cm\n", largo);
    printf("Ancho del contenedor: %.6lf cm\n", ancho);
    printf("Alto del contenedor: %.6lf cm\n", alto);
    printf("Volumen calculado: %.6lf cm³\n", volumen);
    
    // Verificar si la solución es físicamente válida
    if (x > 0 && x < 5 && largo > 0 && ancho > 0) {
        printf("✅ SOLUCIÓN FÍSICAMENTE VÁLIDA\n");
    } else {
        printf("❌ SOLUCIÓN NO VÁLIDA FÍSICAMENTE\n");
        if (x <= 0) printf("   - El lado x debe ser positivo\n");
        if (x >= 5) printf("   - x debe ser menor que 5 cm (mitad del ancho)\n");
        if (largo <= 0) printf("   - El largo resultante es negativo o cero\n");
        if (ancho <= 0) printf("   - El ancho resultante es negativo o cero\n");
    }
}

// Método de Newton-Raphson para el problema del contenedor
double newton_raphson_contenedor(double x0, double tolerancia, int max_iter) {
    double x1, error;
    int i = 0;
    
    printf("\n=== MÉTODO DE NEWTON-RAPHSON ===\n");
    printf("Ecuación: f(x) = 4x³ - 52x² + 160x - 100 = 0\n");
    printf("Derivada: f'(x) = 12x² - 104x + 160\n");
    printf("Tolerancia: %.2e\n", tolerancia);
    printf("Valor inicial x0: %.6lf\n\n", x0);
    
    do {
        i++;
        
        // Verificar si la derivada es muy pequeña
        if (fabs(df(x0)) < 1e-12) {
            printf("Derivada muy pequeña. El método puede no converger.\n");
            return x0;
        }
        
        // Fórmula de Newton-Raphson
        x1 = x0 - f(x0) / df(x0);
        
        // Calcular error absoluto
        error = fabs(x1 - x0);
        
        printf("Iteración %d: x = %.12lf, f(x) = %.3e, error = %.3e\n", 
               i, x1, f(x1), error);
        
        // Actualizar x0 para la siguiente iteración
        x0 = x1;
        
    } while (error > tolerancia && i < max_iter);
    
    if (i >= max_iter) {
        printf("\n⚠️  Se alcanzó el número máximo de iteraciones.\n");
    } else {
        printf("\n✅ Convergencia alcanzada en %d iteraciones.\n", i);
    }
    
    return x1;
}

// Función para mostrar el menú
void mostrar_menu() {
    printf("\n=================================================\n");
    printf("           PROBLEMA DEL CONTENEDOR SIN TAPA\n");
    printf("=================================================\n");
    printf("Hoja metálica: 10 × 16 cm | Volumen objetivo: 100 cm³\n");
    printf("Precisión requerida: 1 × 10⁻⁹ cm\n");
    printf("=================================================\n");
    printf("1. Buscar solución cerca de x = 1 (contenedor alto)\n");
    printf("2. Buscar solución cerca de x = 2 (contenedor mediano)\n");
    printf("3. Buscar solución cerca de x = 4 (contenedor bajo)\n");
    printf("4. Buscar todas las soluciones automáticamente\n");
    printf("5. Introducir valor inicial personalizado\n");
    printf("6. Mostrar información del problema\n");
    printf("0. Salir del programa\n");
    printf("=================================================\n");
    printf("Seleccione una opción: ");
}

// Función para mostrar información del problema
void mostrar_info_problema() {
    printf("\n=================================================\n");
    printf("                    INFORMACIÓN DEL PROBLEMA\n");
    printf("=================================================\n");
    printf("\n📋 DESCRIPCIÓN:\n");
    printf("Se construye un contenedor sin tapa a partir de una hoja\n");
    printf("metálica rectangular de 10 × 16 cm recortando cuadrados\n");
    printf("de lado x en cada esquina.\n\n");
    
    printf("📐 GEOMETRÍA:\n");
    printf("• Largo del contenedor: 16 - 2x\n");
    printf("• Ancho del contenedor: 10 - 2x\n");
    printf("• Alto del contenedor: x\n");
    printf("• Volumen: V(x) = x(10-2x)(16-2x) = 100 cm³\n\n");
    
    printf("🔢 ECUACIÓN MATEMÁTICA:\n");
    printf("f(x) = 4x³ - 52x² + 160x - 100 = 0\n");
    printf("f'(x) = 12x² - 104x + 160\n\n");
    
    printf("⚠️  RESTRICCIONES FÍSICAS:\n");
    printf("• 0 < x < 5 cm (para que el contenedor tenga sentido)\n");
    printf("• x > 0 (no se puede recortar longitud negativa)\n");
    printf("• x < 5 (x debe ser menor que la mitad del ancho)\n\n");
    
    printf("🎯 OBJETIVO:\n");
    printf("Encontrar el valor de x con precisión 1 × 10⁻⁹ cm\n");
}

// Función para continuar
void continuar_programa() {
    printf("\nPresione Enter para continuar...");
    getchar();
    getchar(); // Para capturar el Enter pendiente
}

int main() {
    int opcion;
    double x0, resultado, tolerancia = 1e-9;
    int max_iter = 100;
    
    printf("=================================================\n");
    printf("    BIENVENIDO AL SOLUCIONADOR DEL PROBLEMA DEL CONTENEDOR\n");
    printf("=================================================\n");
    
    do {
        mostrar_menu();
        scanf("%d", &opcion);
        
        if (opcion == 0) {
            printf("\n¡Gracias por usar el programa!\n");
            break;
        }
        
        switch(opcion) {
            case 1: {
                printf("\nBuscando solución cerca de x = 1...\n");
                x0 = 1.0;
                resultado = newton_raphson_contenedor(x0, tolerancia, max_iter);
                calcular_dimensiones(resultado);
                continuar_programa();
                break;
            }
            
            case 2: {
                printf("\nBuscando solución cerca de x = 2...\n");
                x0 = 2.0;
                resultado = newton_raphson_contenedor(x0, tolerancia, max_iter);
                calcular_dimensiones(resultado);
                continuar_programa();
                break;
            }
            
            case 3: {
                printf("\nBuscando solución cerca de x = 4...\n");
                x0 = 4.0;
                resultado = newton_raphson_contenedor(x0, tolerancia, max_iter);
                calcular_dimensiones(resultado);
                continuar_programa();
                break;
            }
            
            case 4: {
                printf("\n===========================================\n");
                printf("        BÚSQUEDA AUTOMÁTICA DE TODAS LAS SOLUCIONES\n");
                printf("===========================================\n");
                
                double valores_iniciales[] = {0.5, 1.5, 2.5, 3.5, 4.5};
                int num_valores = 5;
                
                for (int i = 0; i < num_valores; i++) {
                    printf("\n--- Prueba %d: x0 = %.1lf ---\n", i+1, valores_iniciales[i]);
                    resultado = newton_raphson_contenedor(valores_iniciales[i], tolerancia, max_iter);
                    calcular_dimensiones(resultado);
                    printf("\n");
                }
                continuar_programa();
                break;
            }
            
            case 5: {
                printf("\nIngrese el valor inicial x0: ");
                scanf("%lf", &x0);
                
                if (x0 <= 0 || x0 >= 5) {
                    printf("⚠️  Advertencia: x0 = %.3lf está fuera del rango físico válido (0, 5)\n", x0);
                    printf("El método puede converger a una solución no física.\n");
                }
                
                resultado = newton_raphson_contenedor(x0, tolerancia, max_iter);
                calcular_dimensiones(resultado);
                continuar_programa();
                break;
            }
            
            case 6: {
                mostrar_info_problema();
                continuar_programa();
                break;
            }
            
            default:
                printf("\n❌ Opción inválida. Por favor seleccione una opción del 0 al 6.\n");
                continuar_programa();
                break;
        }
        
    } while (opcion != 0);
    
    return 0;
}