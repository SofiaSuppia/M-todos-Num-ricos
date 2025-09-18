#include <iostream>
#include <math.h>
using namespace std;

#define GRADO 10 //grado del polinomio
#define COL 2 //columnas de la matriz de datos
#define FIL 10 //filas de la matriz de datos

double f(double x){
    return x*x; //funcion a interpolar
}
 
int main(int argc, char const *argv[])
{
    int i=0, j=0, k=0, opcion=0;
    double x[GRADO]={0};
    double y[GRADO]={0};
    double xsombrero = 0, prod=0, sum=0, error=0, arreglo[FIL][COL]={0}, b[FIL]={0};
    
    do
    {
        printf("\n---MENU---");
        printf("\n1. Lagrange");
        printf("\n2. Salir");
        printf("\nOPCION: ");
        scanf("%d",&opcion);

        switch(opcion)
        {
            case 1:
                printf("\n---LAGRANGE---");
                printf("\nIngrese los datos de x:\n");
                for (size_t i = 0; i < GRADO; i++)
                {
                    printf("x[%d]: ",i+1);
                    scanf("%lf",&x[i]);
                }

                printf("\n Calcular los datos de y\n");
                for (size_t i = 0; i < GRADO; i++)
                {
                    printf("\ny[%d]: %lf", i+1, f(x[i]));
                    scanf("%lf",&y[i]);
                }
                
                printf("\nIngrese el valor a interpolar: ");
                scanf("%lf",&xsombrero);
                for (size_t k = 0; k < GRADO; k++)
                {
                    prod=1;
                    for (size_t j = 0; j < GRADO; j++)
                    {
                        if (j != k)
                        {
                            prod = prod * ((xsombrero - x[j]) / (x[k] - x[j]));
                        }
                    }
                    sum = sum + y[k] * prod;
                }
                printf("\nEl valor interpolado es: %lf", sum);
                error = fabs(sum - f(xsombrero));
                printf("\nEl error absoluto es: %lf", error);
                printf("\nEl error relativo es: %lf", error/f(xsombrero));
                
                break;

            case 2:
                printf("\nSaliendo...");
                break;

            default:
                printf("\nOpción no válida");
                break;
        }

    } while (opcion != 2);
    
    return 0;
}
