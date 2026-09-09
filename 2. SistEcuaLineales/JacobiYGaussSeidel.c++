#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <iostream>

using namespace std;

#define FIL 2 //cambiar para ingresar los "b"
#define COL 2

// Programa principal que resuelve un sistema de ecuaciones lineales
// mediante los métodos iterativos de Jacobi y Gauss-Seidel.
int main(int argc, char const *argv[])
{
    
    int i=0,j=0,opcion=0,iteracion=0;
    double b[FIL]={0},x0[FIL]={0},xv[FIL]={0},e=1,tol=0,suma=0;

    // Archivo que contiene la matriz del sistema lineal
    FILE *fp;
    char ch;

    // Se abre la matriz desde un archivo de texto
    fp = fopen("SistEcuaLineales/matriz.dat","r");
    if ( fp == NULL )
    {
        perror( "No se puede abrir el archivo");
		return 1;
    }

    // Se cuenta cuántas filas tiene la matriz para dimensionar el arreglo
    int filas=0;
    char c;
    int maxValues = 0;
	int columnas;
    
    while((c = fgetc(fp)) != EOF)
    {
        if(c == '\n'){
            filas++;
		}
    }

    printf("numero de filas=%i\n",filas);

    // Se vuelve a abrir el archivo para reiniciar el puntero de lectura
    // y luego guardar los datos en un arreglo bidimensional
    fclose(fp);
    fp = fopen("SistEcuaLineales/matriz.dat","r");
    double arreglo[filas][100]={0};


    // Se cargan los números leídos del archivo en la matriz 'arreglo'
    
	for(i = 0; i < filas; i++) {
		j = 0;
		do {
        fscanf(fp, "%lf", &(arreglo[i][j]));
        j++;
			} while((c = fgetc(fp)) != '\n');
	}

	columnas = j;
	
	printf("numero de columnas= %i\n",columnas);

    // Se imprime la matriz para verificar que la lectura fue correcta
    
        for(i=0;i<filas;i++)
    {
        for(j=0;j<columnas;j++)
        {
            printf("%lf ",arreglo[i][j]);
        }
        printf("\n");
    }

    // Se solicitan los valores del vector b del sistema Ax = b
    printf("\n---Carga de los b---");
        for(i= 0; i < FIL; i++){
            printf("\nb[%d]= ",i);
            scanf("%lf",&b[i]);
        }

    // Menú principal para elegir el método iterativo a aplicar
    do
    {
        
        printf("\n---MENU---");
        printf("\n1- Jacobi");
        printf("\n2- Gauss - Seidel");
        printf("\n----------");
        printf("\nOPCION: ");
        scanf("%d",&opcion);

        switch(opcion){

            case 1:

            // Método iterativo de Jacobi
            printf("\nIngrese tolerancia: ");
            scanf("%lf",&tol);

            // Verifica si la matriz puede converger: diagonal dominante y diagonal no nula
            for(i=0;i<FIL;i++){

                suma=0;
                for(j=0;j<COL;j++){

                    if(j!=i){
                        suma=suma+fabs(arreglo[i][j]);
                    }
                }
                if(fabs(arreglo[i][i]) < suma){
                    printf("\nLa matriz no es diagonalmente dominante");
                }
                if(fabs(arreglo[i][i])==0){
                    printf("\nCeros en la diagonal");
                    return 1;
                }
            }//verificacion de posible convergencia
            
            iteracion=0;

            // Repetición de las iteraciones del método de Jacobi hasta alcanzar la tolerancia
            do
            {
                iteracion++;
                double error=0;
                for(i=0;i<FIL;i++){
                    suma=0;
                    for(j=0;j<COL;j++){
                        if(j != i){
                            suma=suma+arreglo[i][j]*xv[j];
                        }
                    }
                    // Fórmula del método de Jacobi para calcular la nueva aproximación
                    x0[i]=(b[i]-suma)/arreglo[i][i];
                    error=error+pow((x0[i]-xv[i]),2);
                    xv[i]=x0[i];
                }//calculo de xnuevo (x0)
                e=sqrt(error);
                printf("\niter %d",iteracion);
                printf("\nerror: %.12lf",e);
                printf("\n");
            } while (e > tol && iteracion < 10000 );
            
            // Muestra la solución obtenida por Jacobi
            for(i=0;i<FIL;i++){

                printf("\nX0[%d]: %lf",i+1,x0[i]);

            }
            printf("\nError: %.12lf",e);
            printf("\nIteraciones: %d",iteracion);





            break;

            case 2:
            
            // Método iterativo de Gauss-Seidel
            printf("\nIngrese tolerancia: ");
            scanf("%lf",&tol);

            // Verifica la convergencia de la matriz antes de iterar
            for(i=0;i<FIL;i++){

                suma=0;
                for(j=0;j<COL;j++){

                    if(j!=i){
                        suma=suma+fabs(arreglo[i][j]);
                    }
                }
                if(fabs(arreglo[i][i]) < suma){
                    printf("\nLa matriz no es diagonalmente dominante");
                }
                if(fabs(arreglo[i][i])==0){
                    printf("\nCeros en la diagonal");
                    return 1;
                }
            }//verificacion de posible convergencia
            
            iteracion=0;

            // Repetición de las iteraciones del método de Gauss-Seidel
            do
            {
                iteracion++;
                double error=0;
                for(i=0;i<FIL;i++){
                    suma=0;
                if(i==0){
                    // Se usan valores actuales de las demás incógnitas para la primera fila
                    for(j=1;j<COL;j++){
                        suma=suma+arreglo[i][j]*xv[j];
                    }
                    x0[i]=(b[i]-suma)/(arreglo[i][i]);
                }else{
                    // Se aprovechan los valores ya actualizados de x0 para acelerar la convergencia
                    for(j=0;j<i;j++){
                        suma=suma+arreglo[i][j]*x0[j];
                    }
                    for(j=i+1;j<COL;j++){
                        suma=suma+arreglo[i][j]*xv[j];
                    }

                    x0[i]=(b[i]-suma)/(arreglo[i][i]);
                
                    error=error+pow((x0[i]-xv[i]),2);
                    xv[i]=x0[i];
                }//calculo de xnuevo (x0)
                }
                e=sqrt(error);
                printf("\niter %d",iteracion);
                printf("\nerror: %.12lf",e);
                printf("\n");
            } while (e > tol && iteracion < 10000 );
            
            // Muestra la matriz y la solución final obtenida por Gauss-Seidel
            printf("\n--Matriz resuelta--\n");
					
				for(i=0;i<filas;i++)
			{
				for(j=0;j<columnas;j++)
				{
					printf("%lf ",arreglo[i][j]);
				}
				printf("\n");
			} 

            for(i=0;i<FIL;i++){

                printf("\nX0[%d]: %lf",i+1,x0[i]);

            }
            printf("\nError: %.12lf",e);
            printf("\nIteraciones: %d",iteracion);


            break;

        }




    } while (opcion != 3);
    

    // Fin del programa
    return 0;
}