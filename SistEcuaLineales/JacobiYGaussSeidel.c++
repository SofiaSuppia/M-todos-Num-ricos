#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <iostream>

using namespace std;

#define FIL 2 //cambiar para ingresar los "b"
#define COL 2

int main(int argc, char const *argv[])
{
    
    int i=0,j=0,opcion=0,iteracion=0;
    double b[FIL]={0},x0[FIL]={0},xv[FIL]={0},e=1,tol=0,suma=0;

    FILE *fp;
    char ch;
    fp = fopen("SistEcuaLineales/matriz.dat","r");
    if ( fp == NULL )
    {
        perror( "No se puede abrir el archivo");
		return 1;
    }

    //contador de filas
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

    //cerramos y reabrimos el archivo para resetear el puntero
    fclose(fp);
    fp = fopen("SistEcuaLineales/matriz.dat","r");
    double arreglo[filas][100]={0};


    //Cargo los datos leidos en el array
    
	for(i = 0; i < filas; i++) {
		j = 0;
		do {
        fscanf(fp, "%lf", &(arreglo[i][j]));
        j++;
			} while((c = fgetc(fp)) != '\n');
	}

	columnas = j;
	
	printf("numero de columnas= %i\n",columnas);
    //imprimo la matriz para verificar que lo leyo correctamente
    
        for(i=0;i<filas;i++)
    {
        for(j=0;j<columnas;j++)
        {
            printf("%lf ",arreglo[i][j]);
        }
        printf("\n");
    }

    printf("\n---Carga de los b---");
        for(i= 0; i < FIL; i++){
            printf("\nb[%d]= ",i);
            scanf("%lf",&b[i]);
        }

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

            
            printf("\nIngrese tolerancia: ");
            scanf("%lf",&tol);

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
                    x0[i]=(b[i]-suma)/arreglo[i][i];
                    error=error+pow((x0[i]-xv[i]),2);
                    xv[i]=x0[i];
                }//calculo de xnuevo (x0)
                e=sqrt(error);
                printf("\niter %d",iteracion);
                printf("\nerror: %.12lf",e);
                printf("\n");
            } while (e > tol && iteracion < 10000 );
            
            for(i=0;i<FIL;i++){

                printf("\nX0[%d]: %lf",i+1,x0[i]);

            }
            printf("\nError: %.12lf",e);
            printf("\nIteraciones: %d",iteracion);





            break;

            case 2:
            
            printf("\nIngrese tolerancia: ");
            scanf("%lf",&tol);

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

            do
            {
                iteracion++;
                double error=0;
                for(i=0;i<FIL;i++){
                    suma=0;
                if(i==0){
                    for(j=1;j<COL;j++){
                        suma=suma+arreglo[i][j]*xv[j];
                    }
                    x0[i]=(b[i]-suma)/(arreglo[i][i]);
                }else{
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
    


    return 0;
}