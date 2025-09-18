#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <iostream>

#define FIL 4 //cambiar para ingresar los "b"
#define COL 4

using namespace std;

int main(int argc, char const *argv[])
{
	
	int opcion=0,iteracion=0,i=0,j=0,p=0,l=0,m=0,k=0;
	double /*arreglo[FIL][COL]={0},*/dato=0,swap=0,b[FIL]={0},factor=0,det=1,x[FIL]={0},suma=0;
	
	/*printf("\n---CARGA DE LA MATRIZ---");
	for (i = 0; i < FIL; i++)
	{
		for (j = 0; j < COL; j++)
		{
			printf("\narreglo[%d][%d]= ",i,j);
			scanf("%lf",&dato);
			arreglo[i][j]=dato;
		}
		
	}
	
	for(i=0;i<FIL;i++){
		printf("\n");
		for(j=0;j<COL;j++){
			printf("%lf\t",arreglo[i][j]);
		}
	}
	*/

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
    fp = fopen("SistEcuaLineales\\matriz.dat","r");
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
		
		printf("\n----MENU----\n");
		printf("1- Eliminacion gaussiana\n");
		printf("2- Salir\n");
		printf("------------\n");
		printf("OPCION: ");
		scanf("%d",&opcion);
		
		switch (opcion)
		{
		case 1:
			for(j=0;j<FIL;j++){
				p=j;
				
				if(fabs(arreglo[j][j]) < 0.00001){
					l=j+1;
					do
					{
						if(fabs(arreglo[l][j]) > fabs(arreglo[j][j])){
							p=l;
						}
					} while (l < FIL);
					
					m=0;
					do
					{
						swap=arreglo[p][m];
						arreglo[p][m]=arreglo[j][m];
						arreglo[j][m]=swap;
					} while (m <= FIL);
					
					swap=b[p];
					b[p]=b[j];
					b[j]=swap;
				}
				for(i=j+1;i<FIL;i++){
					factor=(-1)*(arreglo[i][j]/arreglo[j][j]);
					
					for(k=j;k<COL;k++){
						arreglo[i][k]=arreglo[i][k]+factor*arreglo[j][k];
					}
					b[i]=b[i]+factor*b[j];
				}
			}//esta en forma triangular
			
			for(i=0;i<FIL;i++){
				det=det*arreglo[i][i];
			}
			
			
			if(det==0){
				printf("\nMatriz singular");
				return 1;
			}
			else{
				for(i=FIL-1;i>=0;i--){
					suma=b[i];
					for(k=i+1;k<FIL;k++){
						suma=suma-arreglo[i][k]*x[k];
					}
					x[i]=suma/arreglo[i][i];
				}
				
			}
			
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
				printf("\nx[%d]= %lf",i,x[i]);
			}
			
			printf("\nDeterminante: %lf",det);
			break;
		}
	} while (opcion != 2);
	return 0;
}