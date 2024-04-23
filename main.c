#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void getDataF(float ** x , float ** y , int dim);
float * allocFloat1D(float * tab , int dim);
void display(float * x , float * y , int dim);
float distance(float *  x,int indice );
float distance_1(float * x , int indice );
float lamda(float  distance,float distance_1);
float ro(float  lamda );
float psy(float  distance,float distance_1 ,int indice , float * y);
float ** matrice(float * x , float * y,float ** B ,int dim);
void writeOnFile(float ** A , float * B , int dim);
float * Cholesky();

int main(){
	float * x = NULL;
	float * y = NULL;
	int dim = 7;
	float * B = NULL;
	float ** A = NULL;

	getDataF(&x , &y , dim);
	A = matrice(x , y , &B , dim);
	display(x, y , dim);
	writeOnFile(A , B , dim);
	Cholesky();
	return 0;
}

void writeOnFile(float ** A , float * B , int dim){
	FILE * fp = fopen("matrice.txt" , "w");
	if(fp){
		fprintf(fp , "%d\n" , dim-2);
		for(int i = 1 ; i < dim-1 ; i++){
			for(int j = 1 ; j < dim-1 ; j++){
				fprintf(fp , "%f " , A[i][j]);
			}
			fprintf(fp , "\n");
		}
		for(int i = 1 ; i < dim-1 ; i++){
			fprintf(fp , "%f\n" , B[i]);
		}
	}else{
		printf("Impossible d'ouvrir ce fichier !! \n");
	}
	fclose(fp);
}

float ** allocFloat2D(float ** tab , int dim){
	tab = (float **)calloc(dim , sizeof(float *));
	for(int i = 0 ; i < dim ; i++){
		tab[i]=(float *)calloc(dim , sizeof(float));
	}
	return tab;
}

float ** matrice(float * x , float * y,float ** B ,int dim){
	float ** A = allocFloat2D(A , dim);
	float * b = allocFloat1D(b, dim);

	for(int i = 1 ; i < dim ; i++){
		A[i-1][i-1] = 2;
		A[i-1][i] = lamda(distance(x,i) , distance_1(x , i));
		A[i][i-1] = ro(A[i-1][i]);
		b[i-1] = psy(distance(x , i) , distance_1(x , i) , i , y);
	}
	A[dim-1][dim-1] = 2;

	*B = b;

	return A;
}

float psy(float distance , float distance_1 , int indice , float * y){
	return (6/(distance_1 + distance)) * ((((y[indice] - y[indice-1])/distance))-((y[indice - 1] - y[indice-2])/distance_1)); 
}

float ro(float lamda){
	return 1 - lamda;
}

float lamda(float distance ,float distance_1){
	return distance /(distance_1 + distance);
}

float distance_1(float * x , int indice ){
	return fabs(x[indice - 1] - x[indice]);
}

float distance(float * x , int indice){
	return fabs(x[indice] - x[indice - 1]);
}

void display(float * x , float * y , int dim){
	printf("=== Voici les données === \n");
	for(int i = 0 ; i < dim ; i++){
		printf("x[%d] = %f \t y[%d] = %f\n" , i , x[i] , i , y[i]);
	}
}

float * allocFloat1D(float * tab , int dim){
	tab = (float *)malloc(dim * sizeof(float));
	for(int i = 0 ; i < dim ; i++){
		tab[i] = 0;
	}
	return tab;
}

void getDataF(float ** x , float ** y , int dim){
	float * x_p = allocFloat1D(x_p  , dim);
	float * y_p = allocFloat1D(y_p , dim);

	FILE * fp = fopen("data.txt" , "r");
	if(fp){
		for(int i = 0 ; i < dim ; i++){
			fscanf(fp , "%f %f" , &x_p[i] , &y_p[i]);
		}
	}else{
		printf("Impossible d'ouvrir le fichier !!\n");
	}
	*x = x_p;
	*y = y_p;
}

float * Cholesky(){
	FILE * pf = fopen("matrice.txt" , "r");
	
	float ** A = NULL;
	float * B = NULL;
	float * X = NULL;
	int dim = 0;

	fscanf(pf , "%d" , &dim);
	A = allocFloat2D(A , dim);
	for(int i = 0 ; i < dim ; i++){
			for(int j = 0 ; j < dim ; j++){
					fscanf(pf , "%f" , &A[i][j]);
			}
	}

	B = allocFloat1D(B , dim);
	for(int i=0; i<dim; i++) {
			fscanf(pf, "%f", &B[i]);
	}

	X = allocFloat1D(X , dim);

	/// Factoriser A sachant qu'on travaille sur place
	float *y = malloc(dim * sizeof(float));
	for(int i = 0 ; i < dim ; i++){
		y[i] = 0;
	}
	/// Calcul de la matrice B
	// En utilisant A
	for(int i = 0 ; i < dim ; i++){
		for(int j = 0 ; j <= i ;j++){
			float sum = 0;
			if(j == i){
				for(int k = 0 ;k < i ; k++){
					sum += pow(A[j][k] , 2); 
				}
				A[j][j] = sqrt( A[j][j] - sum );
			}
			else {
				for(int k = 0 ; k < j ; k++){
					sum += A[i][k] * A[j][k];
				}
				A[i][j] = (A[i][j] - sum) / A[j][j];
			}
		}
	}	
	printf("\n");
	/// Resolution B * y = B
	// En utilisant A
	for(int i = 0 ; i < dim ; i++){
		float s = 0 ;
		for(int j = 0 ; j < i ; j++){
			s+=  y[j] * A[i][j];
		}
		y[i] =(B[i] - s) / A[i][i]; 
	}
	
	/// Resolution Bt x = y
	// En utilisant A
	for(int i = dim - 1 ; i >= 0 ; i--){
		float s = 0;
		for(int j = i + 1 ; j < dim ; j++){
			s+= A[j][i] * X[j];
		}
		X[i] = (y[i] - s ) / A[i][i];
		printf("%f\n" , X[i]);
	}
	return X;
}

