#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv)
{
	printf("hello, world\n");
	int m;
	int n;
	scanf("%d %d", &m, &n);
	
	printf("m = %d, n = %d\n", m, n);
	
	fflush(stdout);
	
	double** make_matrix(int m, int n)
	{
		double** a;
		int i;
	
		a = calloc(m, sizeof(double*));
		for (i=0; i<m; i += 1)
			a[i] = calloc(n, sizeof(double));
		return a;
	}

	

	double* c_values = calloc(n, sizeof(double));
	double** a_values = make_matrix(m,n);	
	double* b_values = calloc(n,sizeof(double));
	
	int i;
	double* c;

	//read c
	for (i = 0; i < n; i+=1)
			
		scanf( "%lf", &c_values[i]);
		
	//read a
	int y;
	for (i = 0; i < m; i+=1)
		for(y= 0; y<n; y+=1)
			scanf("%lf",&a_values[i][y]);
	
	//read b
	for (i=0; i<n; i+=1)
		scanf("%lf",&b_values[i]);

	//Print result
	for (i=0; i <n; i+=1)	
		printf("%10.3lf ",c_values[i]);	
		fflush(stdout);
	
	printf("\n");
	
	fflush(stdout);
	
	for(i=0; i<m;i+=1)
		for (y=0; y<n; y+=1)
			printf("%10.3lf ",a_values[i][y]);							
			fflush(stdout);
		printf("\n");	
		fflush(stdout);
	
	for(i=0; i<n;i+=1)
		printf("%10.3lf ",b_values[i]);
		fflush(stdout);
	
	printf("\n");
	fflush(stdout);
	
}	



