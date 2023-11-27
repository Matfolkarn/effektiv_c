#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv)
{
	printf("hello, world\n");
	int m;
	int n;
	scanf("%d %d", &m, &n);
	
	printf("m = %d, n = %d\n", m, n);
	
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
	double* b_values = calloc(m,sizeof(double));
	
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
	printf("max z = ");
	for (i=0; i <n; i+=1)
	{	
		printf("%10.3lf * x%d",c_values[i], i);	
		fflush(stdout);
		if (i + 1 < n){
			printf(" + ");
		}
	}
	printf("\n");
	
	fflush(stdout);
	
	for(i=0; i<m;i+=1)
	{
		printf("A%d: ",i);
		for (y=0; y<n; y+=1)
		{
			printf("%+10.3lf * x%d ",a_values[i][y], y);	
			
			if (y+1<n) {
				printf(" + ");
			}
		}
		printf(" <=  %10.3lf\n", b_values[i]);	
	}
	printf("\n");
	fflush(stdout);
	free(c_values);
	for(i = 0;i <m;i++ ){
		free(a_values[i]);
	}
	free(a_values);
	free(b_values);
	
}	



