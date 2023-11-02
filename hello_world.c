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
	double** mat = make_matrix(m,n);
}
