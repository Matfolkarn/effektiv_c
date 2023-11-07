#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>


void print_values(int m, int n,double* c_values, double** a_values, double* b_values){
        int i;
        int y;

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
}


struct simplex_t{
	int	m;
	int	n;
	
	int* 	var;
	double**	a;
	double*	b;
	double*	x;
	double*	c;
	double	y;
};
//F02 /4


int init(struct simplex_t *s,int m, int n, double** a, double* b, double* c, double* x, double y, int* var){
	int i;
	int k;

	s->m = m;
	s->n = n;
	s->a = a;
	s->b = b;
	s->c = c;
	s->x = x;
	s->y = y;
	s->var = var;
	

	if (s->var == NULL){
		int* new_var;
		new_var = calloc(m+n+1, sizeof(int));
		s->var = new_var;
		for(i=0;i <m+n; i++) {
			s->var[i] = i;
		}
	}

	//Finding the lowest value in b
	for(k=0, i=1; i<m;i=i+1){
		if(s->b[i] < s->b[k]){
			k=i;
		}
	}

	return k;

}


//Finds out if there is a positive value of constant c, if so max z is not yet found.
int select_nonbasic(struct simplex_t *s) {
       int 	i;
	for(i=0; i<s->n; i++){
		
		//set epsilon to zero? or 10^-9
		if(s->c[i] > 0) {
			return i;
		}
	}
	return -1;
}	

int initial(struct simplex_t *s ,int m, int n, double** a, double* b, double* c, double* x, double y, int* var){
	int i;
	int k;
	int j;
	double w;
	
	k = init(s,m,n,a,b,c,x,y,var);

	if (s->b[k] > 0) {
		return 1;
	}	
	else {
		//Case if some b value is negative, then this will not work.
		return 0;
	}	
}

void pivot(struct simplex_t *s, int row, int col){
	double** a = s->a;
	double* b = s->b;
	double* c = s->c;
	int m = s->m;
	int n = s->n;
	int i;
	int j;
	int t;

	t = s->var[col];
	s->var[col] = s->var[n+row];
	s->var[n+row] = t;
	s->y = s->y + c[col]*b[row]/a[row][col];
	for(i=0; i<n;i++){
		if(i != col){
			c[i] = c[i] - c[col] * a[row][i] / a[row][col];
		}
	}
	c[col] = -c[col]/a[row][col];
	for(i=0; i<m; i++) {
		if (i != row) {
			b[i] = b[i] - a[i][col]*b[row]/a[row][col];
		}
	}
	for (i=0; i<m;i++) {
		if (i != row) {
			for (j=0; j<n; j++) {
				if (j != col){
					a[i][j] = a[i][j] - a[i][col] * a[row][j]/a[row][col];
				}
			}
		}
	}
	for (i = 0; i<m; i++){
		if (i != row) {
			a[i][col] = -a[i][col] /a[row][col];
		}
	}
	for (i=0; i<n; i++){
		if( i != col) {
			a[row][i] = a[row][i] / a[row][col];
		}
	}
	b[row] = b[row]/a[row][col];
	a[row][col] = 1/a[row][col];
}	

int xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h){
	struct simplex_t s;
	int i;
	int row;
	int col;

	int initiall;
	initiall = initial(&s, m,n,a,b,c,x,y,var);



	print_values(m,n,s.c,s.a, s.b);	
	
	if (initial == 0){
		//delete s.var free memory of it, why not all values? 
		free(s.var);
		return nan("NaN");
	}
	while((col = select_nonbasic(&s)) >=0) {
		row = -1;
		for(i=0; i<m; i++){
			//Can epsilon = 0?
			if((s.a[i][col] >0) && ((row < 0) || ((s.b[i]/s.a[i][col] )< (b[row] /a[row][col])))) {
				row = i;
			}				
		}
		if (row < 0) {
			free(s.var);
			return INT_MAX;
		}
		pivot(&s, row,col);
	
		print_values(m,n,s.c,s.a, s.b);	
	}
	if (h = 0) {
		for(i=0; i<n;i++){
			if (s.var[i] < n){
				x[s.var[i]] = 0;
			}
		}
		for (i=0; i<m;i++){
			if(s.var[n+1] < n){
				x[s.var[n+1]] = s.b[i];
			}
		}
		free(s.var);
	
	}	
	else {
		for(i=0; i<n; i++){
			x[i] = 0;
		}
		for(i=n; i<n+m;i++){
			x[i] = s.b[i-n];
		}
	}		
	return s.y;
}

int simplex(int m, int n, double** a, double* b, double* c,double* x,double y) {
	return xsimplex(m,n,a,b,c,x,y,NULL,0);
};

double** make_matrix(int m, int n)
        {
                double** a;
                int i;

                a = calloc(m, sizeof(double*));
                for (i=0; i<m; i += 1)
                        a[i] = calloc(n, sizeof(double));
                return a;
	}

void read_input(int *m, int *n, double* c_values, double** a_values, double* b_values) {
        
        scanf("%d %d", m, n);
	
	int i;
        
	
	c_values = calloc(*m, sizeof(double));
	//read c
        for (i = 0; i < *n; i+=1)
                scanf( "%lf", &c_values[i]);

	a_values = make_matrix(*m,*n);
        //read a
        int y;
        for (i = 0; i < *m; i+=1)
                for(y= 0; y<*n; y+=1)
                        scanf("%lf",&a_values[i][y]);
        
	b_values = calloc (*m,sizeof(double));
	//read b
        for (i=0; i<*n; i+=1)
                scanf("%lf",&b_values[i]);
}
int main() { 

//TEST
	int m1;
	int n1;
	double* c1;
	double** a1;
	double* b1;
	read_input(&m1,&n1,c1,a1,b1);
	printf("%d\n",m1);
//	print_values(m1,n1,c1,a1,b1);	
//TEST


	int res;
	
	//Allocatin input, should use function for this
	int m = 2;
	int n = 2;
	double** a = make_matrix(m,n);
	a[0][0] = -0.5;
	a[0][1] = 1;
	a[1][0] = 3;
	a[1][1] = 1;
	double* b = calloc(2, sizeof(double));
	b[0] = 4;
	b[1] = 18;
	double* c = calloc(2, sizeof(double));
	c[0] = 1;
	c[1] = 2;

	double* x = calloc(m, sizeof(double));//???
	double y = 0;//?????
	
	//Test print
	res = simplex(m,n,a,b,c,x,y);
	printf("\nOptimal value is: %d\n",res);
	//result print
}

