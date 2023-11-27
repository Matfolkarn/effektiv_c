#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>


int xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

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

int glob;
void print_values(struct simplex_t *s){
        int i;
        int y;
	glob += 1;
        printf("max z = ");	
        for (i=0; i <s->n; i+=1)
        {
                printf("%10.3lf * x%d",s->c[i], i);
                
                if (i + 1 < s->n){
                        printf(" + ");
                }
        }
        printf("\n");

        for(i=0; i<s->m;i+=1)
        {
                printf("A%d: ",i);
                for (y=0; y<s->n; y+=1)
                {
                        printf("%+10.3lf * x%d ",s->a[i][y], y);

                        if (y+1<s->n) {
                                printf(" + ");
                        }
                }
                printf(" <=  %10.3lf\n", s->b[i]);
        }
        printf("\n");
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
       

	int i;
	
	for(i=0; i<s->n; i++){
		
		//set epsilon to zero? or 10^-9
		if(s->c[i] > 1E-6) {
			return i;
		}
	}
	return -1;
}

void prepare(struct simplex_t *s, int k) {
	printf("\nPREPARE\n");
	int m = s->m;
	int n = s->n;
	int i;
	for (i = m+n; i >n; i--){
		s->var[i] = s->var[i-1];
	
	}
	s->var[n] = m+n;
	n=n+1;
	for(i=0; i <m; i++) {
		s->a[i][n-1] = -1;
	}

	s->x = calloc(m+n,sizeof(double));
	s->c = calloc(n,sizeof(double));
	s->c[n-1] = -1;
	s->n = n;
	pivot(s,k, n-1);
}

int initial(struct simplex_t *s ,int m, int n, double** a, double* b, double* c, double* x, double y, int* var){
	int i;
	int k;
	int j;
	double w;
	
	k = init(s,m,n,a,b,c,x,y,var);

	if (s->b[k] >= 0) {
		return 1; //feasible
	}	
	else {
		//Case if some b value is negative, then this will not work.
		prepare(s, k);
		n = s->n;
		s->y = xsimplex(m,n,s->a, s->b, s->c, s->x, 0,s->var, 1);
		for(i=0; i<m+n; i++) {
			if(s->var[i] == m+n-1) {	
				if(s->x[i] > 1E-6) {
					free(s->x);
					free(s->c);
					return 0; //Går ej att lösa
				}	
				else {
					break;
				}
			}
		}
		if(i>=n){
			for(j = 0, k = 0; k <n; k++) {
				if( fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])) {
				       j = k;
				}
			}
			pivot(s,i-n,j);
			i = j;
		if(i<n-1) {
	 		k = s->var[i];
			s->var[i] = s->var[n-1];
			s->var[n-1] = k;
			for(k=0; k<m;k++) {
				w = s->a[k][n-1];
				s->a[k][n-1] = s->a[k][i];
				s->a[k][i] = w;		
			}
		} else {
			//// xn+m is nonbasic and last. forget it. ????
		}
		free(s->c);
		s->c = c;
		s->y = y;
		for(k=n-1; k< n+m-1; k++) {
			s->var[k] = s->var[k+1];
		}
		n = s->n = s->n-1;
		double* t;
	       	t = calloc(n,sizeof(double));
		for(k=0;k<n; k++) {
			for(j=0; j<n;j++) {
				if( k == s->var[j]){
					t[j] = t[j] + s->c[k];
					goto outer_loop_continue; // Does this really work?
				}
			
			}
			for(j=0; j<m; j++) {
				if(s->var[n+j] == k) {
					break;
				}
			}
			s->y = s->y + s->c[k]*s->b[j];
			for(i=0; i<n; i++) {
				t[i] = t[i] - s->c[k] * s->a[j][i];
			}
			outer_loop_continue:
			continue;

		}	
		for(i=0; i<n;i++) {
			s->c[i] = t[i];
		}
		free(t);
		free(s->x);
		return 1;
		
		}
	}	
}

int xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h){
struct simplex_t s;
	int i;
	int row;
	int col;

	int initiall;
	initiall = initial(&s, m,n,a,b,c,x,y,var);



	print_values(&s);	
	
	if (initiall == 0){
		//delete s.var free memory of it, why not all values? 
		free(s.var);
		return (int)nan("0");
	}
	while((col = select_nonbasic(&s)) >=0) {
		row = -1;
		for(i=0; i<m; i++){
			//Can epsilon = 0?
			if((s.a[i][col] >1E-6) && ((row < 0) || ((s.b[i]/s.a[i][col] )< (b[row] /a[row][col])))) {
				row = i;
			}				
		}
		if (row < 0) {
			free(s.var);
			return INT_MAX;
		}
		pivot(&s, row,col);
	
		print_values(&s);	
	}
	if (h == 0) {
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
	
	}	else {
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

	int	local_array[10];
	int i;
	for (i = 0; i < 11; i += 1)
		local_array[i] = i;

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

void read_input(int *m, int *n, double** c_values, double*** a_values, double** b_values) {
        
        scanf("%d %d", m, n);
	
	int i;
        
	
	*c_values = calloc(*n, sizeof(double));
	//read c
        for (i = 0; i < *n; i+=1)
                scanf( "%lf", &(*c_values)[i]);

	*a_values = make_matrix(*m,*n);
        //read a
        int y;
        for (i = 0; i < *m; i+=1)
                for(y= 0; y<*n; y+=1)
                        scanf("%lf",&(*a_values)[i][y]);
        
	*b_values = calloc (*m,sizeof(double));
	//read b
        for (i=0; i<*m; i+=1)
                scanf("%lf",&(*b_values)[i]);
}
int main() { 

	int m;
	int n;
	double* c;
	double** a;
	double* b;
	read_input(&m,&n,&c,&a,&b);



	int res;
	
	double* x = calloc(n+1, sizeof(double));//Varför n+1, i xsimplex loopas n+m igenom ???
	double y = 0;
	
	//Test print
	res = simplex(m,n,a,b,c,x,y);
	printf("\nOptimal value is: %d\n",res);
	//result print
	free(c);
	int i;
	for(i=0;i<m;i++){
		free(a[i]);
	}
	free(a);
	free(b);
	free(x);
}

