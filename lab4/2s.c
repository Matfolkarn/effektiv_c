#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include<float.h>
#include <string.h>



double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

int nbr_simplex = 1;

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

void print_values(struct simplex_t *s){
        int i;
        int y;
        printf("max z = ");	
        for (i=0; i <s->n; i+=1)
        {
                printf("%10.3lf * x%d",s->c[i], i);
                
                if (i + 1 < s->n){
                        printf(" + ");
                }else {
			printf( " + %lf", s->y);
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
		if(s->c[i] > 1e-6) {
			return i;
		}
	}
	return -1;
}

void prepare(struct simplex_t *s, int k) {
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

	if (b[k] >= 0) {
	
		return 1; //feasible
	}	
	else {
		//Case if some b value is negative, then this will not work.
		prepare(s, k);
		n = s->n;
		s->y = xsimplex(m,n,s->a, s->b, s->c, s->x, 0,s->var, 1);
		for(i=0; i<m+n; i++) {
			if(s->var[i] == m+n-1) {	
				if(fabs(s->x[i]) > 1e-6) {
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
		}
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

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h){
struct simplex_t s;
	int i;
	int row;
	int col;

	int initiall;
	initiall = initial(&s, m,n,a,b,c,x,y,var);



	//print_values(&s);	
	
	if (initiall == 0){
		//delete s.var free memory of it, why not all values? 
		free(s.var);
		return NAN;
	}
	while((col = select_nonbasic(&s)) >=0) {
		row = -1;
		for(i=0; i<m; i++){
			//Can epsilon = 0?
			if((s.a[i][col] >1e-6) && ((row < 0) || ((s.b[i]/s.a[i][col] )< (b[row] /a[row][col])))) {
				row = i;
			}				
		}
		if (row < 0) {
			free(s.var);
			return INFINITY;
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
			if(s.var[n+i] < n){
				x[s.var[n+i]] = s.b[i];
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

double simplex(int m, int n, double** a, double* b, double* c,double* x,double y) {
	printf("\n------------\nSimplex: %d \n-------------------\n", nbr_simplex);
	nbr_simplex++;
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



struct node_t {
	int m;
	int n;
	int k;
	int h;
	double xh;
	double ak;
	double bk;
	double* min;
	double* max;
	double** a;
	double* b;
	double* x;
	double* c;
	double z;
};


///Set
struct linked_list {
	struct list_item *first_item;
};

struct list_item {
	struct list_item *next;
	struct node_t *elem;
};

struct node_t *linked_list_pop(struct linked_list *list) {
	struct list_item *pop = list->first_item;
	list->first_item = pop->next;
	struct node_t *p = pop->elem;
	free(pop);
	return p;
}

void free_linked_list(struct linked_list *list) {
	if (list->first_item != NULL) {
		struct list_item *l = list->first_item;
		struct list_item *next;
		while (l != NULL) {
			next = list->first_item->next;
			free(l);
			l = next;
		}
	}
	free(list);
}

int linked_list_size(struct linked_list *list) {
	if (list->first_item == NULL)
		return 0;
	struct list_item *l = list->first_item;
	int i = 0;

	while (l != NULL) {
		i++;
		l = l->next;
	}
	return i;
}

void linked_list_add(struct linked_list *list, struct node_t *p) {
	struct list_item *item = calloc(1, sizeof(struct list_item));
	item->elem = p;
	if (list->first_item != NULL)
		item->next = list->first_item;
	list->first_item = item;
}


//Set



void free_node(struct node_t* h) {
	free(h->min);
	free(h->max);
	free(h->b);
	free(h->x);
	free(h->c);
	int i;
	for(i=0;i<h->m+1;i++) {
		free(h->a[i]);
	}
	free(h->a);
	free(h);
}


int size = 1;
int presize;


struct node_t *initial_node(int m, int n, double** a, double* b, double* c) {
	struct node_t *p = calloc(1, sizeof(struct node_t));
	
	int i;
	p->a = calloc(m+1, sizeof(double*));
	for(i=0 ; i<m+1; i++)
	{
		p->a[i] = calloc(n+1,sizeof(double));
	}

	p->b = calloc(m+1, sizeof(double));
	p->c = calloc(n+1, sizeof(double));
	p->x = calloc(n+1, sizeof(double));
	p->min = calloc(n, sizeof(double));
	p->max = calloc(n, sizeof(double));
	p->m = m;
	p->n = n;	
	
	//copy a
	int y;
	for (i = 0; i< m;i++) {
		memcpy(p->a[i], a[i], n * sizeof(double));
	}
	memcpy(p->b, b, m * sizeof(double));
    memcpy(p->c, c, n * sizeof(double));

	for(i=0; i< n; i++){
		p->min[i] = -INFINITY;
		p->max[i] = INFINITY;
	}
	return p;
}

struct node_t *extend(struct node_t *p,int m,int n ,double** a,double* b,double* c,int k,double ak,double bk) 
{
	struct node_t* q = calloc(1, sizeof(struct node_t));
	int i;
	int j;
	q->k = k;
	q->ak = ak;
	q->bk = bk;
	if ((ak > 0) && (p->max[k] < INFINITY) ) {
		q->m = p->m;
	}
	else if((ak < 0 )&& (p->min[k] > 0)) {
	       q->m = p->m;
	}
 	else {
		q->m = p->m +1;
	}
	q->n = p->n;
	q->h = -1;
	q->a = (double**)calloc(q->m + 1, sizeof(double*));
	for(i=0; i< q->m+1; i++) {
		q->a[i] = (double*)calloc(q->n+1, sizeof(double));
		
	}
	q->b = calloc(q->m+1, sizeof(double));
	q->c = calloc(q->n+1,sizeof(double));
	q->x = calloc(q->n+1, sizeof(double));
	q->min = calloc(n, sizeof(double));
	q->max = calloc(n, sizeof(double));

	for(i=0; i <n; i++) {
		q->min[i] = p->min[i];
		q->max[i] = p->max[i];
	}
	for(i=0; i<m;i++){
		for(j=0; j<q->n;j++){
			q->a[i][j] = a[i][j];
		}
		q->b[i] = b[i];
	}
	for(i=0; i < n; i++){
		q->c[i] = c[i];
	}
	if (ak > 0){
		if ((q->max[k] == INFINITY) || (bk < q->max[k])) {
			q->max[k] = bk;
		}
	}
	else if ((q->min[k] == -INFINITY) || (-bk > q->min[k])) {
		q->min[k] = -bk;
	}
	for (i=m, j=0; j<n; j++) {
		if(q->min[j] > -INFINITY) {
			q->a[i][j] = -1;
			q->b[i] = -q->min[j];
			i +=1;
		}
		if(q->max[j] < INFINITY) {
			q->a[i][j] = 1.0;
			q->b[i] = q->max[j];
			i += 1;
		}
	}
	return q;
}

int is_integer(double* xp) {
	double x = *xp;
	
	double r = lround(x);
	if (fabs(r-x) < 1e-6) {
		*xp = r;
		return 1;
	}
	else {
		return 0;
	}
}

int integer(struct node_t *p) {
	int i;
	for(i = 0; i<p->n; i++) {
		if(!is_integer(&(p->x[i]))) {
			return 0;
		}
	}
	return 1;
}


void bound(struct node_t *p,struct linked_list *h, double* zp,double* x) {
	int i;
	if (p->z > *zp) {

		*zp = p->z;
		for(i=0; i<p->n;i++) {
			x[i] = p->x[i];
			//H SHOULD BE SET
		}
		//remove and delete all nodes q in h with q.z< p.z

		struct list_item *q, *prev, *next;
		if(h->first_item == NULL)
			return;


		while (h->first_item->elem->z < p->z)
		{
			q = h->first_item;
			h->first_item = q->next;
			//FREE?
			//free_node_t(q->elem);
			//free(q);
			if (h->first_item == NULL)
				return;
		}

		prev = h->first_item;
		q = prev->next;
		while (q != NULL)
		{
			next = q->next;
			if (q->elem->z < p->z)
			{
				prev->next = q->next;
				//FREE??
				// free_node_t(q->elem);
				// free(q);
			} else
			{
				prev = q;
			}
			q = next;
		}
		

	}
}

int branch (struct node_t * q,double z) {
	double min;
	double max;
	if (q->z< z) {
		return 0;
	}
	int h;
	for (h = 0; h < q->n; h++) {
		if (is_integer(&q->x[h]) == 0) {
				if(q->min[h] == -INFINITY) {
					min=0;
				} else{
					min = q->min[h];
				}
				max = q->max[h];
				if (floor(q->x[h] < min) || ceil(q->x[h] > max)) {
					continue;
				}
				q->h = h;
				q->xh = q->x[h];

				int i;
				for(i=0; i<q->m;i++) {
					free(q->a[i]);
				}
				free(q->a);
				free(q->b);
				free(q->c);
				free(q->x);
				return 1;
		}
	}
	return 0;
	
}



void succ(struct node_t *p, struct linked_list* h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {
	struct node_t *q = extend(p,m,n,a,b,c,k,ak,bk);
	if (q == NULL) {
		return;
	}
	q->z = simplex(q->m,q->n,q->a,q->b,q->c,q->x,0);
	if (isfinite(q->z)) {
		if(integer(q) == 1) {
			bound(q,h,zp,x);
		} else if (branch(q, *zp) == 1) {
			linked_list_add(h,q);	
			return;
		}
	}
	free_node(q);
}



//:int simplex(int m, int n, double** a, double* b, double* c,double* x,double y) {
double intopt(int m, int n, double** a, double* b, double* c, double* x) {
	struct node_t *p = initial_node(m,n,a,b,c);
	struct linked_list* h = calloc(1, sizeof(struct linked_list));
	h->first_item = calloc(1,sizeof(struct list_item));
	h->first_item->elem = p;
	double z = -INFINITY;
	p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
	
	if(integer(p) || !isfinite(p->z)) {
		z = p->z;
		if (integer(p)) {
			int i;
			//x size = n+1
			for( i=0; i < p->n +1;i++) {
				x[i] = p->x[i];
				printf("\nx_%d: = %lf", i, x[i]);
			}
			free_node(p);
			//free(h);
			return z;
		}
	}

	branch(p,z);
	//0 should be empty
	while(h->first_item != NULL) {
		//Take p from h
		p = linked_list_pop(h); 

		succ(p,h,m,n,a,b,c,p->h,1,floor(p->xh),&z,x);

		succ(p,h,m,n,a,b,c,p->h,-1,-ceil(p->xh),&z,x);
		free(p->min);
		free(p->max);
		free(p);
	}
	if (z == -INFINITY) {
		//free(h);
		return NAN;
	} else {
		int i;
		//for( i=0; i < p->n +1;i++) {
		//		x[i] = p->x[i];
		//		printf("\nx_%d: = %lf", i, x[i]);
		//	}
		//free(h);
		
		return z;
	}
	
}


double main() { 

	int m;
	int n;
	double* c;
	double** a;
	double* b;
	read_input(&m,&n,&c,&a,&b);



	double res;
	
	double* x = calloc(n+1, sizeof(double));//Varför n+1, i xsimplex loopas n+m igenom ???
	double y = 0;
	
	//Test print
	//res = simplex(m,n,a,b,c,x,y);
	res = intopt(m,n,a,b,c,x);
	printf("\nOptimal value is: %lf\n",res);
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
