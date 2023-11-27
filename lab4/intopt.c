#include<stdio.h>
#include<float.h>
#include<stdlib.h>
#include<math.h>
#include"simplex_t.h"

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


struct node_t initial_node(int m, int n, double** a, double* b, double* c) {
	struct node_t p;
	
	int i;
	p.a = calloc(m+1, sizeof(double*));
	for(i=0 ; i<m+1; i++)
	{
		p.a[i] = calloc(n+1,sizeof(double));
	}

	p.b = calloc(m+1, sizeof(double));
	p.c = calloc(n+1, sizeof(double));
	p.x = calloc(n+1, sizeof(double));
	p.min = calloc(n, sizeof(double));
	p.max = calloc(n, sizeof(double));
	p.m = m;
	p.n = n;	
	
	//copy a
	int y;
	for (i = 0; i< m+1;i++) {
		for(y=0; y<n+1; y++) {
			p.a[i][y] = a[i][y];
		}
	}
	for(i=0; i<m+1;i++){
		p.b[i] = b[i];
	}
	for(i=0; i<n+1; i++){
		p.c[i] = c[i];
	}

	for(i=0; i< n; i++){
		p.min[i] = DBL_MIN;
		p.max[i] = DBL_MAX;
	}
	return p;
}

struct node_t extend(struct node_t p,int m,int n ,double** a,double* b,double* c,int k,double ak,double bk) 
{
	struct node_t q;
	int i;
	int j;
	q.k = k;
	q.ak = ak;
	q.bk = bk;
	if ((ak > 0) && (p.max[k] < DBL_MAX) ) {
		q.m = p.m;
	}
	else if((ak < 0 )&& (p.min[k] > 0)) {
	       q.m = p.m;
	}
 	else {
		q.m = p.m +1;
	}
	q.n = p.n;
	q.h = -1;
	q.a = calloc(q.m + 1, sizeof(double*));
	for(i=0; i< m+1; i++) {
		q.a[i] = calloc(q.n+1, sizeof(double));
		
	}
	q.b = calloc(q.m+1, sizeof(double));
	q.c = calloc(q.n+1,sizeof(double));
	q.x = calloc(q.n+1, sizeof(double));
	q.min = calloc(n, sizeof(double));
	q.max = calloc(n, sizeof(double));

	for(i=0; i <n; i++) {
		q.min[i] = p.min[i];
		q.max[i] = p.max[i];
	}
	for(i=0; i<m;i++){
		for(j=0; j<q.n+1;j++){
			q.a[i][j] = a[i][j];
		}
		q.b[i] = b[i];
	}
	for(i=0; i < q.n+1; i++){
		q.c[i] = c[i];
	}
	if (ak > 0){
		if ((q.max[k] == DBL_MAX) || (bk < q.max[k])) {
			q.max[k] = bk;
		}
	}
	else if ((q.min[k] == DBL_MIN) || (-bk > q.min[k])) {
		q.min[k] = -bk;
	}
	for (i=m, j=0; j<n; j++) {
		if(q.min[j] > DBL_MIN) {
			q.a[i][j] = -1;
			q.b[i] = -q.min[j];
			i +=1;
		}
		if(q.max[j] < DBL_MAX) {
			q.a[i][j] = 1;
			q.b[i] = q.max[j];
			i += 1;
		}
	}
	return q;
}

int is_integer(double* xp) {
	double x = *xp;
	double r = lround(x);
	if (abs(r-x) < 1E-6) {
		*xp = r;
		return 1;
	}
	else {
		return 0;
	}
}

int integer(struct node_t p) {
	int i;
	for(i = 0; i<p.n; i++) {
		if(is_integer(&p.x[i]) == 0) {
			return 0;
		}
	}
	return 1;
}

void bound(struct node_t p,struct node_t* h, double* zp,double* x) {
	int i;
	if (p.z > *zp) {
		*zp = p.z;
		for(i=0; i<p.n;i++) {
			x[i] = p.x[i];
			//H SHOULD BE SET
			//remove and delete all nodes q in h with q.z p.z
		}
	}
}

int branch (struct node_t q,double z) {
	double min;
	double max;
	if (q.z< z) {
		return 0;
	}
	int h;
	for (h = 0; h < q.n; h++) {
		if (is_integer(&q.x[h]) == 0) {
				if(q.min[h] == DBL_MIN) {
					min=0;
				} else{
					min = q.min[h];
				}
				max = q.max[h];
				if ((q.x[h] < min) || (q.x[h] > max)) {
					continue;
				}
				q.h = h;
				q.xh = q.x[h];

				int i;
				for(i=0; i<q.m;i++) {
					free(q.a[i]);
				}
				free(q.a);
				free(q.b);
				free(q.c);
				free(q.x);
				return 1;
		}
		return 0;
	}
}

void free_node(struct node_t* h) {
	free(h->min);
	free(h->max);
	free(h->b);
	free(h->x);
	free(h->c);
	int i;
	for(i=0;i<h->m;i++) {
		free(h->a[i]);
	}
	free(h->a);
}

int size = 1;
void succ(struct node_t p, struct node_t* h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {
	struct node_t q = extend(p,m,n,a,b,c,k,ak,bk);
	if (&q == NULL) {
		return;
	}
	q.z = simplex(q.m,q.n,q.a,q.b,q.c,q.x,0);
	if (isfinite(q.z)) {
		if(integer(q) == 1) {
			bound(q,h,zp,x);
		} else if (branch(q, *zp) == 1) {
			//add q to h
			size++;
			struct node_t* h = realloc(h, size* sizeof(struct node_t));
			h[size-1] = q;		
			return;
		}
	}
	free_node(&q);
}

double intopt(int m, int n, double** a, double* b, double* c, double* x) {
	struct node_t p = initial_node(m,n,a,b,c);

	struct node_t* h = calloc(size, sizeof(struct node_t));
	h[0] = p;
	double z = DBL_MIN;
	p.z = simplex(p.m, p.n, p.a, p.b, p.c, p.x, 0);
	if(integer(p) == 1) {
		z = p.z;
		if (integer(p) == 1) {
			int i;
			//x size = n+1
			for( i=0; i < p.n +1;i++) {
				x[i] = p.x[i];
			}
			free_node(&p);
			return z;
		}
		branch(p,z);
		//0 should be empty
		while(h != 0) {
			//Take p from h
			struct node_t p = h[size -1];
			size--;	
			struct node_t* h = realloc(h,size* sizeof(struct node_t));	
			succ(p,h,m,n,a,b,c,p.h,1,p.xh,&z,x);

			succ(p,h,m,n,a,b,c,p.h,1,-p.xh,&z,x);
			free_node(&p);			
		}
		if (z == DBL_MIN) {
			return nan("0");
		} else {
			return z;
		}
	}
}

