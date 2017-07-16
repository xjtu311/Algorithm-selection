/* fh_random_regtreefit_big_leaves_twofeaturetypes
Builds a regression tree from the input, picking features partly at random.
*/
 
#include "mex.h"
#include <stdlib.h> 
#include <math.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* Useful macro to swap 2 elements, using
 * the standard trick do { } while(0) to
 * make this one statement (otherwise it
 * could be wrong to use like 'else SWAP(A,B);').
 */
#define SWAP(A,B) \
do { int tmp = A; A = B; B = tmp; } while(0)

#define D_SWAP(A,B) \
do { double tmp = A; A = B; B = tmp; } while(0)

double Rand()
{
    return rand() / (double) (RAND_MAX);
}

/* Returns a sample from Normal(0,1) -- taken from Tom Minka's lightspeed
 */
double RandN(void)
{
  static double previous;
  static int usePrevious = 0;
  double x,y,radius;
  if(usePrevious) {
    usePrevious = 0;
    return previous;
  }
  /* Generate a random point inside the unit circle */
  do {
    x = 2*Rand()-1;
    y = 2*Rand()-1;
    radius = (x*x)+(y*y);
  } while((radius >= 1.0) || (radius == 0.0));
  /* Box-Muller formula */
  radius = sqrt(-2*log(radius)/radius);
  x *= radius;
  y *= radius;
  previous = y;
  usePrevious = 1;
  return x;
}

/* Returns a sample from Gamma(a, 1). -- taken from Tom Minka's lightspeed
 * For Gamma(a,b), scale the result by b.
 */
double GammaRand(double a)
{
  /* Algorithm:
   * G. Marsaglia and W.W. Tsang, A simple method for generating gamma
   * variables, ACM Transactions on Mathematical Software, Vol. 26, No. 3,
   * Pages 363-372, September, 2000.
   * http://portal.acm.org/citation.cfm?id=358414
   */
  double boost, d, c, v;
  if(a < 1) {
    /* boost using Marsaglia's (1961) method: gam(a) = gam(a+1)*U^(1/a) */
    boost = exp(log(Rand())/a);
    a++;
  } 
  else boost = 1;
  d = a-1.0/3; c = 1.0/sqrt(9*d);
  while(1) {
    double x,u;
    do {
      x = RandN();
      v = 1+c*x;
    } while(v <= 0);
    v = v*v*v;
    x = x*x;
    u = Rand();
    if((u < 1-.0331*x*x) || 
       (log(u) < 0.5*x + d*(1-v+log(v)))) break;
  }
  return( boost*d*v );
}


double mean(const double *array, const int N){
    int i;
    double sum_x;
    sum_x = 0;
    for(i=0; i<N; i++){
        sum_x += array[i];
    }
    return sum_x/N;
}

double geomean(const double *array, const int N){
    int i;
    double sum_log_x;
    sum_log_x = 0;
    for(i=0; i<N; i++){
        sum_log_x += log(array[i]);
    }
    return exp(sum_log_x/N);
}


double variance_n(const double *array, const int N){
    int i;
    double sum_x, sum_x_squared;
    if(N<=1){
        return 0;
    }
    sum_x = 0;
    sum_x_squared = 0;
    for(i=0; i<N; i++){
        sum_x += array[i];
        sum_x_squared += array[i]*array[i];
    }
/*
    printMatrixDouble(array, N, 1);
    printf("%lf, %lf, %d, %lf\n", sum_x, sum_x_squared, N, 1.0/(N-1) * (sum_x_squared - 1.0/N*sum_x*sum_x));
*/
    return fabs(1.0/(N) * (sum_x_squared - (1.0/N) * sum_x*sum_x));
}

double variance(const double *array, const int N){
    int i;
    double sum_x, sum_x_squared;
    if(N<=1){
        return 0;
    }
    sum_x = 0;
    sum_x_squared = 0;
    for(i=0; i<N; i++){
        sum_x += array[i];
        sum_x_squared += array[i]*array[i];
    }
/*
    printMatrixDouble(array, N, 1);
    printf("%lf, %lf, %d, %lf\n", sum_x, sum_x_squared, N, 1.0/(N-1) * (sum_x_squared - 1.0/N*sum_x*sum_x));
*/
    return fabs(1.0/(N-1) * (sum_x_squared - (1.0/N) * sum_x*sum_x));
}

/*=== Global variable for using built-in qsort since qsort_r can't be linked here.*/
double* ptr_to_double_array_for_qsort;

void quick(int *a, int min, int max) {
  if (max - min > 1) {

    int i = min;
    int j = max;
    /* indices are positive */
    int pivot = a[(i+j) >> 1];
    do {
      while(a[i] < pivot) i++;
      while(a[j] > pivot) j--;
      if (i > j) break;
      SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && a[j] == pivot) j--;
    if (min < j) quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && a[i] == pivot) i++;
    if (i < max) quick(a, i, max);

  } else if (a[min] > a[max])
    SWAP(a[min], a[max]);
}


void d_quick(double *a, int min, int max) {
  if (max - min > 1) {
    int i = min;
    int j = max;
    /* indices are positive */
    double pivot = a[(i+j) >> 1];
    do {
      while(a[i] < pivot) i++;
      while(a[j] > pivot) j--;
      if (i > j) break;
      D_SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && a[j] == pivot) j--;
    if (min < j) d_quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && a[i] == pivot) i++;
    if (i < max) d_quick(a, i, max);
  } else if (a[min] > a[max])
    D_SWAP(a[min], a[max]);
}


void dp_quick(int *a, int min, int max) {
  if (max - min > 1) {

    int i = min;
    int j = max;
    /* indices are positive */
    double pivot = ptr_to_double_array_for_qsort[a[(i+j) >> 1]];
    do {
      while(ptr_to_double_array_for_qsort[a[i]] < pivot) i++;
      while(ptr_to_double_array_for_qsort[a[j]] > pivot) j--;
      if (i > j) break;
      SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && ptr_to_double_array_for_qsort[a[j]] == pivot) j--;
    if (min < j) dp_quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && ptr_to_double_array_for_qsort[a[i]] == pivot) i++;
    if (i < max) dp_quick(a, i, max);

  } else if (ptr_to_double_array_for_qsort[a[min]] > ptr_to_double_array_for_qsort[a[max]])
    SWAP(a[min], a[max]);
}

/* Quick sort for sorting an int array a.
   Example: quick_sort(a, size);
// Equivalent: qsort(a, size, sizeof(int), compare_ints);
*/
void quick_sort(int *a, int size) {
  if (size > 1) {
    quick(a, 0, size-1);
  }
}

/* Quick sort for sorting a double array b.
   Example: d_quick_sort(b, size);
// Equivalent: qsort(b, size, sizeof(double), compare_doubles);
*/
void d_quick_sort(double *a, int size) {
  if (size > 1) {
    d_quick(a, 0, size-1);
  }
}

/* Quick sort for sorting the indices a of a double array b.
   Example: (int array a, double arrays b, c)
	for(i=0; i<size; i++) a[i]=i;
	ptr_to_double_array_for_qsort = b;
    dp_quick_sort(a, size);
//	Equivalent: qsort(a, size, sizeof(int), compare_idxs_into_double_array);

    If we also want the double array sorted, get it in array c:
    for(i=0; i<size; i++) c[i]=b[a[i]];*/
void dp_quick_sort(int *a, int size) {
  if (size > 1) {
    dp_quick(a, 0, size-1);
  }
}

/* For randperm in C */
/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(int *array, size_t n)
{
    if (n > 1) {
        size_t i;
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
    }
}


int compare_ints(const void * a, const void * b)
{
	const int *ia = (const int *) a;
	const int *ib = (const int *) b;
	
	return (*ia > *ib) - (*ia < *ib);
}


int compare_doubles (const void *a, const void *b){
	const double *da = (const double *) a;
	const double *db = (const double *) b;
 
	return (*da > *db) - (*da < *db);
}

/*=== Compare c[a] against c[b], where a and b are just indices. */
int compare_idxs_into_double_array (const void *a, const void *b){
	const int *ia = (const int *) a;
	const int *ib = (const int *) b;
	
	return (ptr_to_double_array_for_qsort[*ia] > ptr_to_double_array_for_qsort[*ib]) - (ptr_to_double_array_for_qsort[*ia] < ptr_to_double_array_for_qsort[*ib]);
}


void printMatrixDouble(double* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%lf ",X[i+nRows*j]);}
        printf(">\n");}
}

void printMatrixInt(int* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%d ",X[i+nRows*j]);}
        printf(">\n");}
}



void Rcritval_cat(const double *x, const double *Ycum, const int* rows, int nX, int nrows, double *critval_res, int *xleft, int *xright, int *numLeftPointer, int *numRightPointer){
	/* Declare Variables */
	int *t, *sorder, *diff_t, *n1, *maxlocs;
	int i, j, n, maxnumlocs, maxloc;
    
    double *B, *catmeans, *Ysplit1, *allx, *mu1, *mu2, *ssx;

	n = nrows + 1;

/*
	printf("x:\n");
    printMatrixDouble(x, nX, 1);  
	printf("\n\n");

	printf("Ycum:\n");
    printMatrixDouble(Ycum, nX, 1);  
	printf("\n\n");

	printf("rows:\n");
    printMatrixInt(rows , nrows, 1);  
	printf("\n\n");
*/

	/*=== Allocate Memory for Auxiliary Arrays. */
    t = mxCalloc(n,sizeof(int));
    sorder = mxCalloc(n,sizeof(int));
    diff_t = mxCalloc(n,sizeof(int));
    n1 = mxCalloc(n-1,sizeof(int));
    maxlocs = mxCalloc(n-1,sizeof(int));
    
    B = mxCalloc(n,sizeof(double));
	catmeans = mxCalloc(n,sizeof(double));
	Ysplit1 = mxCalloc(n-1,sizeof(double));
	allx = mxCalloc(n,sizeof(double));
	mu1 = mxCalloc(n-1,sizeof(double));
	mu2 = mxCalloc(n-1,sizeof(double));
	ssx = mxCalloc(n-1,sizeof(double));

    /*=== First get all possible split points. 
    //=== t are the changepoints + the last index   Matlab: t = [rows; size(Ycum,1)]; */
    for(i=0; i<n-1; i++){
        t[i] = rows[i]-1;
    }
    t[n-1] = nX-1;
/*
	printf("t:\n");
	printMatrixInt(t, n, 1);  
	printf("\n\n");
*/

    /*=== B contains the category sums.     Matlab: B = Ycum(t,:); B(2:end,:) = B(2:end,:) - B(1:end-1,:); */
    B[0] = Ycum[t[0]];
    for(i=1; i<n; i++){
        B[i] = Ycum[t[i]] - Ycum[t[i-1]];
    }
/*
	printMatrixDouble(B, n, 1);  
	printf("\n\n");
*/

    /*=== diff_t are the number of points in a category */
    diff_t[0] = t[0]+1;
    for(i=1; i<n; i++){
        diff_t[i] = t[i]-t[i-1]; 
    }
/*
	printMatrixInt(diff_t, n, 1);  
	printf("\n\n");
*/

	/*=== catmeans contains the means for the categories. */
	for(i=0; i<n; i++){
        catmeans[i] = B[i] / MAX(1, diff_t[i]);
    }
/*
	printMatrixDouble(catmeans, n, 1);  
	printf("\n\n");
*/
	
	/*=== Sort categories by mean response. Matlab: [smeans,sorder] = sort(catmeans); */
    for(i=0; i<n; i++){
        sorder[i] = i;
    }
	ptr_to_double_array_for_qsort = catmeans;
/*
	shuffle(sorder, n);    // shuffle to speed up qsort. 
	qsort(sorder, n, sizeof(int), compare_idxs_into_double_array);
*/	
    dp_quick_sort(sorder, n);
/*	qsort(sorder, n, sizeof(int), compare_ints);
//	qsort(ptr_to_double_array_for_qsort, n, sizeof(double), compare_doubles);
//	qsort(catmeans, n, sizeof(double), compare_doubles); */
	
/*
	printMatrixInt(sorder, n, 1);  
	printf("\n\n");
	
	printMatrixDouble(catmeans, n, 1);  
	printf("\n\n");
*/	
	
    /*=== Ysplit1 is Ycum(rows sorted by mean response)
    //=== n1(i) is the number of points going left when splitting using the ith subset*/
	for(i=0; i<n-1; i++){
        for(j=0; j<=i; j++){
            Ysplit1[i] = Ysplit1[i] + B[sorder[j]];
            n1[i] = n1[i] + diff_t[sorder[j]];
        }
    }
/*    
	printMatrixDouble(Ysplit1, n-1, 1);  
	printf("\n\n");
	printMatrixInt(n1, n-1, 1);  
	printf("\n\n");
*/	
	/*=== Take one x value from each unique set */
    for(i=0; i<n; i++){
        allx[i] = x[t[i]];
    }
/*    
	printf("\nallx:\n");
	printMatrixDouble(allx, n, 1);  
	printf("\n\n");
*/
	
	/*=== Get left/right means (same for cat/cont) */
    for(i=0; i<n-1; i++){
        mu1[i] = Ysplit1[i] / n1[i];
    }
    for(i=0; i<n-1; i++){
        mu2[i] = (Ycum[nX-1] - Ysplit1[i]) / (nX - n1[i]);
    }
/*
	printf("\nmu1:\n");
	printMatrixDouble(mu1, n-1, 1);  
	printf("\n\n");
	
	printf("\nmu2:\n");
	printMatrixDouble(mu2, n-1, 1);  
	printf("\n\n");
*/
	
	/*=== Get best split, val and location. */
    critval_res[0] = -10000000000000.0;
    for(i=0; i<n-1; i++){
        ssx[i] = n1[i]*mu1[i]*mu1[i] + (nX-n1[i])*mu2[i]*mu2[i];
        if(ssx[i] > critval_res[0]-1e-10){
        	if(ssx[i] > critval_res[0]+1e-10){
	            critval_res[0] = ssx[i];
	            maxnumlocs = 0;
	        }
	        maxlocs[maxnumlocs] = i;
        	maxnumlocs = maxnumlocs + 1;
		}
    }
/*
	printf("\nssx:\n");
	printMatrixDouble(ssx, n-1, 1);  
	printf("\n\n");

	printf("\nmaxlocs:\n");
	printMatrixInt(maxlocs, n-1, 1);  
	printf("\n\n");
*/
	maxloc = maxlocs[rand()%maxnumlocs];
/*	maxloc = maxlocs[0]; // TODO: remove this after debugging */
/*
	printf("%maxnumlocs = %d\n", maxnumlocs);
	printf("maxloc %d\n", maxloc);
*/

	numLeftPointer[0] = maxloc+1;
	numRightPointer[0] = n-maxloc-1;

	/* Now we can fill the result arrays xleft and xright as usual. */
	for(i=0; i<maxloc+1; i++){
	    xleft[i] = (int) allx[sorder[i]];
	}
	for(i=maxloc+1; i<n; i++){
	    xright[i-maxloc-1] = (int) allx[sorder[i]];
	}
/*
	printf("\nxleft  inside:\n");
	printMatrixInt(xleft, maxloc+1, 1);  
	printf("\n\n");

	printf("\nxright inside:\n");
	printMatrixInt(xright, n-maxloc-1, 1);  
	printf("\n\n");
*/

	/*=== Sort outputs. */
/*
	qsort(xleft, maxloc+1, sizeof(int), compare_ints);
	qsort(xright, n-maxloc-1, sizeof(int), compare_ints);
*/
    quick_sort(xleft, maxloc+1);
    quick_sort(xright, n-maxloc-1);

/*
	xleft[0] = 1;
	xleft[1] = 2;
	
	xright[0] = 3;
	xright[1] = 4;
	xright[2] = 6;
	xright[3] = 5;
*/
	/*=== Free Memory for Auxiliary Arrays. */
	mxFree(t);
	mxFree(sorder);
	mxFree(diff_t);
	mxFree(n1);
	mxFree(maxlocs);

	mxFree(B);
	mxFree(catmeans);
	mxFree(Ysplit1);
	mxFree(allx);
	mxFree(mu1);
	mxFree(mu2);
	mxFree(ssx);
}

void Rcritval_cont(const double *x, const double *Ycum, const int* rows, int nX, int nrows, double *critval_res, double *cutval_res){
	/* Declare Variables */
	int *maxlocs;
	int i, n, maxnumlocs, maxloc, cutloc;
    
    double u, *Ysplit1, *mu1, *mu2, *ssx;
	
	n = nrows + 1;

/*
    printMatrixDouble(x, nX, 1);  
	printf("\n\n");

    printMatrixDouble(Ycum, nX, 1);  
	printf("\n\n");

    printMatrixInt(rows , nrows, 1);  
	printf("\n\n");
*/

	/*=== Allocate Memory for Auxiliary Arrays. */
    maxlocs = mxCalloc(n-1,sizeof(int));
    
	Ysplit1 = mxCalloc(n-1,sizeof(double));
	mu1 = mxCalloc(n-1,sizeof(double));
	mu2 = mxCalloc(n-1,sizeof(double));
	ssx = mxCalloc(n-1,sizeof(double));

    /*=== Ysplit1 is Ycum(rows sorted by mean response) */
    for(i=0; i<n-1; i++){
        Ysplit1[i] = Ycum[rows[i]-1];
    }
	
	/*=== Get left/right means (same for cat/cont) */
    for(i=0; i<n-1; i++){
        mu1[i] = Ysplit1[i] / rows[i];
    }
    for(i=0; i<n-1; i++){
        mu2[i] = (Ycum[nX-1] - Ysplit1[i]) / (nX - rows[i]);
    }
/*
	printMatrixDouble(mu1, n-1, 1);  
	printf("\n\n");
	printMatrixDouble(mu2, n-1, 1);  
	printf("\n\n");
*/
	
	/*=== Get best split, val and location. */
    critval_res[0] = -10000000000000.0;
    for(i=0; i<n-1; i++){
        ssx[i] = rows[i]*mu1[i]*mu1[i] + (nX-rows[i])*mu2[i]*mu2[i];
        if(ssx[i] > critval_res[0]-1e-10){
        	if(ssx[i] > critval_res[0]+1e-10){
	            critval_res[0] = ssx[i];
	            maxnumlocs = 0;
	        }
	        maxlocs[maxnumlocs] = i;
        	maxnumlocs = maxnumlocs + 1;
		}
    }
/*
	printf("\nssx:\n");
	printMatrixDouble(ssx, n-1, 1);  
	printf("\n\n");

	printf("\nmaxlocs:\n");
	printMatrixInt(maxlocs, n-1, 1);  
	printf("\n\n");
*/

	maxloc = maxlocs[rand()%maxnumlocs];
/*	printf("%maxnumlocs = %d\n", maxnumlocs);
//	printf("maxloc %d\n", maxloc); */

/*	maxloc = maxlocs[0]; // TODO: remove this after debugging */
	

	/*=== Get cutval. */
	cutloc = rows[maxloc]-1;
    u = rand()/(RAND_MAX+0.0);

/*    printf("below: %lf, above: %lf\n", x[cutloc], x[cutloc+1]);   
    printf("u = %lf\n", u); */
/*	cutval_res[0] = (x[cutloc] + x[cutloc+1])/2;*/
    
    if(x[cutloc+1] - x[cutloc] < 1.9*1e-6){
        cutval_res[0] = (x[cutloc] + x[cutloc+1])/2;
    } else {
    	cutval_res[0] = ((1-u)*(x[cutloc]+1e-6) + u*(x[cutloc+1]-1e-6));
        if( cutval_res[0] < x[cutloc]+1e-8 || cutval_res[0] > x[cutloc+1]-1e-8 ){
            printf("below: %lf, above: %lf, u: %lf, chosen: %lf\n", x[cutloc], x[cutloc+1], u, cutval_res[0]);
            mexErrMsgTxt("random splitpoint has to lie in between the upper and lower limit");
        }
    }
    

	/*=== Free Memory for Auxiliary Arrays. */
	mxFree(Ysplit1);
	mxFree(mu1);
	mxFree(mu2);
	mxFree(ssx);
	mxFree(maxlocs);
}


void  buildTheTree(const double* X, const double* y, const double* orig_y, const int Splitmin, const int numFeaturesType1, const double p, const double percentageFeatures, const int* iscat, const mxArray* domains_cat, const int N, const int nvars, const double kappa, const double cutoff_penalty_factor, int* nodenumber, int* parent, double* yfitnode, mxArray* ysub, int* cutvar, double* cutpoint, int* leftchildren, int* rightchildren, double* resuberr, int* nodesize, mxArray* catsplit, double* leaf_g, double* leaf_m, int* leaf_n, int* numNodesPointer, int* numNcatsplitPointer){
	int i, j, k, offset, nRest, ncatsplit, tnode, nextunusednode, Nnode, bestvar, nRandom, numVarsToConsider, jvar, numrows, xcat, numBestLeft, numBestRight, nleft, nright, currnode, parent_node, catsplit_index, num_compatible, num_missing_to_left, num_missing_to_right;
	int *noderows, *leftside, *rightside, *assignednode, *randomPermutation, *idx, *rows;
	int *xleft, *xright, *numLeftPointer, *numRightPointer, *bestleft, *bestright, *xleftForResult, *xrightForResult, *compatible_values, *missing_values_for_left, *missing_values_for_right;
	double ybar, sst, mincost, bestcrit, bestcut, probForceSplitOnFeatureType1;
	double *xnoderow, *ynode, *x, *ycum, *critvalPointer, *cutvalPointer, *ysub_for_result;
	bool impure, ismember;
	int dims_left[2], dims_right[2], dims[2];
	mxArray *mx_xleft, *mx_xright, *mx_ysub, *mx_to_get_compatible_values;
    int K = 100;

	/*=== Allocate Memory for Auxiliary Arrays.*/
    noderows = mxCalloc(N,sizeof(int));
    leftside = mxCalloc(N,sizeof(int));
    rightside = mxCalloc(N,sizeof(int));
    assignednode = mxCalloc(N,sizeof(int));
    idx = mxCalloc(N,sizeof(int));
    rows = mxCalloc(N,sizeof(int));
    bestleft = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
    bestright = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */

    randomPermutation = mxCalloc(nvars,sizeof(int));
    
	xnoderow = mxCalloc(N,sizeof(double));
	ynode = mxCalloc(N,sizeof(double));
	x = mxCalloc(N,sizeof(double));
	ycum = mxCalloc(N,sizeof(double));
    
	/*=== Allocate Memory for Auxiliary Arrays for calling Rcritval_cat & Rcritval_cont.*/
	xleft = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
	xright = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
	numLeftPointer = mxCalloc(1,sizeof(int));
	numRightPointer = mxCalloc(1,sizeof(int));
	critvalPointer = mxCalloc(1,sizeof(double));
	cutvalPointer = mxCalloc(1,sizeof(double));

/*
	printf("X:\n");
	printMatrixDouble(X, N, nvars);
	printf("\n\n");
*/	
	/*=== Initialize variables. */
	ncatsplit = 0;
	nextunusednode = 2;
	nodenumber[0] = 1;
	resuberr[0] = 0;
	for(i=0; i<N; i++){
		assignednode[i] = 1;
	}
    
	/*=== Keep processing nodes until done. */
	tnode = 1;
	while(tnode < nextunusednode){ 
        /*=== Record information about this node.*/
		/* Matlab: noderows = find(assignednode==tnode); */
		Nnode = 0;
		for(i=0; i<N; i++){
	       if( assignednode[i] == tnode ){
	           noderows[Nnode++] = i;
	       }
		}
/*		printf("tnode: %d, Nnode is %d (must be >0)\n", tnode, Nnode);*/
        if( Nnode == 0 ){
            mexErrMsgTxt("Nnode is 0 !!!");
        }
            

		/* Matlab: ynode = y(noderows); */
		for(i=0; i<Nnode; i++){
			ynode[i] = y[noderows[i]];
/*            printf("ynode[%d]=%lf\n",i,ynode[i]);*/
		}
		
/*
		printf("ynode:\n");
		printMatrixDouble(ynode, Nnode, 1);
		printf("\n\n");
*/		
		/*=== Compute mean, variance and related statistics for this node. */
		/* Matlab: ybar = mean(ynode); */
		ybar = 0;
		for(i=0; i<Nnode; i++){
            		ybar += ynode[i]; 
		}
		ybar = ybar/Nnode;

		yfitnode[tnode-1] = ybar;

		/* Matlab: sst = norm(ynode-ybar)^2;   % total sum of squares at this node */
		sst = 0;
		for(i=0; i<Nnode; i++){
			sst = sst + (ynode[i]-ybar)*(ynode[i]-ybar);
		}

        if (Nnode > 1) {
            mincost = sqrt(sst / (Nnode-1)); /* stddev of ynode */
        } else {
            mincost = 0;
        }
		impure = (mincost > 1e-10 * resuberr[0]);

/*		printf("sst=%lf, mincost = %lf, impure=%d \n\n", sst, mincost, impure?1:0);*/
		
		/*=== Initialize variables before looping over possible split vars. */
		bestcrit          = -1e12;
		nodesize[tnode-1]   = Nnode;
		resuberr[tnode-1]   = mincost;
		cutvar[tnode-1]     = 0;
		cutpoint[tnode-1]   = 0;
		leftchildren[tnode-1] = 0;
		rightchildren[tnode-1] = 0;

        /*=== Consider splitting this node. */
		if ( (Nnode>=Splitmin) && impure ){     /* split only large impure nodes */
			/* Matlab: Xnode = X(noderows,:);  I don't want to deal with temporary matrices, so I'll work around by indexing. */
			
			bestvar = -1;
			bestcut = 0;

			/*=== First decision: force split on algorithm parameter? */
            probForceSplitOnFeatureType1 = pow( (log(N)/log(2) - log(Nnode)/log(2)) / (log(N)/log(2) - 1), p);
/*            probForceSplitOnFeatureType1 = 0.5; /*MAX(probForceSplitOnFeatureType1, 0.5);*/
            /*probForceSplitOnFeatureType1 = 0.8;*/
			/*      //probForceSplitOnFeatureType1= (p ^(1/(N-2.0)))^(Nnode-2); */

           /*
            printf("Nnode=%d, N=%d, probForceSplitOnFeatureType1=%lf, p_r=%lf\n",Nnode, N, probForceSplitOnFeatureType1, p);
            */
			if (probForceSplitOnFeatureType1*RAND_MAX > rand()){
				/*=== Force split on algorithm parameter. */
				nRandom = numFeaturesType1;
                nRest = nvars-numFeaturesType1;
 /*               printf("forceSplitOnFeatureType1 on\n");*/
			} else {
				/*=== Don't force any type of split variable. */
				nRandom = nvars;
                nRest=0;
/*                printf("forceSplitOnFeatureType1 off\n");*/
			}
			
			/* Matlab: randomPermutation = randperm(nRandom); */
			for(i=0; i<nvars; i++){
				randomPermutation[i] = i;
			}
			shuffle(randomPermutation, nRandom);    /* take out for debugging */ 

            if (nRest>0){
                shuffle(randomPermutation+nRandom, nRest); 
/*                printf("randperm forced:\n");*/
            } else {
/*                printf("randperm not forced:\n");*/
            }
/*            printMatrixInt(randomPermutation, nvars, 1); */


/*			printf("percentageFeatures=%lf\n",percentageFeatures);*/
			numVarsToConsider = MAX(1, (int) floor(percentageFeatures*nRandom));
/*            printf("nvars=%d, numVarsToConsider=%d\n",nvars,numVarsToConsider);*/

			/*=== Second decision: which one? */
			for(i=0; i<nvars; i++){ /* we allow numVarsToConsider options, but must have at least one real option */
				jvar=randomPermutation[i];
				xcat = iscat[jvar];
			
				/* Matlab: Xnoderow = Xnode(:,jvar) */
                offset = jvar*N; /* index into matrix: row + column*numRows */
				for(j=0; j<Nnode; j++){
					xnoderow[j] = X[noderows[j] + offset];  
				}
/*
				printf("jvar=%d, xnoderow:\n", jvar);
				printMatrixDouble(xnoderow, Nnode, 1);
				printf("\n\n"); 
*/				
		        /* Matlab: xnoderow=xnoderow(1:Nnode); %only sort first Nnode elements.
				// Matlab: [x,idx] = sort(xnoderow);          % get sorted jth x variable */
				ptr_to_double_array_for_qsort = xnoderow;
			    for(j=0; j<Nnode; j++){
			        idx[j] = j;
			    }
/*
				shuffle(idx, Nnode);    // shuffle to speed up qsort. 
				qsort(idx, Nnode, sizeof(int), compare_idxs_into_double_array);
*/

                dp_quick_sort(idx, Nnode);
                for(j=0; j<Nnode; j++){
			        x[j] = xnoderow[idx[j]];
			    }
/*
                printf("jvar=%d, x:\n", jvar);
				printMatrixDouble(x, Nnode, 1);
				printf("\n\n");
*/
                /*=== Determine if there's anything to split along this variable. */
				if (x[Nnode-1]-x[0] < 1e-10){
					continue;
				}
				
				/* Matlab: rows = find(x(1:end-1)+maxeps < x(2:end));
				// WATCH OUT: rows holds the indices as original in Matlab, not C style (but it itself is referenced standard C style starting at 0) */
				numrows = 0;
				for (j=0; j<Nnode-1; j++){
					if (x[j]+1e-10 < x[j+1]){ 
						rows[numrows++] = j+1; /* the +1 here is to make this compatible with calling Rcritval_cat & cont from Matlab directly. */
					}
				}
				if (numrows==0){
					continue;
				}

				/* Matlab: ycum = cumsum(ynode(idx,:) - ybar);  % centered response cum. sum */
				ycum[0] = ynode[idx[0]] - ybar;
				for (j=1; j<Nnode; j++){
					ycum[j] = ycum[j-1] + ynode[idx[j]] - ybar; /* % centered response cum. sum */
				}

				/*=== Do the core work: get the best split of the variable and its quality. */
				if (xcat>0){
					Rcritval_cat(x, ycum, rows, Nnode, numrows, critvalPointer, xleft, xright, numLeftPointer, numRightPointer);
/*
					printf("numleft=%d, xleft:\n", numLeftPointer[0]);
					printMatrixInt(xleft, numLeftPointer[0], 1);
					printf("\n\n");

					printf("numright=%d, xright:\n", numRightPointer[0]);
					printMatrixInt(xright, numRightPointer[0], 1);
					printf("\n\n");
*/
				} else {
					Rcritval_cont(x, ycum, rows, Nnode, numrows, critvalPointer, cutvalPointer);
				}

				/*=== Change best split if this one is best so far. */
				if (critvalPointer[0] > bestcrit + 1e-10){
					bestcrit = critvalPointer[0];
					bestvar = jvar;
					if (xcat>0){
						numBestLeft = numLeftPointer[0];
						numBestRight = numRightPointer[0];
						for(j=0; j<numBestLeft; j++){
							bestleft[j] = xleft[j];
						}
						for(j=0; j<numBestRight; j++){
							bestright[j] = xright[j];
						}
					} else {
						bestcut = cutvalPointer[0];
					}
				}
/*
                printf("i=%d, jvar = %d, xcat = %d, bestcrit = %lf, bestvar = %d, critval = %lf\n", i, jvar, xcat, bestcrit, bestvar, critvalPointer[0]);
*/
/*
				printf("jvar=%d, bestvar=%d, bestcrit=%lf:\n", jvar, bestvar, bestcrit);
				if (iscat[bestvar]){
					printf("numBestLeft=%d, bestleft:\n", numBestLeft);
					printMatrixInt(bestleft, numBestLeft, 1);
					printf("\n\n");

					printf("numBestRight=%d, bestright:\n", numBestRight);
					printMatrixInt(bestright, numBestRight, 1);
					printf("\n\n");
				} else {
					printf("bestcut=%lf\n", bestcut);
				}
*/				
				if (i >= numVarsToConsider-1 && bestcrit > -1e11){
/*                    printf("Have looked at %d variables, bestcrit=%lf, breaking\n",i+1,bestcrit);*/
					break;
				}
/*                printf("i=%d, bestcrit=%lf\n",i,bestcrit);*/
			}
/*
            printf("\n");
*/			
			/*=== Split this node using the best rule found. */
			if (bestvar == -1){
				/* Terminal node */
		        /*printf("Terminal node with %d data points and impure=%d\n", Nnode, impure?1:0);*/
        	} else {
				for (j=0; j<Nnode; j++){
					x[j] = X[noderows[j] + bestvar*N];
				}
				
				/*printf("\nnoderows:\n");
				printMatrixInt(noderows, Nnode, 1);  
				printf("\n\n");
				
				printf("\nx to be split:\n");
				printMatrixDouble(x, Nnode, 1);  
				printf("\n\n");*/
				
    
				if (iscat[bestvar]){
/*				printf("splitting on cat %d\n", bestvar);*/
					cutvar[tnode-1] = -(bestvar+1);          /* negative indicates cat. var. split */
					ncatsplit = ncatsplit + 1;  	   /* index into catsplit cell array */
					cutpoint[tnode-1] = ncatsplit;
					                    
                    /* 1: To get all compatible values, walk up the tree, looking
                     * for a split on the same parameter. If none is found
                     * take the initial domain of that parameter. */
                    currnode = tnode;
                    while (currnode > 1){
/*                        printf("currnode = %d\n", currnode);*/
                        parent_node = parent[currnode-1];
/*                        printf("parent_node = %d, cutvar[parent_node-1]=%d\n", parent_node, cutvar[parent_node-1]);*/
                        if (cutvar[parent_node-1] == -(bestvar+1)){
                            /* Take values from there, depending on whether which child we are */
                            catsplit_index = cutpoint[parent_node-1];
/*                            printf("catsplit_index = %d\n", catsplit_index);*/

                            if (leftchildren[parent_node-1] == currnode){
                                mx_to_get_compatible_values = mxGetCell(catsplit, catsplit_index-1);
                            } else {
                                if (! (rightchildren[parent_node-1] == currnode)){
                                    mexErrMsgTxt("currnode must either be left or right child of its parent.");
                                }
                                mx_to_get_compatible_values = mxGetCell(catsplit, catsplit_index-1+N);
                            }
                            break;
                        }
                        currnode = parent_node;
                    }
/*                    printf("done loop\n");*/
                    if (currnode == 1){
                        /* Get compatible values from initial domain. */
/*                        printf("bestvar=%d\n",bestvar);*/
                        mx_to_get_compatible_values = mxGetCell(domains_cat, bestvar);
                    }
                    /* Get compatible values from mx_to_get_compatible_values. */
                    num_compatible = mxGetNumberOfElements(mx_to_get_compatible_values);
                    compatible_values = (int*) mxGetData(mx_to_get_compatible_values);
                    
/*                    printf("num_compatible=%d\n",num_compatible);*/
                    /* 2: For each compatible but missing value choose a side u.a.r. */
                    missing_values_for_left = mxCalloc(num_compatible,sizeof(int));
                    missing_values_for_right = mxCalloc(num_compatible,sizeof(int));
                    num_missing_to_left = 0;
                    num_missing_to_right = 0;
                    for (i=0; i<num_compatible; i++){
                        for (j=0; j<numBestLeft; j++){
                            if (compatible_values[i] == bestleft[j]) break;
                        }
                        if (j == numBestLeft){
                            for (j=0; j<numBestRight; j++){
                                if (compatible_values[i] == bestright[j]) break;
                            }
                            if (j == numBestRight){
                                /* Missing but compatible value: choose side u.a.r. */
                                if (rand()%2 == 0){
                                    missing_values_for_left[num_missing_to_left++] = compatible_values[i];
                                } else {
                                    missing_values_for_right[num_missing_to_right++] = compatible_values[i];
                                }
                            }
                        }
                    }
/*                    printf("num_missing_to_left=%d\n",num_missing_to_left);
                    printf("num_missing_to_right=%d\n",num_missing_to_right);*/
                    
                    /* 3: Merge the determined and the randomly assigned missing values */
                    for (i=num_missing_to_left; i<num_missing_to_left+numBestLeft; i++){
                        missing_values_for_left[i] = bestleft[i-num_missing_to_left];
                    }
                    quick_sort(missing_values_for_left, num_missing_to_left+numBestLeft);

                    for (i=num_missing_to_right; i<num_missing_to_right+numBestRight; i++){
                        missing_values_for_right[i] = bestright[i-num_missing_to_right];
                    }
                    quick_sort(missing_values_for_right, num_missing_to_right+numBestRight);

                    /* 4: Put that information into the cell array. */
                                        
                    /*=== Set up the structures to fill the cell array output. */
					dims_left[0] = 1;
					dims_left[1] = num_missing_to_left+numBestLeft; /* numBestLeft; */
					
					dims_right[0] = 1;
					dims_right[1] = num_missing_to_right+numBestRight; /* numBestRight; */
					
					mx_xleft = mxCreateNumericArray(2, dims_left, mxINT32_CLASS, mxREAL);
					mx_xright = mxCreateNumericArray(2, dims_right, mxINT32_CLASS, mxREAL);
					
					mxSetCell(catsplit, ncatsplit-1, mx_xleft);
					mxSetCell(catsplit, ncatsplit-1+N, mx_xright);
			
					xleftForResult = (int*) mxGetData(mx_xleft);
					xrightForResult = (int*) mxGetData(mx_xright);

					/*=== Copy result from our temporary arrays to the ones associated with the output. */
					/*
                    for(i=0; i<numBestLeft; i++){
						xleftForResult[i] = bestleft[i];
					}
					for(i=0; i<numBestRight; i++){
						xrightForResult[i] = bestright[i];
					}
                    */
                    for(i=0; i<num_missing_to_left+numBestLeft; i++){
						xleftForResult[i] = missing_values_for_left[i];
					}
					for(i=0; i<num_missing_to_right+numBestRight; i++){
						xrightForResult[i] = missing_values_for_right[i];
					}
                    mxFree(missing_values_for_left);
                    mxFree(missing_values_for_right);
                    	
					/* Matlab: leftside = ismember(x,bestleftrightcell{1});
					// Matlab: rightside = ismember(x,bestleftrightcell{2}); */
					nleft = 0;
					nright = 0;

/*					printf("numBestLeft=%d\n", numBestLeft);
					printf("\nbestleft:\n");
					printMatrixInt(bestleft, numBestLeft, 1);  
					printf("\n\n");

					printf("numBestRight=%d\n", numBestRight);
					printf("\nbestright:\n");
					printMatrixInt(bestright, numBestRight, 1);  
					printf("\n\n");

					printf("\nx:\n");
					printMatrixDouble(x, Nnode, 1);  
					printf("\n\n");*/

					for (j=0; j<Nnode; j++){
                        
						ismember = false;
						for (k=0; k<numBestLeft; k++){
							if ( ((int) floor(x[j]+0.5)) == bestleft[k]) ismember = true;
						}
						
						if (ismember){
							leftside[nleft] = j;
							nleft = nleft+1;
						} else {
							rightside[nright] = j;
							nright = nright+1;
						}
					}
				} else {
/*					printf("splitting on cont %d at splitpoint %lf\n", bestvar, bestcut); */
					cutvar[tnode-1] = bestvar+1;
					cutpoint[tnode-1] = bestcut;
					
					/* Matlab: leftside = x<=bestcut; %logical   
					// Matlab: rightside = ~leftside; */
					nleft = 0;
					nright = 0;
					for (j=0; j<Nnode; j++){
						if (x[j] <= bestcut){
							leftside[nleft] = j;
							nleft = nleft+1;
						} else {
							rightside[nright] = j;
							nright = nright+1;
						}
					}
				}
/*				
				printf("\nleftside:\n");
				printMatrixInt(leftside, nleft, 1);  
				printf("\n\n");

				printf("\nrightside:\n");
				printMatrixInt(rightside, nright, 1);  
				printf("\n\n");
*/			

				leftchildren[tnode-1] = nextunusednode;
				rightchildren[tnode-1] = nextunusednode+1;
				for (j=0; j<nleft; j++){
					assignednode[noderows[leftside[j]]] = nextunusednode;
				}
				for (j=0; j<nright; j++){
					assignednode[noderows[rightside[j]]] = nextunusednode+1;
				}
	
/*				
				printf("\nassignednode:\n");
				printMatrixInt(assignednode, N, 1);  
				printf("\n\n");
*/				

				nodenumber[nextunusednode-1] = nextunusednode;
				nodenumber[nextunusednode-1+1] = nextunusednode+1;
				parent[nextunusednode-1] = tnode;
				parent[nextunusednode-1+1] = tnode;
				nextunusednode = nextunusednode+2;
			}
		} 

		if (leftchildren[tnode-1] == 0){
            /* Leaf => store results falling here (don't store them everywhere to avoid O(N^2) storage)*/
/*            printf("Leaf\n");*/            

            /*=== Set up the structures to fill the cell array output. */
            dims[0] = 1;
            dims[1] = Nnode;
					
            mx_ysub = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
            mxSetCell(ysub, tnode-1, mx_ysub);
            ysub_for_result = mxGetPr(mx_ysub);
            
            /* Save *runtimes*, not losses. */
            for(i=0; i<Nnode; i++){
                ynode[i] = orig_y[noderows[i]];
    		}

            for(i=0; i<Nnode; i++){
                ysub_for_result[i] = ynode[i];
    		}

            leaf_g[tnode-1] = geomean(ynode, Nnode);
            leaf_m[tnode-1] = mean(ynode, Nnode);
            leaf_n[tnode-1] = Nnode;
            
            /* Different meaning in leaves than elsewhere, b/c leaf preds are used by fillMatricesandComputeObjByLeafWalk*/
            yfitnode[tnode-1] = leaf_m[tnode-1];

            
/*            printf("leaf %d: leaf_g=%lf, leaf_m=%lf, leaf_n=%d\n", tnode-1, leaf_g[tnode-1], leaf_m[tnode-1], leaf_n[tnode-1]);*/
        }
        tnode = tnode+1;
	}

	numNodesPointer[0] = nextunusednode - 1;
	numNcatsplitPointer[0] = ncatsplit;
	
	/*=== Free Memory for Auxiliary Arrays. */
	mxFree(noderows);
	mxFree(leftside);
	mxFree(rightside);
	mxFree(assignednode);
	mxFree(idx);
	mxFree(rows);
	mxFree(bestleft);
	mxFree(bestright);

	mxFree(randomPermutation);

	mxFree(xnoderow);
	mxFree(ynode);
	mxFree(x);
	mxFree(ycum);
	
	mxFree(xleft);
	mxFree(xright);
	mxFree(numLeftPointer);
	mxFree(numRightPointer);
	mxFree(critvalPointer);
	mxFree(cutvalPointer);
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	double *X, *y, *tmp_double_ptr, *yfitnode, *cutpoint, *resuberr, *adjusted_y, *leaf_g, *leaf_m;
	double p, percentageFeatures, kappa, cutoff_penalty_factor;
	int *tmp_int_ptr, *iscat, *nodenumber, *parent, *cutvar, *leftchildren, *rightchildren, *nodesize, *numNodesPointer, *numNcatsplitPointer, *leaf_n, dim[1], dims[2];
	int numFeaturesType1, N, nvars, mrows, ncols, seed, Splitmin, i, buildLog10tree;

  /* Check for proper number of arguments. */
  if(nrhs!=12 || nlhs != 16) {
    mexErrMsgTxt("USAGE: [nodenumber, parent, yfitnode, ysub, cutvar, cutpoint, leftchildren, rightchildren, resuberr, nodesize, catsplit, leaf_g, leaf_m, leaf_n, numNodes, ncatsplit] = fh_random_regtreefit_big_leaves_twofeaturetypes_dist(X, y, Splitmin, numFeaturesType1, p, percentageFeatures, iscat, domains_cat, kappa, cutoff_penalty_factor, seed, buildLog10tree).");
  }
  
  /* Check each argument for proper form and dimensions. */
  N = mxGetM(prhs[0]);
  nvars = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("X must be a noncomplex double matrix.");
  }
  X = mxGetPr(prhs[0]);

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mrows==N) || !(ncols==1) ) {
    mexErrMsgTxt("y must be a noncomplex double column vector of the same length as size(X,1).");
  }
  y = mxGetPr(prhs[1]);

  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsInt32(prhs[2]) || mxIsComplex(prhs[1]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("Splitmin must be a noncomplex int scalar (cast it to int!).");
  }
  tmp_int_ptr = (int*) mxGetPr(prhs[2]);
  Splitmin = tmp_int_ptr[0];
  
  
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsInt32(prhs[3]) || mxIsComplex(prhs[3]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("numFeaturesType1 must be a noncomplex int scalar (cast it to int!).");
  }
  tmp_int_ptr = (int*) mxGetPr(prhs[3]);
  numFeaturesType1 = tmp_int_ptr[0];

  mrows = mxGetM(prhs[4]);
  ncols = mxGetN(prhs[4]);
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("p must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[4]);
  p = tmp_double_ptr[0];

  mrows = mxGetM(prhs[5]);
  ncols = mxGetN(prhs[5]);
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("percentageFeatures must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[5]);
  percentageFeatures = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[6]);
  ncols = mxGetN(prhs[6]);
  if( !mxIsInt32(prhs[6]) || mxIsComplex(prhs[6]) || !(mrows==nvars) || !(ncols==1) ) {
    mexErrMsgTxt("iscat must be a noncomplex int column vector of the same length as size(X,2).");
  }
  iscat = (int*) mxGetPr(prhs[6]);

  mrows = mxGetM(prhs[7]);
  ncols = mxGetN(prhs[7]);
  if( !mxIsCell(prhs[7]) || !(ncols==1) || (!mrows==nvars) ) {
    mexErrMsgTxt("domain_sizes_cat must be a nvars x 1 cell array (empty entries for cont. dimensions)");
  }

  mrows = mxGetM(prhs[8]);
  ncols = mxGetN(prhs[8]);
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("kappa must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[8]);
  kappa = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[9]);
  ncols = mxGetN(prhs[9]);
  if( !mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("cutoff_penalty_factor must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[9]);
  cutoff_penalty_factor = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[10]);
  ncols = mxGetN(prhs[10]);
  if( !mxIsInt32(prhs[10]) || mxIsComplex(prhs[10]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("seed must be a noncomplex int scalar.");
  }
  tmp_int_ptr = (int*) mxGetData(prhs[10]);
  seed = tmp_int_ptr[0];
  srand ( seed );

  mrows = mxGetM(prhs[11]);
  ncols = mxGetN(prhs[11]);
  if( !mxIsInt32(prhs[11]) || mxIsComplex(prhs[11]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("buildLog10tree must be a noncomplex int scalar.");
  }
  tmp_int_ptr = (int*) mxGetData(prhs[11]);
  buildLog10tree = tmp_int_ptr[0];
  
  
  /* Create vectors for return arguments and assign pointers. */
  /* These have to be of size 2*N since the number of nodes can be that big (well, 2N-1) */
  dim[0] = 2*N;

  plhs[0] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  nodenumber = (int*) mxGetData(plhs[0]); 

  plhs[1] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  parent = (int*) mxGetPr(plhs[1]);

  plhs[2] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  yfitnode = mxGetPr(plhs[2]);
  
  plhs[3] = mxCreateCellArray(1, dim); /*ysub*/

  plhs[4] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  cutvar = (int*) mxGetData(plhs[4]);

  plhs[5] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  cutpoint = mxGetPr(plhs[5]);

  plhs[6] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  leftchildren = (int*) mxGetData(plhs[6]);

  plhs[7] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  rightchildren = (int*) mxGetData(plhs[7]);
  
  plhs[8] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  resuberr = mxGetPr(plhs[8]);

  plhs[9] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  nodesize = (int*) mxGetData(plhs[9]);

  dims[0] = N;
  dims[1] = 2;
  plhs[10] = mxCreateCellArray(2, dims); /* catsplit */

  plhs[11] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  leaf_g = mxGetPr(plhs[11]);
  
  plhs[12] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  leaf_m = mxGetPr(plhs[12]);

  plhs[13] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  leaf_n = (int*) mxGetData(plhs[13]);
  
  dim[0] = 1;
  plhs[14] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  numNodesPointer = (int*) mxGetData(plhs[14]);

  plhs[15] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  numNcatsplitPointer = (int*) mxGetData(plhs[15]);

  /*=== Objective meanX, where X=cutoff_penalty_factor */
  /* We use this information directly to *build* the tree, so the tree
   * is already built using the metric we care about */
/*  printf("kappa = %lf\n", kappa);
  printMatrixDouble(y, 1, N);*/
/*  
  printf("y:\n");
  printMatrixDouble(y, 1, N);
*/
  /*=== Allocate Memory for Auxiliary Arrays. */
  adjusted_y = mxCalloc(N,sizeof(double));
  for(i=0; i<N; i++){
      if(y[i] > kappa-1e-4){
        adjusted_y[i] = kappa*cutoff_penalty_factor; /* kappa; /*y[i]; /*kappa*cutoff_penalty_factor;*/
      } else {
        adjusted_y[i] = y[i];
      }
      if (buildLog10tree){
          adjusted_y[i] = log10(adjusted_y[i]);
      }
  }
/*  
  printf("kappa=%lf, buildLog10tree=%d, adjusted_y:\n", kappa, buildLog10tree);
  printMatrixDouble(adjusted_y, 1, N);
  printMatrixDouble(adjusted_y, 1, N); 
  printMatrixDouble(y, 1, N); 
*/
  
  buildTheTree(X, adjusted_y, y, Splitmin, numFeaturesType1, p, percentageFeatures, iscat, prhs[7], N, nvars, kappa, cutoff_penalty_factor, nodenumber, parent, yfitnode, plhs[3], cutvar, cutpoint, leftchildren, rightchildren, resuberr, nodesize, plhs[10], leaf_g, leaf_m, leaf_n, numNodesPointer, numNcatsplitPointer);
  mxFree(adjusted_y);
}