#include "mex.h"
#include <math.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(ABS)
#define	ABS(A)	((A) > (0) ? (A) : (-A))
#endif

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

double corr_coeff(const double *x, const double *y, const int N){
    double sum_x, sum_y, sum_sq_x, sum_sq_y, sum_coproduct, numerator;
    int i;
    
    sum_x = 0;
    sum_y = 0;
    sum_sq_x = 0;
    sum_sq_y = 0;
    sum_coproduct = 0;
    for (i=0; i<N; i++){
        sum_x += x[i];
        sum_y += y[i];
        sum_sq_x += x[i]*x[i];
        sum_sq_y += y[i]*y[i];
        sum_coproduct += x[i]*y[i];
    }
/*
    printf("%lf, %lf, %lf, %lf, %lf\n", sum_x, sum_y, sum_sq_x, sum_sq_y, sum_coproduct);
    printf("%lf   %lf   %lf\n", N*sum_coproduct - sum_x*sum_y, sqrt(N*sum_sq_x - sum_x*sum_x), sqrt(N*sum_sq_y - sum_y*sum_y));
 **/
    numerator = (N*sum_coproduct - sum_x*sum_y);
    if (ABS(numerator)<1e-10) {
        return 1;
    }
    return numerator / (sqrt(N*sum_sq_x - sum_x*sum_x) * sqrt(N*sum_sq_y - sum_y*sum_y));
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x, *y, *cc;
  int N, M, nY, mY, mrows, ncols, i, dims[2];
  dims[0] = 1;
  dims[1] = 1;
  
  /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Usage: cc = corr_coeff(x,y).");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /* Check each argument for proper form and dimensions. */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  N = MAX(mrows, ncols);
  M = MIN(mrows, ncols);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(M==1)) {
    mexErrMsgTxt("var must be a noncomplex double vector of length N.");
  }
  x = mxGetPr(prhs[0]);

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  nY = MAX(mrows, ncols);
  mY = MIN(mrows, ncols);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mY==1) || !(nY==N)) {
    mexErrMsgTxt("y must be a noncomplex double vector of length N (where N=length(x)).");
  }
  y = mxGetPr(prhs[1]);
  
  /* Create matrix for the second return argument (a scalar) and assign pointer. */
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  cc = mxGetPr(plhs[0]);
  
  cc[0] = corr_coeff(x, y, N);
}