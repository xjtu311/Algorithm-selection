#include "mex.h"
#include <math.h>

/* Basic test of how fast nested loops are.
 * Result: 10 x 100 x 1000 x 100, each accessing a 100 x 1000 x 100 matrix takes 6 seconds 
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    int i, o_idx, s, p, t, N_s, nPi, nTheta;
    int *AllP_Prime;
/*    double one_over_kappa=0.1;*/
    int one_over_kappa=10;
    nTheta=100;
    nPi=1000;
    N_s=100; 
    
    /*=== Allocate Memory for Auxiliary Arrays. */
    printf("Allocating int array of size %d\n", nTheta*nPi*N_s);
    AllP_Prime = mxCalloc(nTheta*nPi*N_s, sizeof(int));

    for(i=0; i<nTheta*nPi*N_s; i++){
        AllP_Prime[i] = 0;
    }
    
    printf("Starting loop\n");
    for( o_idx=0; o_idx<10; o_idx++ ){
        for( s=0; s<N_s; s++ ){
            for( p=0; p<nPi; p++ ){
                for( t=0; t<nTheta; t++ ){
                    AllP_Prime[t*nPi*N_s + p*N_s + s] = (int) (o_idx+s+p+t) * one_over_kappa;
                }
            }
        }
    }
    printf("Done.\n");
}