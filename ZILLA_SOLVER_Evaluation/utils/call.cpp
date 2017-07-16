#include <mex.h>
#include <stdlib.h>
#include <errno.h>


void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  /* A pointer to the commandline we want to run. */
  char *sText = 0; 
  
  /* Check that we are properly prototyped. */
  if ( nlhs > 0 ) mexErrMsgTxt ( "No output arguments" );
  if ( nrhs > 1 ) mexErrMsgTxt ( "Not more than one input argument." );
  
  /* Set the commandline to what we want. */
  if ( nrhs == 1 ) {
    
    /*Get the commandline.*/
    sText = mxArrayToString ( prhs[0] );    
    
    /* Do the actual program call, and check for errors. */
    int retVal = system ( sText );
    if ( retVal == 0 ) {
      /* No error. Continue as normal. */
    } else {
      mexPrintf ( "System call returned: %d\n", retVal );
      mexErrMsgTxt( "Error calling program." );
    }   
  } else if ( nrhs == 0 ) { /* No input means we do a test run. */
    mexPrintf( "No input command given. Checking for command interpreter:\n" );
    if ( system(0) != 0 ) {
      mexPrintf ( "Command processor found\n" );
    } else {
      mexErrMsgTxt( "No command processor found" );
    }
  }
  return;
}
