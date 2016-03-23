#include "mex.h"

/*
 * xtimesy.c - example found in API guide
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 */
/* $Revision: 1.10.6.4 $ */

void xtimesy(double *xr, double *xi, double *yr, double *yi,
        double *zr, double *zi, size_t m)
{
  mwSize i,j;
  
  for (i=0; i<m; i++) {
    for (j=0; j<i; j++) {
      *(zr+i+m*j) = *(xr+j) * *(yr+i) - *(xi+j) * *(yi+i);
      *(zi+i+m*j) = *(xr+j) * *(yi+i) + *(xi+j) * *(yr+i);
      *(zr+m*i+j) = *(zr+i+m*j);
      *(zi+m*i+j) = *(zi+i+m*j);
    }
    *(zr+i+m*i) = *(xr+i) * *(yr+i) - *(xi+i) * *(yi+i);
    *(zi+i+m*i) = *(xr+i) * *(yi+i) + *(xi+i) * *(yr+i);
  }
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *yr,*zr,*xr;
  double *yi,*zi,*xi;
  size_t mrowsx,mrowsy,ncols;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=2) 
    mexErrMsgIdAndTxt( "resp_mex:invalidNumInputs",
            "Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "resp_mex:invalidNumOutputs",
            "One output required.");
  
  if( !mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) )
      mexErrMsgIdAndTxt( "MATLAB:convec:inputsNotComplex",
              "Inputs must be complex.\n");
  
  
  /*  create a pointer to the input matrix x */
  xr = mxGetPr(prhs[0]);
  xi = mxGetPi(prhs[0]);
  
  /*  get the dimensions of the matrix input x */
  mrowsx = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if(ncols!=1) 
    mexErrMsgIdAndTxt( "resp_mex:invalidinput",
            "column vectors required.");
  
  /*  create a pointer to the input matrix y */
  yr = mxGetPr(prhs[1]);
  yi = mxGetPi(prhs[1]);
  
  /*  get the dimensions of the matrix input y */
  mrowsy = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  
  if(ncols!=1) 
    mexErrMsgIdAndTxt( "resp_mex:invalidinput",
            "column vectors required.");
  
  if(mrowsx!=mrowsy) 
    mexErrMsgIdAndTxt( "resp_mex:invalidinput",
            "x and y must be the same size.");
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix( (mwSize)mrowsx, (mwSize)mrowsy, mxCOMPLEX);
  
  /*  create a C pointer to a copy of the output matrix */
  zr = mxGetPr(plhs[0]);
  zi = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  xtimesy(xr,xi,yr,yi,zr,zi,mrowsx);
  
}
