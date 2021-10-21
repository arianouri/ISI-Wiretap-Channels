#include <mex.h>
#include <math.h>

#define cube3D prhs[0]
#define cube3D_out plhs[0]

void in_invCube(double *pD, double *pinD, mwSize s1, mwSize s2, mwSize s3);

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
      mwSize           nDimNum;
      const mwSize     *pDims;
      double           *pD;
      double           *pinD;
      mwSize           s1, s2, s3;
      //mxArray          *Data;
      
      nDimNum = mxGetNumberOfDimensions(cube3D);
      pDims = mxGetDimensions(cube3D);
      cube3D_out = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
      
      pinD = (double *) mxGetPr(cube3D);
      pD = (double *) mxGetPr(cube3D_out);
      s1 = pDims[0];
      s2 = pDims[1];
      s3 = pDims[2];
            
      in_invCube(pD,pinD,s1,s2,s3);
}

void in_invCube(double *pD, double *pinD, mwSize s1, mwSize s2, mwSize s3)
{
      int i, j, k;
      
      for(i=0; i<s1; i++)
        { for(j=0; j<s2; j++) 
         { for(k=0; k<s3; k++) 
          { pD[i + j * s1 + k * (s1 * s2)]= pinD[j + i * s2 + k * (s1 * s2)];}
         }
        }
}


