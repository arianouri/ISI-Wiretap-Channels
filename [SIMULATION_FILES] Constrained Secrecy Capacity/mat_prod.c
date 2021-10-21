#include <mex.h>
#include <math.h>

#define A_in prhs[0]
#define B_in prhs[1]
#define C_out plhs[0]


void in_mat_prod(double *A, double *B, double *C, int M_B, int M_A, int N_B, int N_A);

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	double *C,*A,*B;
    int M_A, N_A, M_B, N_B;

    A = mxGetPr(A_in);
	B = mxGetPr(B_in);
    
    M_A = mxGetM(A_in);
    N_A = mxGetN(A_in);
    
    M_B = mxGetM(B_in);
    N_B = mxGetN(B_in);
    
    if(N_A-M_B != 0)
        mexErrMsgTxt("Incorrect dimensions for matrix multiplication.");
    
	C_out = mxCreateDoubleMatrix(M_A, N_B, mxREAL);
	C = mxGetPr(C_out);

	in_mat_prod(A ,B, C, M_B, M_A, N_B, N_A);
	
}


void in_mat_prod(double *A, double *B, double *C, int M_B, int M_A, int N_B, int N_A)
{    
    int i, j, k;
    for(i=0;i<M_A;i++){    
        for(j=0;j<N_B;j++){    
            C[i+M_A*j]=0;    
            for(k=0;k<N_A;k++){    
                C[i+M_A*j]+=A[i+M_A*k]*B[k+M_B*j];    
            }    
        }    
    }   

}


