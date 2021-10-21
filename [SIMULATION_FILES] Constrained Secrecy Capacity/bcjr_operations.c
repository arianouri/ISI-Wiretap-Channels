#include <mex.h>
#include <math.h>
#define OUTPUT_LLR plhs[0]

void in_BCJR_operations(double *ouTpuT,
                        double isiLength,
                        int CodeLength,
                        int WindowSize,
                        double *D,
                        double *T,
                        mwSize s1,
                        mwSize s2,
                        mwSize s3,
                        double *alphaWin,
                        double *betaWin,
                        double *alphaZer,
                        double *betaZer,
                        double *b,
                        double *TrellisIndexPos,
                        double *TrellisIndexNeg,
                        int sizeTrellisIndexPosNeg);  //<<BCJR
void in_eye(double, double **);
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{   
    
    double *output;
    output = mxGetPr(plhs[0]);
    
    double isiLength;
    int CodeLength, WindowSize;
    isiLength=(double)mxGetScalar(prhs[0]);
    CodeLength=(int)mxGetScalar(prhs[1]);
    WindowSize=(int)mxGetScalar(prhs[2]);
    
    double *D, *T;
    D = (double *) mxGetPr(prhs[3]);
    T = (double *) mxGetPr(prhs[4]);
   
    const mwSize *cubeDims;
    mwSize s1, s2, s3;
    cubeDims = mxGetDimensions(prhs[4]);
    s1 = cubeDims[0];
    s2 = cubeDims[1];
    s3 = cubeDims[2];
    
    double *alphaWin, *betaWin, *alphaZer, *betaZer, *b;
    alphaWin = (double *) mxGetPr(prhs[5]);
    betaWin = (double *) mxGetPr(prhs[6]);
    alphaZer = (double *) mxGetPr(prhs[7]);
    betaZer = (double *) mxGetPr(prhs[8]);
    b = (double *) mxGetPr(prhs[9]);
    
	double *TrellisIndexPos;
	double *TrellisIndexNeg;
    int sizeTrellisIndexPosNeg;
    
    TrellisIndexPos = (double *)mxGetPr(prhs[10]);
    TrellisIndexNeg = (double *)mxGetPr(prhs[11]);
    sizeTrellisIndexPosNeg = (int)mxGetScalar(prhs[12]);
	
    double *ouTpuT;
	OUTPUT_LLR = mxCreateDoubleMatrix(1, CodeLength, mxREAL);
    ouTpuT = (double *)mxGetPr(OUTPUT_LLR);
	
    in_BCJR_operations(ouTpuT,isiLength,CodeLength,WindowSize,D,T,s1,s2,s3,alphaWin,betaWin,alphaZer,betaZer,b,TrellisIndexPos,TrellisIndexNeg,sizeTrellisIndexPosNeg);  //<<BCJR
}
//***********************************************************************************
//***********************************************************************************
void in_BCJR_operations(double *ouTpuT,
                        double isiLength,
                        int CodeLength,
                        int WindowSize,
                        double *D,
                        double *T,
                        mwSize s1,
                        mwSize s2,
                        mwSize s3,
                        double *alphaWin,
                        double *betaWin,
                        double *alphaZer,
                        double *betaZer,
                        double *b,
                        double *TrellisIndexPos,
                        double *TrellisIndexNeg,
                        int sizeTrellisIndexPosNeg) //<<BCJR
{    
	double sum_AuxAlphaZer_Vec, sum_AuxBetaZer_Vec;
	double b_Num, b_Den;
	double temp;
	int itemp;
    double StateCardinality;
    StateCardinality = pow(2,isiLength);

	double ** AuxAlphaZer = (double **)malloc(StateCardinality * sizeof(double *));
	for (int ii = 0; ii< StateCardinality; ii++) AuxAlphaZer[ii] = (double *)malloc(StateCardinality * sizeof(double));

	double ** AuxBetaZer = (double **)malloc(StateCardinality * sizeof(double *));
	for (int ii = 0; ii< StateCardinality; ii++) AuxBetaZer[ii] = (double *)malloc(StateCardinality * sizeof(double));

	double * auxCube_f;
	auxCube_f = (double *)malloc(s1 * s2 * s3 * sizeof(double));
	double * auxCube_b;
	auxCube_b = (double *)malloc(s1 * s2 * s3 * sizeof(double));

	double ** aux_AuxAlphaZer = (double **)malloc(StateCardinality * sizeof(double *));
	for (int ii = 0; ii< StateCardinality; ii++) aux_AuxAlphaZer[ii] = (double *)malloc(StateCardinality * sizeof(double));

	double ** aux_AuxBetaZer = (double **)malloc(StateCardinality * sizeof(double *));
	for (int ii = 0; ii< StateCardinality; ii++) aux_AuxBetaZer[ii] = (double *)malloc(StateCardinality * sizeof(double));

	double * AuxAlphaZer_Vec;
	AuxAlphaZer_Vec = (double *)malloc(s1 * sizeof(double));
	double * AuxBetaZer_Vec;
	AuxBetaZer_Vec = (double *)malloc(s1 * sizeof(double));



    for(int i=WindowSize; i<(CodeLength-WindowSize); i++){

        in_eye(StateCardinality, AuxAlphaZer);
		in_eye(StateCardinality, AuxBetaZer);

            for(int j=0; j<WindowSize; j++){
                
                            for(int a=0;a<s1;a++){    
                                for(int b=0;b<s2;b++){
                                    auxCube_b[a + b * s1 + (i-j-1) * (s1 * s2)]=0;
                                    for(int c=0;c<s1;c++){    
                                        auxCube_b[a + b * s1 + (i-j-1) * (s1 * s2)]
                                                +=D[a + c * s1 + (i-j-1) * (s1 * s2)] * T[b + c * s1 + (i-j-1) * (s1 * s2)];
                                    }    
                                }    
                            }

							
                            
                            for(int a=0;a<s1;a++){    
                                for(int b=0;b<s2;b++){
                                    (*(*(aux_AuxAlphaZer + a) + b))=0;
                                    for(int c=0;c<s1;c++){                                            
                                        (*(*(aux_AuxAlphaZer + a) + b))
                                                += (*(*(AuxAlphaZer + a) + c)) * auxCube_b[c + b * s1 + (i-j - 1) * (s1 * s2)];
                                    }    
                                }    
                            }

							for (int a = 0; a<s1; a++) {
								for (int b = 0; b<s2; b++) {
									temp = (*(*(aux_AuxAlphaZer + a) + b));
									(*(*(AuxAlphaZer + a) + b))=temp;
									}
								}
							//free(aux_AuxAlphaZer);
            

                            
                            for(int a=0;a<s1;a++){    
                                for(int b=0;b<s2;b++){
                                    auxCube_f[a + b * s1 + (i+j + 1) * (s1 * s2)]=0;
                                    for(int c=0;c<s1;c++){    
                                        auxCube_f[a + b * s1 + (i+j + 1) * (s1 * s2)]
                                                +=D[a + c * s1 + (i+j + 1) * (s1 * s2)] * T[b + c * s1 + (i+j + 1) * (s1 * s2)];
                                    }    
                                }    
                            }

							
                            
                            for(int a=0;a<s1;a++){    
                                for(int b=0;b<s2;b++){
                                    (*(*(aux_AuxBetaZer + a) + b))=0;
                                    for(int c=0;c<s1;c++){                                            
                                        (*(*(aux_AuxBetaZer + a) + b))
                                                += (*(*(AuxBetaZer + a) + c)) * auxCube_f[b + c * s1 + (i+j + 1) * (s1 * s2)];
                                    }    
                                }    
                            }
                      
							for (int a = 0; a<s1; a++) {
								for (int b = 0; b<s2; b++) {
									temp = (*(*(aux_AuxBetaZer + a) + b));
									(*(*(AuxBetaZer + a) + b)) = temp;
								}
							}
							//free(aux_AuxBetaZer);
            }
    //free(auxCube_b); 
    //free(auxCube_f);
    
    for(int a=0;a<s1;a++){    
            AuxAlphaZer_Vec[a]=0;
            for(int c=0;c<s1;c++){                                            
                AuxAlphaZer_Vec[a] 
                        += (*(*(AuxAlphaZer + a) + c)) * alphaWin[c];
            }    
    }
    for(int a=0;a<s1;a++){    
            AuxBetaZer_Vec[a]=0;
            for(int c=0;c<s1;c++){                                            
                AuxBetaZer_Vec[a] 
                        += (*(*(AuxBetaZer + a) + c)) * betaWin[c];
            }
    }
    
    sum_AuxAlphaZer_Vec=0;
    sum_AuxBetaZer_Vec=0;
    
    for(int a=0;a<s1;a++){
        sum_AuxAlphaZer_Vec += AuxAlphaZer_Vec[a];
        sum_AuxBetaZer_Vec += AuxBetaZer_Vec[a];
    }
    
    for(int a=0;a<s1;a++){
       alphaZer[a+i*s1] = AuxAlphaZer_Vec[a] / sum_AuxAlphaZer_Vec;
       betaZer[a+i*s1] = AuxBetaZer_Vec[a] / sum_AuxBetaZer_Vec;
    }
    
    for(int a=0;a<s1;a++){    
            b[a + i * s1]=0;
            for(int c=0;c<s1;c++){
                        b[a + i * s1] += T[c + a * s1 + (i) * (s1 * s2)] * alphaZer[c + i * s1];
            }    
    }
    
    for(int a=0;a<s1;a++){
       b[a + i * s1]=b[a + i * s1] * betaZer[a + i * s1];
    }
    b_Num = 0;
    b_Den = 0;
    
    for(int index=0;index<sizeTrellisIndexPosNeg;index++){
		itemp = (int) TrellisIndexPos[index] - 1;
        b_Num += b[itemp + i * s1];
		itemp = (int)TrellisIndexNeg[index] - 1;
        b_Den += b[itemp + i * s1];
    }
        
        ouTpuT[i]= b_Num / b_Den;
    }
}
//***********************************************************************************
//***********************************************************************************
/*void in_invCube(double *cubeIN, double *cubeOUT, mwSize s1, mwSize s2, mwSize s3)
{
      int i, j, k;
      for(i=0; i<s1; i++)
        { for(j=0; j<s2; j++) 
         { for(k=0; k<s3; k++) 
          { cubeOUT[i + j * s1 + k * (s1 * s2)]= cubeIN[j + i * s2 + k * (s1 * s2)];}
         }
        }
}*/
//***********************************************************************************
void in_eye(double N, double **eye_ar)
{
	//int **in_eye(int);                 HEADER
	//int** auxeye;						 OUTPUT VARIABLE
	//auxeye=in_eye(N);					 FUNCTION USAGE
	//eye[i][j] = *(*(auxeye + i) + j);	 2-D ARRAY INDEXING
	int M; M = N;
	//double ** eye_ar = (double **)malloc(N * sizeof(int *));
	//for (int i = 0; i< N; i++) eye_ar[i] = (int *)malloc(M * sizeof(int));

	int i, j;
	for (i = 0; i<N; i++) {
		for (j = 0; j<N; j++) {
			if (i == j)
				*(*(eye_ar + i) + j) = 1;
			//eye[i][j] = 1;
			else
				*(*(eye_ar + i) + j) = 0;
			//eye[i][j] = 0;
		}
	}
	//return eye_ar;
}
//***********************************************************************************
/*void in_mat_prod(double *A, double *B, double *C, int M_B, int M_A, int N_B, int N_A)
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
}*/
