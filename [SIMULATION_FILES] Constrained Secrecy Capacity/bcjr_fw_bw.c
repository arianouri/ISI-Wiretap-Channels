#include <mex.h>
#include <math.h>

#define ALPHA plhs[0]
#define BETA plhs[1]
#define SIGMA plhs[2]

void in_bcjr_fw_bw(double *alpha,
				   double *beta,
				   double *alpha_init,
				   double *beta_init,
				   double *post_state_MAR,
				   double *past_state_MAR,
				   mwSize state_Cardinal,
				   mwSize emState_length,
				   mwSize emState,
				   double *gamma,
				   int Window_Size,
				   double *sigma,
				   int *past_s_map);
                   
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    
//     double *alpha, *beta;
//     alpha=mxGetPr(plhs[0]);
//     beta=mxGetPr(plhs[1]);
    
    double *alpha_init, *beta_init;
    alpha_init=(double *)mxGetPr(prhs[0]);
    beta_init=(double *)mxGetPr(prhs[1]);

    double *post_state_MAR, *past_state_MAR;
    post_state_MAR=(double *)mxGetPr(prhs[2]);
	past_state_MAR=(double *)mxGetPr(prhs[3]);

    const mwSize *matDIMs_al_bt;
    mwSize state_Cardinal;
    matDIMs_al_bt=mxGetDimensions(prhs[0]);
    state_Cardinal=matDIMs_al_bt[0];
    
    const mwSize *matDIMs_post_state;
    mwSize emState;
    matDIMs_post_state=mxGetDimensions(prhs[2]);
    emState=matDIMs_post_state[1];
    
    double *gamma;
    gamma=(double *)mxGetPr(prhs[4]);
	const mwSize *matDIMS_gamma;
	mwSize emState_length;
	matDIMS_gamma= mxGetDimensions(prhs[4]);
	emState_length = matDIMS_gamma[1];

	int Window_Size;
	Window_Size = (int)mxGetScalar(prhs[5]);

	int *past_s_map;
	past_s_map = (int *)mxGetPr(prhs[6]);

	double *alpha, *beta;
	ALPHA = mxCreateDoubleMatrix(state_Cardinal, Window_Size+1, mxREAL);
	BETA = mxCreateDoubleMatrix(state_Cardinal, Window_Size+1, mxREAL);
	alpha = (double *)mxGetPr(ALPHA);
	beta = (double *)mxGetPr(BETA);
	
	double *sigma;
	SIGMA=mxCreateDoubleMatrix(emState*state_Cardinal, emState_length, mxREAL);
	sigma=(double *)mxGetPr(SIGMA);

    
    in_bcjr_fw_bw(alpha,beta,alpha_init,beta_init,post_state_MAR,past_state_MAR,state_Cardinal,emState_length,emState,gamma,Window_Size,sigma,past_s_map);
}

void in_bcjr_fw_bw(double *alpha,
                   double *beta,
                   double *alpha_init,
                   double *beta_init,
                   double *post_state_MAR,
				   double *past_state_MAR,
                   mwSize state_Cardinal,
				   mwSize emState_length,
                   mwSize emState,
                   double *gamma,
				   int Window_Size,
				   double *sigma,
				   int *past_s_map)
{
    int emBranch_tot=emState*state_Cardinal; 
    int scale_factor, scale_re_factor;
	int state_i, state_j, auxState_j, auxState_i, stIndex, i;
    double cumsum_alpha, cumsum_beta;
	int tb, ta;

	for (int k = Window_Size + 1; k < emState_length - Window_Size; k++) {

		for(i=0;i<state_Cardinal*(Window_Size+1);i++){
			alpha[i]=alpha_init[i];
			//beta[i+state_Cardinal*(emState_length-1)]=beta_init[i+state_Cardinal*(emState_length-1)];
			beta[i]=beta_init[i];
		}

		for(int wIndex=1; wIndex<=Window_Size; wIndex++){


			// {\alpha} and {\beta} calculation

				// Forward recursion
				ta = k - Window_Size - 1 + wIndex;
				scale_factor = 0;
				for (state_j = 0; state_j<state_Cardinal; state_j++) {
					for (auxState_i = 0; auxState_i<emState; auxState_i++) {

						state_i = past_state_MAR[state_j + state_Cardinal*auxState_i];
						scale_re_factor = past_s_map[scale_factor];
						alpha[state_j + wIndex*state_Cardinal] = alpha[state_j + wIndex*state_Cardinal] + alpha[state_i + (wIndex - 1)*state_Cardinal] * gamma[scale_re_factor + ta*(emBranch_tot)]; //FW

						scale_factor++;
					}
				}

				// Backward recursion
				tb = k + Window_Size - wIndex;
				scale_factor = 0;
				for(state_i=0; state_i<state_Cardinal; state_i++){
					for(auxState_j=0; auxState_j<emState; auxState_j++){

						state_j=post_state_MAR[state_i+state_Cardinal*auxState_j];

						beta[state_i + (Window_Size - wIndex)*state_Cardinal] = beta[state_i + (Window_Size - wIndex)*state_Cardinal] + beta[state_j + (Window_Size - wIndex + 1)*state_Cardinal] * gamma[scale_factor + (tb + 1)*(emBranch_tot)]; //BW

						scale_factor++;
					}
				}
				//Normalization
				cumsum_alpha = 0;
				cumsum_beta = 0;
				for (stIndex = 0; stIndex < state_Cardinal; stIndex++){
					cumsum_alpha = cumsum_alpha + alpha[stIndex + wIndex*state_Cardinal];
					cumsum_beta = cumsum_beta + beta[stIndex + (Window_Size - wIndex)*state_Cardinal];
				}
				for (stIndex = 0; stIndex < state_Cardinal; stIndex++){
					alpha[stIndex+ wIndex*state_Cardinal]=alpha[stIndex+ wIndex*state_Cardinal]/cumsum_alpha;
					beta[stIndex+ (Window_Size - wIndex)*state_Cardinal]=beta[stIndex+ (Window_Size - wIndex)*state_Cardinal]/cumsum_beta;
				}
		}

		//isgma calculation
		scale_factor = 0;
		for (state_i = 0; state_i<state_Cardinal; state_i++) {
			for (auxState_j = 0; auxState_j<emState; auxState_j++) {
				state_j = post_state_MAR[state_i + state_Cardinal*auxState_j];

				sigma[scale_factor + k*emBranch_tot] = alpha[state_i + (Window_Size)*state_Cardinal] * beta[state_j] * gamma[scale_factor + k*(emBranch_tot)];

				scale_factor++;
			}
		}
	
	}
}