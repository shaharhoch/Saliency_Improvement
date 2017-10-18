/**
 * A best exemplar finder.  Scans over the entire image (using a
 * sliding window) and finds the exemplar which minimizes the sum
 * squared error (SSE) over the to-be-filled pixels in the target
 * patch.
 *
 * @author Sooraj Bhat
 */

#include "mex.h"
#include <limits.h>
#include <ppl.h>
#include <vector>

void bestexemplarhelper(const int mm, const int nn, const int m, const int n, const int n_imgs,
	const double * img, const uint32_T *targetPatch, //const unsigned char *Ip, const mxLogical *toFill, 
    const mxLogical *sourceRegion, const uint32_T *vp, const double *gradInt,
	double *errors, double *intensities)
{
	//unsigned int err=0, patchErr = 0, bestErr = UINT_MAX;
	register unsigned int mmnn = mm*nn;
	
	//std::vector<unsigned int> errors(vpsize);
    //std::vector<double> intensities(vpsize);
	//memset(&errors[0], 0, vpsize * sizeof(errors[0]));
    //memset(&intensities[0], 0, vpsize * sizeof(intensities[0]));

	Concurrency::parallel_for<const unsigned int>(0, n_imgs, [&](unsigned int i_img) {
		unsigned int threeDOffset = mmnn * 3 * i_img;
		unsigned int oneDOffset = mmnn * i_img;

		Concurrency::parallel_for<const unsigned int>(0, mmnn, [&](unsigned int i) {
			//for(unsigned int i=0;i<(unsigned int)(vpsize); ++i) {	
			unsigned int vpidx = mmnn * 2 * i_img + i * 2;
			errors[oneDOffset + i] = std::numeric_limits<double>::max();
			intensities[oneDOffset + i] = std::numeric_limits<double>::max();
			if (vp[vpidx] != 0) // is patch valid?
			{
				unsigned int validStart = vp[vpidx] - 1;
				unsigned int afterValidEnd = vp[vpidx + 1];
				unsigned int targetIdx = 0, ri = 0, j = 0, sourceNdx = 0, targetNdx = 0;
				unsigned int err = 0, patchErr = 0;
				double patchGradInt = 0.0;
				for (j = validStart; j < afterValidEnd; ++j)
				{
					if (ri == m)
					{
						ri = 0;
						j += mm - m;
					}
					targetNdx = targetPatch[targetIdx] - 1;
					sourceNdx = j;
					if (sourceRegion[targetNdx])//toFill[fillIdx])
					{



						// 				patchErr += abs(img[ndx] - Ip[ndx2]);
						// 				patchErr += abs(img[ndx += mmnn] - Ip[ndx2 += mn]);
						// 				patchErr += abs(img[ndx += mmnn] - Ip[ndx2 += mn]);

						err = img[threeDOffset + sourceNdx] - img[targetNdx];
						patchErr += err*err;
						sourceNdx += mmnn;
						targetNdx += mmnn;
						err = img[threeDOffset + sourceNdx] - img[targetNdx];
						patchErr += err*err;
						sourceNdx += mmnn;
						targetNdx += mmnn;
						err = img[threeDOffset + sourceNdx] - img[targetNdx];
						patchErr += err*err;
					}
					else
					{
						patchGradInt += gradInt[oneDOffset + sourceNdx];
					}
					++ri;
					++targetIdx;
				}
				errors[oneDOffset + i] = patchErr;
				intensities[oneDOffset + i] = patchGradInt;
			}
		});
	});
//	unsigned int bestErr = UINT_MAX;
//    double bestIntensity = 0;
//	for (int i = 0; i < vpsize; ++i)
//	{
//		/*** Update ***/
//		unsigned int error = errors[i];
//		if (error < bestErr)
//		{
//			bestErr = error;
//			*best = vp[i*3+2];
//            bestIntensity = intensities[i];
//		}
//        else if(error == bestErr && intensities[i] < bestIntensity)
//        {
//            bestIntensity = intensities[i];
//            *best = vp[i*3+2];
//        }
//	}
}

/* best = bestexemplarhelper(mm,nn,m,n,img,Ip,toFill,sourceRegion); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int mm, nn, m, n, n_imgs;
	double *img;
    double *gradInt;
	double *errors;
    double *intensities;
	uint32_T *vp, *targetPatch;
	mxLogical *sourceRegion;

	/* Extract the inputs */
	mm = (int)mxGetScalar(prhs[0]);
	nn = (int)mxGetScalar(prhs[1]);
	m = (int)mxGetScalar(prhs[2]);
	n = (int)mxGetScalar(prhs[3]);
	n_imgs = (int)mxGetScalar(prhs[4]);
	img = (double *)mxGetData(prhs[5]);
	//Ip = (unsigned char *)mxGetData(prhs[6]);
	//toFill = mxGetLogicals(prhs[7]);
    targetPatch = (uint32_T*)mxGetData(prhs[6]);
	sourceRegion = mxGetLogicals(prhs[7]);
	vp = (uint32_T*)mxGetData(prhs[8]);
    gradInt = mxGetPr(prhs[9]);

	plhs[0] = mxCreateDoubleMatrix(mm*nn, n_imgs, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(mm*nn, n_imgs, mxREAL);
    errors = mxGetPr(plhs[0]);
    intensities = mxGetPr(plhs[1]);
    
	//memset(errors, 0, mm*nn * sizeof(errors[0]));
    //memset(intensities, 0, mm*nn * sizeof(intensities[0]));

	///* Do the actual work */
	bestexemplarhelper(mm, nn, m, n, n_imgs, img, targetPatch, sourceRegion, vp, gradInt, errors, intensities);

	///* Setup the output */
	//plhs[0] = mxCreateDoubleScalar(best);
}
