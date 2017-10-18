#include "mex.h"
#include <limits.h>
#include <ppl.h>


void GetValidPatches(const int bigr, const int bigc, const int smallr, const int smallc, const int validSz,
	const mxLogical *sourceRegion, const unsigned int * validStarts, unsigned int *validPatches)
{
	//Concurrency::parallel_for(0, validSz, [&](unsigned int i)
	for (unsigned int i = 0; i < validSz; ++i)
	{
		unsigned int validStart = validStarts[i] - 1;
		unsigned int patchEnd = validStart + (bigr*(smallc - 1) + smallr - 1);
		if (validStart % bigr < bigr - smallr && // make sure patch doesnt start near the end of the img
			patchEnd < bigr * bigc) // make sure patch isnt passing the right side of the img
		{
			unsigned int rowIdx = 0;
			bool isValid = true;
			for (unsigned int j = validStart; j <= patchEnd; ++j)
			{
				if (rowIdx == smallr)
				{
					rowIdx = 0;
					j += bigr - smallr;
				}
				if (!sourceRegion[j])
				{
					isValid = false;
					break;
				}
				++rowIdx;
			}
			if (isValid)
			{
				validPatches[i * 2] = validStart + 1;
				validPatches[i * 2 + 1] = patchEnd + 1;
				//validPatches[i * 3 + 2] = (validStart + 1 + patchEnd + 1) >> 1; // center pixel of patch
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int mm, nn, m, n, vssz;
	unsigned int *validStarts, *validPatches;
	mxLogical *sourceRegion;

	mm = (int)mxGetScalar(prhs[0]);
	nn = (int)mxGetScalar(prhs[1]);
	m = (int)mxGetScalar(prhs[2]);
	n = (int)mxGetScalar(prhs[3]);
	vssz = (int)mxGetScalar(prhs[4]);
	sourceRegion = mxGetLogicals(prhs[5]);
	validStarts = (unsigned int*)mxGetData(prhs[6]);


	plhs[0] = mxCreateNumericMatrix(2, vssz, mxUINT32_CLASS, mxREAL);
	validPatches = (unsigned int*)mxGetData(plhs[0]);

	memset(validPatches, 0, 2 * vssz * sizeof(validPatches[0]));
	GetValidPatches(mm, nn, m, n, vssz, sourceRegion, validStarts, validPatches);
}

