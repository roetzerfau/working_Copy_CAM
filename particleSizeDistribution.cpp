#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include "mex.h"

class Pair  {
	public:
		int blockSize, number;
		Pair(int _blockSize, int _number) : blockSize(_blockSize), number(_number)  {}
		void increment()  { ++number; }
		bool operator== (const int bSize) const  { return blockSize == bSize; }
};

bool myLess(Pair a, Pair b)  { return a.blockSize < b.blockSize; }

int NX, numCells;

int upIndex(int i)  { return (i + NX) % numCells; }

int downIndex(int i)  { return ( (i - NX) + numCells ) % numCells; }

int leftIndex(int i)  { return (i - 1 + NX ) % NX + floor(i / NX) * NX; };

int rightIndex(int i)  { return (i + 1) % NX + floor(i / NX) * NX; };

int evaluateParticleSize(int i, int size, std::vector<bool>& alreadyDone)  {
	alreadyDone[i] = true;
	++size;
	if(!alreadyDone[upIndex(i)])	size = evaluateParticleSize(upIndex(i), size, alreadyDone);
	if(!alreadyDone[downIndex(i)])	size = evaluateParticleSize(downIndex(i), size, alreadyDone);
	if(!alreadyDone[leftIndex(i)])	size = evaluateParticleSize(leftIndex(i), size, alreadyDone);
	if(!alreadyDone[rightIndex(i)])	size = evaluateParticleSize(rightIndex(i), size, alreadyDone);
	return size;
}

void particleSizeDistribution(double* inMatrix, std::vector<Pair>& results)  {

	NX = floor(sqrt(numCells));

	int particleSize;
	std::vector<bool> alreadyDone(numCells, false);
	std::vector<Pair>::iterator it;

	for(int i = 0; i < numCells; ++i)
		if(inMatrix[i] != 1)  alreadyDone[i] = true;

	for(int i = 0; i < numCells; ++i)  {
		if (alreadyDone[i])  continue;
		particleSize = evaluateParticleSize(i, 0, alreadyDone);
		it = std::find(results.begin(), results.end(), particleSize);
		if (it == results.end()) results.push_back(Pair(particleSize,1));
		else  it->increment();
	}

	std::sort(results.begin(), results.end(), myLess);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if ( nrhs != 1 )
		mexErrMsgIdAndTxt("MATLAB:arrayFillGetPr:rhs","This function takes exactly one input argument.");

	clock_t begin = clock();

	double *inMatrix;
	double *outMatrix;
	std::vector<Pair> results;

	inMatrix = mxGetPr(prhs[0]);
	numCells = mxGetN(prhs[0]);

	particleSizeDistribution(inMatrix, results); 

	plhs[0] = mxCreateNumericMatrix(results.size(), 2, mxDOUBLE_CLASS, mxREAL);
	outMatrix = mxGetPr(plhs[0]);

	for (int index = 0; index < results.size(); index++ ) {
		outMatrix[index] = results[index].blockSize;
		outMatrix[results.size() + index] = results[index].number;
	}

	clock_t end = clock();

	mexPrintf("Ellapsed time for C++ algorithm is %f seconds.\n", double (end - begin) / CLOCKS_PER_SEC);

}
