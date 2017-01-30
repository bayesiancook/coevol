
#include "GTRSubMatrix.h"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//		 GTRSubMatrix
// ---------------------------------------------------------------------------


GTRSubMatrix::GTRSubMatrix(int inNstate, const double* rr, const double* stat, bool innormalise) : SubMatrix(inNstate, innormalise)	{

	Nrr = Nstate * (Nstate-1) / 2;
	mRelativeRate = rr;
	if (stat)	{
		CopyStationary(stat);
	}
}

void GTRSubMatrix::CopyStationary(const double* instat)	{
	for (int k=0; k<Nstate; k++)	{
		mStationary[k] = instat[k];
	}
}

// ---------------------------------------------------------------------------
//		 ComputeArray
// ---------------------------------------------------------------------------

void	GTRSubMatrix::ComputeArray(int i)	{

	if (mRelativeRate)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				Q[i][j] = RelativeRate(i,j) * mStationary[j];
				total += Q[i][j];
			}
		}

		// should always ensure that the diagonal entry of the matrix Q[i][i] is such that
		// the sum over all entries of the row is equal to 0
		Q[i][i] = - total;
	}
	else	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] = mStationary[j];
		}
		Q[i][i] -= 1;
	}
}

