#ifndef RESAMPLESUB_H
#define RESAMPLESUB_H

#include "ResampleSubMove.h"

class ResampleSubMove : public MCUpdate	{

	public:

	ResampleSubMove(PhyloProcess* inprocess, int innrep, double ins01, double ins10)	{
		process = inprocess;
		s01 = ins01;
		s10 = ins10;
		nrep = innrep;
	}

	// should I register paths ?

	double Move(double tuning)	{

		double nacc = 0;
		for (int rep=0; rep<nrep; rep++)	{
		double loghastings = process->MHProposeMove(this,s01,s10);

		process->Choose();

		process->SwapMatrices();

	}

	protected:

	int nrep;
	double s01;
	double s10;
	PhyloProcess* process;
};

#endif

