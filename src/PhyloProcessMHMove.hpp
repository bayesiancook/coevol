
#ifndef PHYLOMHMOVE_H
#define PHYLOMHMOVE_H

#include "Move.h"
#include "PhyloProcess.h"

class PhyloProcessMHMove : public MCUpdate	{

	public:

	PhyloProcessMHMove(PhyloProcess* inprocess, int innrep, double ins01, double ins10)	{
		process = inprocess;
		nrep = innrep;
		s01 = ins01;
		s10 = ins10;
	}

	double Move(double tuning_modulator){
		return process->MHMove(nrep,s01,s10);
	}

	protected:

	PhyloProcess* process;
	int nrep;
	double s01;
	double s10;
};


#endif


