
#ifndef DOLLOPHYLO_H
#define DOLLOPHYLO_H

#include "RandomBranchSitePath.h"
#include "PhyloProcess.h"
#include "DolloSubMatrix.h"


class DolloRandomBranchSitePath : public RandomBranchSitePath	{

	public:

	DolloRandomBranchSitePath(PhyloProcess* inprocess, Var<PosReal>* inlength, Var<PosReal>* inrate, RandomDolloSubMatrix* inmatrix, Var<Profile>* instationary) :
		RandomBranchSitePath(inprocess,inlength,inrate,inmatrix,instationary)	{

	}

	virtual void Resample()	{

		Reset(stateup);
		if (stateup != statedown)	{
			if (! stateup)	{
				cerr << "error in dollo: cannot undergo transition from 0 to 1\n";
				exit(1);
			}
			double totaltime = GetTime();
			double u = DrawWaitingTimeGivenAtLeastOne(stateup,totaltime);
			Append(0,u/totaltime);
			last->SetRelativeTime(1 - u/totaltime);
		}
		else	{
			last->SetRelativeTime(1);
		}
	}
};

class DolloOneMatrixPhyloProcess : public PhyloProcess	{

	protected:

	public:

	DolloOneMatrixPhyloProcess(LengthTree* intree, RandomDolloSubMatrix* randmatrix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrix = randmatrix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	//
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new DolloRandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(),0);
	}

	RandomDolloSubMatrix* 	GetMatrix() {return matrix;}

	protected:
	RandomDolloSubMatrix* matrix;
};

#endif

