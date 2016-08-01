
#ifndef ONEMATRIXPHYLOPROCESS_H
#define ONEMATRIXPHYLOPROCESS_H

#include "PhyloProcess.h"
#include "RandomSubMatrix.h"
#include "GTRSubMatrix.h"
#include "RandomTypes.h"
#include "ValArray.h"

// all branches and all sites are described by one single General Time Reversible (GTR) substitution process

class OneMatrixPhyloProcess : public PhyloProcess	{

	protected:

	public:

	OneMatrixPhyloProcess(LengthTree* intree, RandomSubMatrix* randmatrix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrix = randmatrix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	//
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(),0);
	}

	RandomSubMatrix* 	GetMatrix() {return matrix;}

	protected:
	RandomSubMatrix* matrix;
};

class SiteMatrixPhyloProcess : public PhyloProcess	{

	protected:

	public:

	SiteMatrixPhyloProcess(LengthTree* intree, RandomSubMatrix** randmatrix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrix = randmatrix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	//
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(site),0);
	}

	RandomSubMatrix* 	GetMatrix(int site) {return matrix[site];}

	protected:
	RandomSubMatrix** matrix;
};

class OneMatrixRASPhyloProcess : public OneMatrixPhyloProcess	{

	protected:

	public:

	OneMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate, RandomSubMatrix* randmatrix,  SequenceAlignment* indata) : OneMatrixPhyloProcess(intree,randmatrix,indata)	{
		rate = inrate;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new RandomBranchSitePath(this, tree->GetBranchLength(link->GetBranch()), GetRate(site), matrix, 0);
	}

	Var<PosReal>* GetRate(int site)	{
		if (! rate)	{
			return 0;
		}
		return rate->GetVal(site);
	}

	protected:

	VarArray<PosReal>* rate;
};


#endif
