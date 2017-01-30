
#ifndef SELPHYLO_H
#define SELPHYLO_H

#include "AllocationTree.h"

class SelectionMatrixTree {

	public:

	// parameters:
	// we need a length tree (because branch lengths encode the mixture)
	// K: the number of classes
	// stats : the array of stationary probability vectors stat[0] ... stat[K-1], each of which is 20-dimensional

	SelectionMatrixTree(AllocationTree* intree, RandomSubMatrix*** inmatrix) {
		tree = intree;
		matrix = inmatrix;
	}
		

	// but return the matrix[k][i] where k is the relevant class and i is the site
	RandomSubMatrix* GetBranchSiteVal(const Branch* branch, const int site)	{
		//cerr << tree->GetBranchAllocation(branch) << '\n';

		return matrix[tree->GetBranchAllocation(branch)][site] ;

	}



	private:

	AllocationTree* tree;
	RandomSubMatrix*** matrix;
	// int site;

};

class SelectionPhyloProcess : public PhyloProcess	{


	public:

	SelectionPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate,  SelectionMatrixTree* inmatrix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{

		rate = inrate;
		matrix = inmatrix;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrix->GetBranchSiteVal(link->GetBranch(), site), 0);
	}

	protected:

	Var<PosReal>* GetRate(int site)	{
	    if (rate)	{
	      return rate->GetVal(site);
	    }
	    return 0;
	}

	VarArray<PosReal>* rate;
	SelectionMatrixTree* matrix;

};

#endif
