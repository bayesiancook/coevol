#ifndef CONJUGATEONEMATRIXPHYLOPROCESS_H
#define CONJUGATEONEMATRIXPHYLOPROCESS_H


#include "OneMatrixPhyloProcess.h"
#include "ConjugatePath.h"


class ConjugateOneMatrixRASPhyloProcess : public OneMatrixRASPhyloProcess	{

	public:

	ConjugateOneMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate, RandomSubMatrix* inmatrix, SequenceAlignment* indata)  : OneMatrixRASPhyloProcess(intree, inrate, inmatrix, indata) {
		relrate = 0;
		stat = 0;
		GTRRandomSubMatrix* gtrmatrix = dynamic_cast<GTRRandomSubMatrix*>(matrix);
		if (gtrmatrix)	{
			relrate = gtrmatrix->GetRandomRelativeRates();
			stat = gtrmatrix->GetRandomStationaries();
			if (! stat)	{
				cerr << "error : gtr matrix does not have random stat\n";
				exit(1);
			}
			if (! relrate)	{
				cerr << "error : gtr matrix does not have random relrate\n";
				exit(1);
			}
		}
	}

	private :

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		/*
		cerr << link << '\t' << site << '\n';
		cerr << stat << '\t' << relrate << '\n';
		cerr << matrix << '\n';
		cerr << GetRate(site) << '\n';
		*/
		return new ConjugateRandomBranchSitePath(this, tree->GetBranchVal(link->GetBranch()), GetRate(site), matrix, stat, relrate);
	}

	Var<PosRealVector>* relrate;
	Var<Profile>* stat;
};

#endif

