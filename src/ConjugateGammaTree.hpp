#ifndef CONJUGATEGAMMATREE_H
#define CONJUGATEGAMMATREE_H

#include "BranchProcess.h"
#include "ConjugatePoisson.h"

class PoissonSemiConjugateTree : public BranchProcess<PosReal> {

	public:

	PoissonSemiConjugate* GetBranchVal(const Branch* branch)	{
		PoissonSemiConjugate* p = dynamic_cast<PoissonSemiConjugate*> (BranchProcess<PosReal>::GetBranchVal(branch));
		if (! p)	{
			if (BranchProcess<PosReal>::GetBranchVal(branch))	{
				cerr << "error in PoissonSemiConjugateTree: invalid dynamic cast\n";
				exit(1);
			}
		}
		return p;
	}

	PoissonSemiConjugateTree(Tree* intree) : BranchProcess<PosReal>(intree) {}

	void InactivateSufficientStatistic()	{
		InactivateSufficientStatistic(GetRoot());
	}

	void InactivateSufficientStatistic(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->InactivateSufficientStatistic();
			InactivateSufficientStatistic(link->Out());
		}
	}
};

class ConjugateGammaTree : public PoissonSemiConjugateTree {

	public:

	ConjugateGamma* GetBranchVal(const Branch* branch)	{
		ConjugateGamma* p = dynamic_cast<ConjugateGamma*> (PoissonSemiConjugateTree::GetBranchVal(branch));
		if (! p)	{
			if (BranchProcess<PosReal>::GetBranchVal(branch))	{
				cerr << "error in ConjugateGammaTree: invalid dynamic cast\n";
				exit(1);
			}
		}
		return p;
	}

	ConjugateGammaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : PoissonSemiConjugateTree(intree) {
		alpha = inalpha;
		beta = inbeta;
		RecursiveCreate(tree->GetRoot());
	}

	~ConjugateGammaTree()	{
		RecursiveDelete(tree->GetRoot());
	}

	void Integrate()	{
		Integrate(GetRoot());
	}

	void Integrate(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->Integrate();
			Integrate(link->Out());
		}
	}

	void Resample()	{
		InactivateSufficientStatistic(GetRoot());
	}

	protected:

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new ConjugateGamma(alpha,beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};

#endif

