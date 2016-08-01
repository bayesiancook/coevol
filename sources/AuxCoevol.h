#ifndef AUXCOEVOL_H
#define AUXCOEVOL_H

/*
template <class R, class D> class DSemiConjugateMove : public MCUpdate	{

	public:

	DSemiConjugateMove(R* inrandom, D* indsemi, double intuning, int inn) : random(inrandom), dsemi(indsemi), tuning(intuning), n(inn) {}

	double Move(double tuning_modulator = 1)	{
		dsemi->ActivateSufficientStatistic();
		double total = 0;
		for (int i=0; i<n; i++)	{
			total += random->Move(tuning* tuning_modulator);
		}
		total /= n;
		dsemi->InactivateSufficientStatistic();
		return total;
	}

	protected:

	R* random;
	D* dsemi;
	double tuning;
	int n;
};

class DSemiConjugateMappingMove : public MCUpdate	{

	public:

	DSemiConjugateMappingMove(PhyloProcess* inprocess, PathConjugateTree* inpathconjtree) : process(inprocess), pathconjtree(inpathconjtree) {}

	double Move(double tuning_modulator=1)	{
		pathconjtree->InactivateSufficientStatistic();
		process->Move(1);
		pathconjtree->ActivateSufficientStatistic();
		return 1;
	}

	protected:

	PhyloProcess* process;
	PathConjugateTree* pathconjtree;
};
*/

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};

#endif

