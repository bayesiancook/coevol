
#include "BranchMatrixMixture.h"

class AAProfileMutSelMixtureRandomMatrixTree : public MixtureRandomMatrixTree<Profile>	{

	public:

	AAProfileMutSelMixtureRandomMatrixTree(Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, BranchVarTree<PosReal>* innefftree, Var<PosReal>* inrootneff)	{
		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		nefftree = innefftree;
		rootneff = inrootneff;
		CreateRandomVariable();
		RecursiveCreate(GetRoot());
		map<DAGnode*,int> tmpmap;
		rvar->FullCorrupt(tmpmap);
		rvar->FullUpdate();
	}

	~AAProfileMutSelMixtureRandomMatrixTree()	{
		RecursiveDelete(GetRoot());
		delete rvar;
	}

	Tree* GetTree() {return nefftree->GetTree();}

	//Var<PosReal>* GetRootNeff() {return rootneff;}

	Rvar<Profile>* CreateRandomVariable()   {
		rvar = new Dirichlet(center, concentration);

		for (int i=0; i<Naa; i++)	{
			if ((*rvar)[i] < 1e-50)	{
				cerr << "i: " << i << "\t" << (*rvar)[i] << "\n";
				cerr << "very small entry in profile draw...exiting\n";
				exit(1);
			}
		}

		//rvar->Sample();  already sampled
		//rvar->SetName("rvar / profile");
		//cerr << "displaying rvar...\n";
		//cerr << (*rvar) << "\n";
		return rvar;
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		RandomSubMatrix* matrix;
		if (link->isRoot())	{
			//matrix = new RandomMGAAProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, rootneff, rvar);
			matrix = new RandomMGAAProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, rootneff);
			//cerr << "root neff pointer : " << nefftree->GetBranchVal(link->GetBranch()) << '\n';
		}
		else	{
			//matrix = new RandomMGAAProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, nefftree->GetBranchVal(link->GetBranch()), rvar);
			matrix = new RandomMGAAProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, nefftree->GetBranchVal(link->GetBranch()));
		}
		matrix->SetName("matrix");
		return matrix;
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
	BranchVarTree<PosReal>* nefftree;
	Var<PosReal>* rootneff;
};
