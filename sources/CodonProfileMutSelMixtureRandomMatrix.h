
#include "Random.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "Chronogram.h"
#include "BranchProcess.h"
// #include "PhyloProcess.h"
#include "GTRSubMatrix.h"
#include "MGCodonProfileMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"

class CodonProfileMutSelMixtureRandomMatrix : public MixtureRandomMatrix<Profile>	{

	public:

	CodonProfileMutSelMixtureRandomMatrix(Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inNeff)	{
		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = inmatrix;
		Neff = inNeff;
		CreateRandomVariable();
		CreateRandomSubMatrix();
		rvar->FullCorrupt();
		rvar->FullUpdate();
	}

	~CodonProfileMutSelMixtureRandomMatrix()	{
	}

	Rvar<Profile>* CreateRandomVariable()	{
		rvar = new Dirichlet(center, concentration);
		for (int i=0; i<rvar->GetDim(); i++)	{
			if ((*rvar)[i] < 1e-50)	{
				cerr << "i: " << i << "\t" << (*rvar)[i] << "\n";
				cerr << "very small entry in profile draw...exiting\n";
				exit(1);
			}
		}
		//rvar->Sample();
		rvar->SetName("rvar / profile");
		return rvar;
	}

	RandomSubMatrix* CreateRandomSubMatrix()	{
		matrix = new RandomMGCodonProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, Neff);
		matrix->SetName("codon mutsel matrix");
		return matrix;
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<PosReal>* Neff;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};


