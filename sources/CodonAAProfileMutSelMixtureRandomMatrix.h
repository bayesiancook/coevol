#ifndef CODONAAPROFILEMUTSELMIXTURERANDOMMATRIX_H
#define CODONAAPROFILEMUTSELMIXTURERANDOMMATRIX_H

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
#include "MGCodonAAProfileMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"

class CodonAAProfileMutSelMixtureRandomMatrix : public MixtureRandomMatrix<Profile>	{

	public:

	CodonAAProfileMutSelMixtureRandomMatrix(Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* incodonprofile, Var<PosReal>* inNeff)	{
		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = inmatrix;
		codonprofile = incodonprofile;
		Neff = inNeff;
		CreateRandomVariable();
		CreateRandomSubMatrix();
		nodemap.clear();
		rvar->FullCorrupt(nodemap);
		rvar->FullUpdate();
	}

	~CodonAAProfileMutSelMixtureRandomMatrix()	{
	}

	Rvar<Profile>* CreateRandomVariable()	{
		rvar = new Dirichlet(center, concentration);
		for (int i=0; i<Naa; i++)	{
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
		matrix = new RandomMGCodonAAProfileMutSelCodonSubMatrix(statespace, nucmatrix, rvar, codonprofile, Neff);
		matrix->SetName("codon aa mutsel matrix");
		return matrix;
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<Profile>* codonprofile;
	Var<PosReal>* Neff;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};


#endif
