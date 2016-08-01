#include "ConjugateMultiVariateTreeProcess.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "Chronogram.h"
#include "BranchProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"

#include "MultiNormalModel.h"

class ConjugateInverseWishartResampleMove : public MCUpdate	{

	public:
	ConjugateInverseWishartResampleMove(ConjugateInverseWishart* insigma)	{
		sigma = insigma;
	}

	double Move(double tuning_modulator = 1)	{
		sigma->Integrate();
		sigma->Resample();
		return 1;
	}

	protected:
	ConjugateInverseWishart* sigma;

};

class ConjugateBranchOmegaMultivariateModel : public BranchOmegaMultivariateModel {

	public:

	ConjugateInverseWishart* GetConjugateInverseWishart() {return dynamic_cast<ConjugateInverseWishart*>(sigma);}
	ConjugateMultivariateNormal* GetConjugateMultivariateNormal(int i) {return dynamic_cast<ConjugateMultivariateNormal*>(data[i]);}

	ConjugateBranchOmegaMultivariateModel(string contdatafile, bool sample=true)	{

		clampdiag = false;

		contdata = new FileContinuousData(contdatafile);
		Ncont = contdata->GetNsite();
		taxonset = contdata->GetTaxonSet();
		Ntaxa = taxonset->GetNtaxa();

		One = new Const<PosReal>(1);

		GammaDiag = new Gamma*[Ncont];
		diag = new Rvar<PosReal>*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(10);
		}
		priorOnSigmaZero = new RvarVec(diag,Ncont);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		sigma = new ConjugateInverseWishart(sigmaZero, Ncont+1);
		sigma->Sample();
		sigma->SetIdentity();

		data = new MultivariateNormal*[Ntaxa];
		for (int i =0; i<Ntaxa; i++)	{
			data[i] = new ConjugateMultivariateNormal(GetConjugateInverseWishart());
		}
		ClampData();

		cerr << "register\n";
		RootRegister(One);
		Register();

		MakeScheduler();

		if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~ConjugateBranchOmegaMultivariateModel() {}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{
		scheduler.Register(new ConjugateInverseWishartResampleMove(GetConjugateInverseWishart()),1,"conjugate sigma");
	}

};

