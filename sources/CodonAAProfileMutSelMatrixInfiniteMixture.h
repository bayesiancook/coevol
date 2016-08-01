#ifndef CODONAAPROFILEMUTSELMATRIXINFINITEMIXTURE_H
#define CODONAAPROFILEMUTSELMATRIXINFINITEMIXTURE_H


class CodonAAProfileMutSelMatrixInfiniteMixture : public MatrixInfiniteMixture<Profile>	{

	public:

	CodonAAProfileMutSelMatrixInfiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, Var<Profile>* incodonprofile, Var<PosReal>* inNeff) :
			MatrixInfiniteMixture<Profile>(insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		codonprofile = incodonprofile;
		Neff = inNeff;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};


	protected:

	Dirichlet* GetDirichlet(int k)    {
		Dirichlet* temp = dynamic_cast<Dirichlet*>(GetComponent(k));
		if (!temp)	{
			cerr << "null pointer...\n";
			exit(1);
		}
                return temp;
	}

	MixtureRandomMatrix<Profile>* CreateComponent(int k)	{
		CodonAAProfileMutSelMixtureRandomMatrix* temp = new CodonAAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix,codonprofile, Neff);
		return temp;
		//return new AAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<PosReal>* Neff;
	Var<Profile>* codonprofile;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};

class CodonAAProfileMutSelMatInfMixValMove : public MCUpdate	{

	public:

	CodonAAProfileMutSelMatInfMixValMove(CodonAAProfileMutSelMatrixInfiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	CodonAAProfileMutSelMatrixInfiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};


#endif
