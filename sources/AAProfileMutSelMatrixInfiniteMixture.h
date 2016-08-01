#ifndef AAPROFILEMUTSELMATRIXINFINITEMIXTURE_H
#define AAPROFILEMUTSELMATRIXINFINITEMIXTURE_H


class AAProfileMutSelMatrixInfiniteMixture : public MatrixInfiniteMixture<Profile>	{

	public:

	AAProfileMutSelMatrixInfiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, Var<PosReal>* inNeff) :
			MatrixInfiniteMixture<Profile>(insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
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
		AAProfileMutSelMixtureRandomMatrix* temp = new AAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
		return temp;
		//return new AAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<PosReal>* Neff;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};

class AAProfileMutSelMatInfMixValMove : public MCUpdate	{

	public:

	AAProfileMutSelMatInfMixValMove(AAProfileMutSelMatrixInfiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	AAProfileMutSelMatrixInfiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};


#endif
