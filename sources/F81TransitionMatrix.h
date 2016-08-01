
#include "RandomSubMatrix.h"

class F81TransitionMatrix : public virtual TransitionMatrix	{

	public:

	F81TransitionMatrix(int inNstate, const double* stat, double inlength)	: TransitionMatrix(inNstate)	{
		length = inlength;
		for (int i=0; i<GetNstate(); i++)	{
			stationary[i] = stat[i];
		}
	}

	protected:

	virtual double GetLength() = 0;

	void ComputeArray()	{

		if (length != GetLength())	{
			cerr << "error : non matching lengths\n";
			cerr << length << '\t' << GetLength() << '\n';
			exit(1);
		}
		double expo = exp(-length);
		// cerr << length << '\n';
		for (int i=0; i<Nstate; i++)	{
			// double total = 0;
			for (int j=0; j<Nstate; j++)	{
				if (i!=j)	{
					R[i][j] = (1 - expo) * stationary[j];
					// total += R[i][j];
				}
			}
			R[i][i] = (1 - expo) * stationary[i] + expo;
		}
	}

	void ComputeStationary()	{
	}

	double length;
};

class RandomF81TransitionMatrix : public virtual F81TransitionMatrix, public virtual RandomTransitionMatrix {

	public:

	RandomF81TransitionMatrix(Var<Profile>* inRandomStat, Var<PosReal>* inRandomLength) :
		TransitionMatrix(inRandomStat->GetDim()),
		F81TransitionMatrix(inRandomStat->GetDim(), inRandomStat->GetArray(), 1.0)	{
		RandomStat = inRandomStat;
		RandomLength = inRandomLength;
		Register(RandomStat);
		Register(RandomLength);
		specialUpdate();
	}

	~RandomF81TransitionMatrix() {}

	double GetLength()	{
		if (! RandomLength)	{
			return 0;
		}
		return RandomLength->val();
	}

	protected:

	void SetParameters()	{
		if (RandomLength)	{
			length = RandomLength->val();
		}
		else	{
			length = 0;
		}
		for (int i=0; i<Nstate; i++)	{
			stationary[i] = (*RandomStat)[i];
		}
	}

	Var<Profile>* RandomStat;
	Var<PosReal>* RandomLength;

};

