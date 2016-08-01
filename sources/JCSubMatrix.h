#ifndef JC_H
#define JC_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"

class JCSubMatrix : public virtual SubMatrix	{

	public:

				JCSubMatrix(int nstate, bool innormalise = false) : SubMatrix(nstate,innormalise) {}
				~JCSubMatrix() {};

	protected:

	void 			ComputeArray(int state);
	void 			ComputeStationary();

	// data members
	int Nstate;

};

void 	JCSubMatrix::ComputeStationary()	{

	for (int i=0; i<Nstate; i++)	{
		mStationary[i] = 1.0 / Nstate;
	}

}

void	JCSubMatrix::ComputeArray(int i)	{

	for (int j=0; j<Nstate; j++)	{
		Q[i][j] = mStationary[j];
	}
	Q[i][i] -= 1;
}

class JCRandomSubMatrix : public RandomSubMatrix, public JCSubMatrix  {

	public:
	JCRandomSubMatrix(Var<PosReal>* inroot, int dim, bool innormalise = false) : SubMatrix(dim,innormalise), RandomSubMatrix(dim,innormalise), JCSubMatrix(dim, innormalise) {

		root = inroot;
		Register(root);
		specialUpdate();
	}

	protected:

	void SetParameters()	{
	}

	private:
	Var<PosReal>* root;
};

#endif

