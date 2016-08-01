

#ifndef RADCONSSUBMATRIX_H
#define RADCONSSUBMATRIX_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"

class RadConsSubMatrix : public virtual SubMatrix	{

	public:

	RadConsSubMatrix(int nstate, const double inr, const double* stat, bool innormalise = false) : SubMatrix(nstate,normalise) {
		r = inr;
		CopyStationary(stat);
	}

	~RadConsSubMatrix() {};

	protected:

	// make a copy of the entries (not of the pointer)
	void CopyStationary(const double* instat)	{
		for (int k=0; k<Nstate; k++)	{
			mStationary[k] = instat[k];
		}
	}

	// copy the pointer
	void SetRatio(double inr) {r = inr;}

	protected:

	void ComputeArray(int i)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i != j)	{
				Q[i][j] = mStationary[j];
				if (Dayhoff6Table[i] != Dayhoff6Table[j])	{
					Q[i][j] *= r;
				}
				total += Q[i][j];
			}
		}
		Q[i][i] = - total;
	}

	void ComputeStationary() {}


	// data members
	double r;
};

class RadConsRandomSubMatrix : public RandomSubMatrix, public RadConsSubMatrix  {

	public:
	RadConsRandomSubMatrix(Var<PosReal>* inratio, Var<Profile>* instat, bool innormalise = false) : SubMatrix(instat->GetDim(), innormalise), RandomSubMatrix(instat->GetDim(), innormalise) , RadConsSubMatrix(instat->GetDim(),inratio->val(),instat->val().GetArray(), innormalise)	{
		ratio = inratio;
		stat = instat;
		Register(ratio);
		Register(stat);
		specialUpdate();
	}

	Var<Profile>* GetRandomStationaries() {return stat;}

	protected:

	void SetParameters()	{
		SetRatio(ratio->val());
		CopyStationary(stat->val().GetArray());
	}

	private:
	Var<PosReal>* ratio;
	Var<Profile>* stat;
};

#endif
