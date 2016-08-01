#ifndef DOLLOMAT_H
#define DOLLOMAT_H

#include "RandomSubMatrix.h"


class DolloSubMatrix: public virtual SubMatrix	{

	public:

	DolloSubMatrix(double inrate, bool innormalise) : SubMatrix(2,innormalise)	{
		rate = inrate;
	}

	protected:

	void SetRate(double inrate)	{
		rate = inrate;
	}

	void ComputeArray(int state)	{

		if (state == 0)	{
			Q[0][1] = 0;
			Q[0][0] = 0;
		}
		else	{
			Q[1][1] = -rate;
			Q[1][0] = rate;
		}
	}

	void ComputeStationary()	{
		mStationary[0] = 0;
		mStationary[1] = 1;
	}

	virtual void 		BackwardPropagate(const double* down, double* up, double length)	{
		double expo = exp(-rate * length);
		up[0] = down[0];
		up[1] = expo * down[1] + (1-expo) * down[0];
	}

	virtual void 		ForwardPropagate(const double* up, double* down, double length)	{
		double expo = exp(-rate * length);
		down[0] = up[0] + up[1] * (1 - expo);
		down[1] = up[1] * expo;
	}

	double rate;

};


class RandomDolloSubMatrix : public virtual DolloSubMatrix, public RandomSubMatrix	{

	public:

	RandomDolloSubMatrix(Var<PosReal>* inRate, bool innormalise) : SubMatrix(2,innormalise), DolloSubMatrix(inRate->val(),innormalise), RandomSubMatrix(2,innormalise)	{
		Rate = inRate;
		Register(Rate);
		specialUpdate();
	}

	void SetParameters()	{
		SetRate(Rate->val());
	}

	protected:

	Var<PosReal>* Rate;

};


#endif

