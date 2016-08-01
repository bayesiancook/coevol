
#ifndef PRMOVE_H
#define PRMOVE_H

#include "Move.h"

class PosRealVectorMove : public MCUpdate	{

	public:

	PosRealVectorMove(Rvar<PosRealVector>* invar, double intuning, int inm);

	double Move(double tuning_modulator = 1);

	private:

	Rvar<PosRealVector>* var;
	double tuning;
	int m;
};

class PosRealVectorTranslationMove : public MCUpdate	{

	public:

	PosRealVectorTranslationMove(Rvar<PosRealVector>* invar, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Rvar<PosRealVector>* var;
	double tuning;
};


PosRealVectorMove::PosRealVectorMove(Rvar<PosRealVector>* invar, double intuning, int inm) : var(invar), tuning(intuning), m(inm) {}

double PosRealVectorMove::Move(double tuning_modulator)	{
	if (! var->isClamped())	{
		var->Corrupt(true);
		double logHastings = var->PosRealVector::ProposeMove(tuning * tuning_modulator,m);
		double deltaLogProb = var->Update();
		double logRatio = deltaLogProb + logHastings;
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			var->Corrupt(false);
			var->Restore();
		}
		return (double) accepted;
	}
	return 1;
}

PosRealVectorTranslationMove::PosRealVectorTranslationMove(Rvar<PosRealVector>* invar, double intuning) : var(invar), tuning(intuning) {}

double PosRealVectorTranslationMove::Move(double tuning_modulator)	{
	if (! var->isClamped())	{
		var->Corrupt(true);
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		var->PosRealVector::ScalarMultiplication(e);
		double logHastings = var->GetDim() * m;
		double deltaLogProb = var->Update();
		double logRatio = deltaLogProb + logHastings;
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			var->Corrupt(false);
			var->Restore();
		}
		return (double) accepted;
	}
	return 1;
}
#endif
