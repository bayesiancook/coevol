
#pragma once

class BetaKappaMove : public MCUpdate, public Mnode	{

	public:

	BetaKappaMove(Rvar<Real>* inbeta, Rvar<Real>* inlogkappa1, double intuning, double ina)	{
		beta = inbeta;
		logkappa1 = inlogkappa1;
		tuning = intuning;
		a = ina;
		beta->Register(this);
		logkappa1->Register(this);
	}

	double Move(double tuning_modulator = 1)	{

		double acc = 1.0;
		if ((!beta->isClamped()) && (! logkappa1->isClamped()))	{
			
			Corrupt(true);
			double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
			beta->setval(beta->val() + u);
			logkappa1->setval(logkappa1->val() + a*u);
			double logratio = Update();
			acc = (log(Random::Uniform()) < logratio);
			if (! acc)	{
				Corrupt(false);
				Restore();
			}
		}
		return acc;

	}

	private:
	Rvar<Real>* beta;
	Rvar<Real>* logkappa1;
	double tuning;
	double a;
};

