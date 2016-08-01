
#ifndef BIPOISSON_H
#define BIPOISSON_H

#include "BaseType.h"
#include "Var.h"
#include "Conjugate.h"
#include "ConjugatePoisson.h"

class Bipoisson : public virtual Rvar<Int>	{

	public:

	Bipoisson(Var<PosReal>* inlength, Var<PosReal>* inrate)	{
		length = inlength;
		rate = inrate;
		Register(length);
		Register(rate);
	}

	void 		drawSample()	{
		double mu = length->val() * rate->val();
		double p = exp(mu) * Random::Uniform();
		int n = 0;
		double ratio = 1;
		double total = 1;
		while (total<p)	{
			n++;
			ratio *= mu / n;
			total += ratio;
		}
		value = n;
	}

	double 		ProposeMove(double tuning)	{
		bkvalue = value;
		flag = false;
		drawSample();
		return 0;
	}

	protected:

	double 		logProb()	{
		double mu = length->val() * rate->val();
		return -mu + value * log(mu) - Random::logGamma(value+1);
	}

	Var<PosReal>* length;
	Var<PosReal>* rate;
};



class ConjugateBipoisson : public ConjugateSampling<Int>, public Bipoisson	{


	public:

	ConjugateBipoisson(ConjugateGamma* inlength, ConjugateGamma* inrate)	: Bipoisson(inlength, inrate) {
		conjugate_up.insert(inlength);
		conjugate_up.insert(inrate);
	}
		
	~ConjugateBipoisson() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		if ((Var<PosReal>*) parent == length)	{
			ConjugateGamma* gammaprior = dynamic_cast<ConjugateGamma*>(parent);
			gammaprior->AddToShape(val());
			gammaprior->AddToScale(rate->val());

		}
		else if ((Var<PosReal>*) parent == rate)	{
			ConjugateGamma* gammaprior = dynamic_cast<ConjugateGamma*>(parent);
			gammaprior->AddToShape(val());
			gammaprior->AddToScale(length->val());
		}
		else	{
			cerr << (Var<PosReal>*) parent << '\t' << parent << '\t' << rate << '\t' << length << '\n';
			cerr << "error in ConjugateBipoisson::AddSuffStat\n";
			exit(1);
		}
	}
};

#endif

