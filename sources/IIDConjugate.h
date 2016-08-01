#ifndef IIDCONJUGATE_H
#define IIDCONJUGATE_H

#include "BaseType.h"
#include "Var.h"
#include "Conjugate.h"

class IIDPoissonSemiConjugate: public virtual SemiConjugatePrior<PosRealVector> {

	public:

	IIDPoissonSemiConjugate()	{
		if (GetDim() == 0)	{
			cerr << "error in IIDPoissonSemiConjugate\n";
			cerr << "classes diamond-inheriting from IIDPoissonSemiConjugate should FIRST inherit from a Rvar<PosRealVector> type\n";
			cerr << "(order of declaration of inheritance matters here)\n";
			exit(1);
		}
		shapestat = new double[GetDim()];
		scalestat = new double[GetDim()];
		bkshapestat = new double[GetDim()];
		bkscalestat = new double[GetDim()];
	}

	virtual ~IIDPoissonSemiConjugate()	{
		delete[] shapestat;
		delete[] scalestat;
		delete[] bkshapestat;
		delete[] bkscalestat;
	}

	void ResetSufficientStatistic()	{
		for (int i=0; i<GetDim(); i++)	{
			shapestat[i] = 0;
			scalestat[i] = 0;
		}
	}

	void SaveSufficientStatistic()	{
		for (int i=0; i<GetDim(); i++)	{
			bkshapestat[i] = shapestat[i];
			bkscalestat[i] = scalestat[i];
		}
	}

	void RestoreSufficientStatistic()	{
		for (int i=0; i<GetDim(); i++)	{
			shapestat[i] = bkshapestat[i];
			scalestat[i] = bkscalestat[i];
		}
	}

	void AddToShape(int * in)	{
		for (int i=0; i<GetDim(); i++)	{
			shapestat[i] += in[i];
		}
	}

	void AddToScale(double* in)	{
		for (int i=0; i<GetDim(); i++)	{
			scalestat[i] += in[i];
		}
	}

	double SuffStatLogProb()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += shapestat[i] * log((*this)[i]) - scalestat[i] * (*this)[i];
		}
		return total;
	}

	protected:

	double* shapestat;
	double* scalestat;
	double* bkshapestat;
	double* bkscalestat;

};

class IIDConjugateGamma : public IIDGamma, public ConjugatePrior<PosRealVector> , IIDPoissonSemiConjugate {

	public:

	IIDConjugateGamma(int indim, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDGamma(indim, inalpha, inbeta) {}

	IIDConjugateGamma(int indim) : IIDGamma(indim) {}

	void GibbsResample()	{
		double a = alpha ? (double) alpha->val() : 1.0;
		double b = beta ? (double) beta->val() : 1.0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = Random::Gamma(a + shapestat[i], b + scalestat[i]);
			if ((*this)[i] < Gamma::GAMMAMIN)	{
				(*this)[i] = Gamma::GAMMAMIN;
			}
		}
	}

	double logProb()	{
		if (isActive())	{
			double a = alpha ? (double) alpha->val() : 1.0;
			double b = beta ? (double) beta->val() : 1.0;
			if (isIntegrated())	{
				double total = 0;
				for (int i=0; i<GetDim(); i++)	{
					total += - Random::logGamma(a) + a * log(b) + Random::logGamma(a + shapestat[i]) - (a + shapestat[i]) * log(b + scalestat[i]);
				}
				return total;
			}
			double total = 0;
			for (int i=0; i<GetDim(); i++)	{
				total += -Random::logGamma(a) + a * log(b) - (b + scalestat[i]) * (*this)[i] + (a + shapestat[i] - 1) * log((*this)[i]);
			}
			return total;
		}
		return IIDGamma::logProb();
	}
};

#endif
