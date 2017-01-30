

#ifndef CONJUGATEMULTINOMIAL_H
#define CONJUGATEMULTINOMIAL_H

#include "Conjugate.h"
#include "RandomTypes.h"


class MultinomialSemiConjugate : public virtual SemiConjugatePrior<Profile> {

	public:

	MultinomialSemiConjugate(int dimension) {
		count = new int[dimension];
		bkcount = new int[dimension];
	}

	virtual ~MultinomialSemiConjugate() {
		delete[] count;
		delete[] bkcount;
	}

	void ResetSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] = 0;
		}
	}

	void SaveSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			bkcount[k] = count[k];
		}
	}

	void RestoreSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] = bkcount[k];
		}
	}

	void AddCounts(const int* incount)	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] += incount[k];
		}
	}

	double SuffStatLogProb()	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += count[k] * log((*this)[k]);
		}
		return total;
	}

	protected:

	int* count;
	int* bkcount;
};


class ConjugateMultinomial : public ConjugateSampling<IntVector>, public Multinomial	{

	public:

	ConjugateMultinomial(MultinomialSemiConjugate* inprobarray, int N) : Multinomial(inprobarray,N) {
		conjugate_up.insert(inprobarray);
	}

	~ConjugateMultinomial() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		MultinomialSemiConjugate* p = dynamic_cast<MultinomialSemiConjugate*>(parent);
		if (! p)	{
			cerr << "cast error in ConjugateMultinomial::AddSuffStat\n";
			exit(1);
		}
		p->AddCounts(Multinomial::val().GetArray());
	}

};

class ConjugateDirichlet : public ConjugatePrior<Profile>, public MultinomialSemiConjugate, public Dirichlet	{

	public:

	ConjugateDirichlet(int dimension) : MultinomialSemiConjugate(dimension), Dirichlet (dimension) {}

	ConjugateDirichlet(Var<Profile>* incenter, Var<PosReal>* inconcentration) : MultinomialSemiConjugate(incenter->GetDim()), Dirichlet(incenter,inconcentration) {}

	void GibbsResample()	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			double tmp = 0;
			if (concentration)	{
				tmp = Random::sGamma(concentration->val() * (center->val()[k] + count[k]));
			}
			else	{
				tmp = Random::sGamma(1 + count[k]);
			}
			if (tmp < Profile::MIN)	{
				tmp = Profile::MIN;
			}
			(*this)[k] = tmp;
			total += tmp;
		}
		for (int k=0; k<GetDim(); k++)	{
			(*this)[k] /= total;
		}
	}

	double logProb()	{
		if (isActive())	{
			if (isIntegrated())	{
				double totalcount = 0;
				for (int i=0; i<GetDim(); i++)	{
					totalcount += count[i];
				}
				double total = 0;
				if (concentration)	{
					total += Random::logGamma(concentration->val()) - Random::logGamma(concentration->val() + totalcount);
					for (int i=0; i<GetDim(); i++)	{
						total -= Random::logGamma(concentration->val() * center->val()[i]) - Random::logGamma(concentration->val() * center->val()[i] + count[i]);
					}
				}
				else	{
					total += Random::logGamma(GetDim()) - Random::logGamma(GetDim() + totalcount);
					for (int i=0; i<GetDim(); i++)	{
						total += Random::logGamma(1 + count[i]);
					}
				}
				return total;
			}
			double total = 0;
			if (concentration)	{
				total = Random::logGamma(concentration->val());
				for (int i=0; i<GetDim(); i++)	{
					total -=  Random::logGamma(concentration->val() * center->val()[i]) - (concentration->val() * center->val()[i] + count[i] - 1) * log((*this)[i]);
				}
			}
			else	{
				total = Random::logGamma(GetDim());
				for (int i=0; i<GetDim(); i++)	{
					total -= count[i] * log((*this)[i]);
				}
			}
			return total;
		}
		return Dirichlet::logProb();
	}
};

class GeneralizedMultinomialSemiConjugate : public virtual SemiConjugatePrior<Profile> {

	public:

	GeneralizedMultinomialSemiConjugate(int dimension) {
		count = new int[dimension];
		bkcount = new int[dimension];
		betastat = new double[dimension];
		bkbetastat = new double[dimension];
	}

	~GeneralizedMultinomialSemiConjugate() {
		delete[] count;
		delete[] bkcount;
		delete[] betastat;
		delete[] bkbetastat;
	}

	void ResetSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] = 0;
			betastat[k] = 0;
		}
	}

	void SaveSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			bkcount[k] = count[k];
			bkbetastat[k] = betastat[k];
		}
	}

	void RestoreSufficientStatistic()	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] = bkcount[k];
			betastat[k] = bkbetastat[k];
		}
	}

	void AddToShape(const int* incount)	{
		for (int k=0; k<GetDim(); k++)	{
			count[k] += incount[k];
		}
	}

	void AddToScale(const double* inbetastat)	{
		for (int k=0; k<GetDim(); k++)	{
			betastat[k] += inbetastat[k];
		}
	}

	double SuffStatLogProb()	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += count[k] * log((*this)[k]) - betastat[k] * (*this)[k];
		}
		return total;
	}

	protected:

	int* count;
	int* bkcount;
	double* betastat;
	double* bkbetastat;

};

class SemiConjugateDirichlet : public GeneralizedMultinomialSemiConjugate, public Dirichlet	{

	public:

	SemiConjugateDirichlet(int dimension) : GeneralizedMultinomialSemiConjugate(dimension), Dirichlet (dimension) {}

	SemiConjugateDirichlet(Var<Profile>* incenter, Var<PosReal>* inconcentration) : GeneralizedMultinomialSemiConjugate(incenter->GetDim()), Dirichlet(incenter,inconcentration) {}

	double logProb()	{
		if (isActive())	{
			double total = 0;
			if (concentration)	{
				total = Random::logGamma(concentration->val());
				for (int i=0; i<GetDim(); i++)	{
					total +=  -Random::logGamma(concentration->val() * center->val()[i]) - betastat[i] * val()[i] + (concentration->val() * center->val()[i] + count[i] - 1) * log(val()[i]);
				}
			}
			else	{
				total = Random::logGamma(GetDim());
				for (int i=0; i<GetDim(); i++)	{
					total += - betastat[i] * val()[i] + count[i] * log(val()[i]);
				}
			}
			return total;
		}
		return Dirichlet::logProb();
	}
};


#endif
