
#ifndef CONJUGATEPATH_H
#define CONJUGATEPATH_H

#include "RandomBranchSitePath.h"
#include "Conjugate.h"
#include "IIDConjugate.h"
#include "ConjugatePoisson.h"
#include "ConjugateMultinomial.h"
#include "GTRSubMatrix.h"

class LengthConjugatePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath {

	public:

	LengthConjugatePath(PhyloProcess* inprocess, Var<PosReal>* inlength) : RandomBranchSitePath(inprocess) {
		conj_length = dynamic_cast<PoissonSemiConjugate*>(inlength);
		if (conj_length)	{
			conjugate_up.insert(conj_length);
			Register(conj_length);
		}
	}

	bool isMatched()	{
		return (conj_length != 0);
	}

	protected:

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (conj_length && (parent == conj_length))	{
			conj_length->AddToShape(GetNsub());
			conj_length->AddToScale(LengthScaleStat());
		}
	}

	double LengthScaleStat()	{
		double total = 0;
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			total += -(*GetSubMatrix())(state,state) * GetRelativeTime(link);
			link = link->Next();
		}
		return total * GetRate();
	}

	private:

	PoissonSemiConjugate* conj_length;
};

class RateConjugatePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath {

	public:

	RateConjugatePath(PhyloProcess* inprocess, Var<PosReal>* inrate) : RandomBranchSitePath(inprocess)	{
		conj_rate = dynamic_cast<PoissonSemiConjugate*>(inrate);
		if (conj_rate)	{
			conjugate_up.insert(conj_rate);
			Register(conj_rate);
		}
	}

	bool isMatched()	{
		return (conj_rate != 0);
	}

	protected:

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (conj_rate && (parent == conj_rate))	{
			conj_rate->AddToShape(GetNsub());
			conj_rate->AddToScale(RateScaleStat());
		}
	}

	double RateScaleStat()	{
		double total = 0;
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			total += -(*GetSubMatrix())(state,state) * GetRelativeTime(link);
			link = link->Next();
		}
		return total * GetTime();
	}

	private:

	PoissonSemiConjugate* conj_rate;
};

class StationaryConjugatePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath {

	public:

	StationaryConjugatePath(PhyloProcess* inprocess, Var<Profile>* instat) : RandomBranchSitePath(inprocess) {
		conj_stat = dynamic_cast<GeneralizedMultinomialSemiConjugate*>(instat);
		if (conj_stat)	{
			conjugate_up.insert(conj_stat);
			Register(conj_stat);
		}
	}

	bool isMatched()	{
		return (conj_stat != 0);
	}

	protected:

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		if (conj_stat && (parent == conj_stat))	{
			int* counts = new int[GetNstate()];
			double* beta = new double[GetNstate()];
			ComputeStat(counts,beta);
			conj_stat->AddToShape(counts);
			conj_stat->AddToScale(beta);
			delete[] counts;
			delete[] beta;
		}
	}

	void ComputeStat(int* counts, double* beta)	{

		for (int k=0; k<GetNstate(); k++)	{
			counts[k] = 0;
			beta[k] = 0;
		}
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			if (isRoot())	{
				counts[state]++;
			}
			else	{
				if (link != Init()) {
					counts[state]++;
				}
				for (int k=0; k<GetNstate(); k++)	{
					if (k != state)	{
						if (GetStationary()[k] > 0)	{
							beta[k] += (*GetSubMatrix())(state,k) / GetStationary()[k] * GetAbsoluteTime(link) * GetRate();
						}
					}
				}
			}
			link = link->Next();
		}
	}

	private:

	GeneralizedMultinomialSemiConjugate* conj_stat;
};


class RelativeRateConjugatePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath {

	public:

	RelativeRateConjugatePath(PhyloProcess* inprocess, Var<PosRealVector>* inrelrate) : RandomBranchSitePath(inprocess)	{
		conj_relrate = dynamic_cast<IIDPoissonSemiConjugate*>(inrelrate);
		if (conj_relrate)	{
			conjugate_up.insert(conj_relrate);
			Register(conj_relrate);
		}
	}

	bool isMatched()	{
		return (conj_relrate != 0);
	}

	protected:

	int Nrr() {return GetNstate() * (GetNstate() - 1) / 2;}
	int rrindex(int i, int j) {return GTRSubMatrix::rrindex(i,j,GetNstate());}

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (conj_relrate && (parent == conj_relrate))	{
			int* counts = new int[Nrr()];
			double* beta = new double[Nrr()];
			ComputeStat(counts,beta);
			conj_relrate->AddToShape(counts);
			conj_relrate->AddToScale(beta);
			delete[] counts;
			delete[] beta;
		}
	}

	void ComputeStat(int* counts, double* beta)	{
		for (int k=0; k<Nrr() ; k++)	{
			counts[k] = 0;
			beta[k] = 0;
		}
		Plink* link = Init();
		while (link)	{
			int state = link->GetState();
			if (link != Last())	{
				counts[rrindex(link->GetState(), link->Next()->GetState())]++;
			}
			for (int k=0; k<GetNstate(); k++)	{
				if (k != state)	{
					beta[rrindex(state,k)] += GetStationary()[k] * GetAbsoluteTime(link) * GetRate();
				}
			}
			link = link->Next();
		}
	}

	private:

	IIDPoissonSemiConjugate* conj_relrate;
};


class ConjugateRandomBranchSitePath : public virtual RandomBranchSitePath, public LengthConjugatePath, public RateConjugatePath, public StationaryConjugatePath, public RelativeRateConjugatePath {

	public:

	ConjugateRandomBranchSitePath(PhyloProcess* inprocess, Var<PosReal>* inlength, Var<PosReal>* inrate, RandomSubMatrix* inmatrix, Var<Profile>* instat, Var<PosRealVector>* inrelrate)	:
			RandomBranchSitePath(inprocess, inlength, inrate, inmatrix, instat),
			LengthConjugatePath(inprocess,inlength),
			RateConjugatePath(inprocess,inrate),
			StationaryConjugatePath(inprocess,instat),
			RelativeRateConjugatePath(inprocess,inrelrate)
			{}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		LengthConjugatePath::AddSufficientStatistic(parent);
		RateConjugatePath::AddSufficientStatistic(parent);
		StationaryConjugatePath::AddSufficientStatistic(parent);
		RelativeRateConjugatePath::AddSufficientStatistic(parent);
	}

};

#endif
