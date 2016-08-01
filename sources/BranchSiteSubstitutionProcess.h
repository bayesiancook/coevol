#ifndef SUBPROCESS_H
#define SUBPROCESS_H

#include "Random.h"
#include "SubMatrix.h"
#include "Var.h"

class BranchSiteSubstitutionProcess	{

	// a matrix
	// a rate
	// a time duration (branch length)
	// a state space?

	public:

	virtual ~BranchSiteSubstitutionProcess() {}

	virtual double GetTime() = 0;
	virtual double GetRate() = 0;
	virtual const double* GetStationary() = 0;
	virtual SubMatrix* GetSubMatrix() = 0;

	virtual int GetNstate() {return GetSubMatrix()->GetNstate();}

	void BackwardPropagate(const double* down, double* up)	{
		GetSubMatrix()->BackwardPropagate(down,up,GetTime()*GetRate());
	}
	void ForwardPropagate(const double* up, double* down)	{
		GetSubMatrix()->ForwardPropagate(up,down,GetTime()*GetRate());
	}

	int DrawFiniteTime(int state);
	int DrawStationary();

	void GetFiniteTimeTransitionProb(int state, double* down);

	double GetFiniteTimeTransitionProb(int stateup, int statedown);

	double DrawWaitingTime(int state);
	double DrawWaitingTimeGivenAtLeastOne(int state, double time);
	int DrawOneStep(int state);

	double WaitingTimeLogProb(int state, double time);
	double ReducedWaitingTimeLogProb(int state, double time);
	double OneStepLogProb(int stateup, int statedown);
	double StationaryLogProb(int state);
	double ReducedOneStepLogProb(int stateup, int statedown);

	int DrawUniformizedTransition(int state, int statedown, int n);
	int DrawUniformizedSubstitutionNumber(int state, int statedown);
	bool CheckUniformizedSubstitutionNumber(int state, int statedown);

};


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


inline double BranchSiteSubstitutionProcess::WaitingTimeLogProb(int state, double time)	{

	double q = -(*GetSubMatrix())(state,state) * GetRate();
	return log(q) - q*time;
}

inline double BranchSiteSubstitutionProcess::ReducedWaitingTimeLogProb(int state, double time)	{

	double q = -(*GetSubMatrix())(state,state) * GetRate();
	return  - q*time;
}

inline double BranchSiteSubstitutionProcess::OneStepLogProb(int stateup, int statedown)	{

	const double* row = GetSubMatrix()->GetRow(stateup);
	cerr << "in BranchSiteSubstitutionProcess::OneStepLogProb : check code\n";
	exit(1);
	return log(row[statedown]) - log(row[stateup]);
}

inline double BranchSiteSubstitutionProcess::ReducedOneStepLogProb(int stateup, int statedown)	{

	const double* row = GetSubMatrix()->GetRow(stateup);
	return log(row[statedown]);
}

inline double BranchSiteSubstitutionProcess::DrawWaitingTime(int state)	{

	double q = -(*GetSubMatrix())(state,state) * GetRate();
	return  -log (1 - Random::Uniform()) / q;
}

inline double BranchSiteSubstitutionProcess::DrawWaitingTimeGivenAtLeastOne(int state, double totaltime)	{

	double q = -(*GetSubMatrix())(state,state) * GetRate();
	return -log(1 - Random::Uniform() * (1 - exp(-q * totaltime))) / q;

}

inline int BranchSiteSubstitutionProcess::DrawOneStep(int state)	{

	const double* row = GetSubMatrix()->GetRow(state);
	double p = -row[state] * Random::Uniform();
	int k = -1;
	double tot = 0;
	do	{
		k++;
		if (k != state)	{
			tot += row[k];
		}
	}
	while ((k<GetNstate()) && (tot < p));
	if (k == GetNstate())	{
		cerr << "error in BranchSiteSubstitutionProcess::DrawOneStep\n";
		exit(1);
	}
	return k;
}

inline void BranchSiteSubstitutionProcess::GetFiniteTimeTransitionProb(int state, double* p)	{

	double* p1 = new double[GetNstate()];
	for (int k=0; k<GetNstate(); k++)	{
		p1[k] = 0;
	}
	p1[state] = 1;
	// GetSubMatrix()->FiniteTime(state,p,GetTime()*GetRate());
	GetSubMatrix()->ForwardPropagate(p1,p,GetTime()*GetRate());
	double tot = 0;
	for (int k=0; k<GetNstate(); k++)	{
		tot += p[k];
	}
	if (fabs(1 - tot) > 1e-5)	{
		cerr << "error in forward propagate: normalization : " << tot << '\t' << fabs(1 - tot) << '\n';
		cerr << GetTime() << '\t' << GetRate() << '\n';
		GetSubMatrix()->ToStream(cerr);
		exit(1);
	}
	delete[] p1;
}

inline double BranchSiteSubstitutionProcess::GetFiniteTimeTransitionProb(int stateup, int statedown)	{

	double* p = new double[GetNstate()];
	double* p1 = new double[GetNstate()];
	for (int k=0; k<GetNstate(); k++)	{
		p1[k] = 0;
	}
	p1[stateup] = 1;
	// GetSubMatrix()->FiniteTime(stateup,p,GetTime()*GetRate());
	GetSubMatrix()->ForwardPropagate(p1,p,GetTime()*GetRate());
	// ForwardPropagate(p1,p);
	double tot = 0;
	for (int k=0; k<GetNstate(); k++)	{
		tot += p[k];
	}
	if (fabs(1 - tot) > 1e-5)	{
		cerr << "error in forward propagate: normalization : " << tot << '\t' << fabs(1 - tot) << '\n';
		GetSubMatrix()->ToStream(cerr);
		exit(1);
	}
	double ret = p[statedown];
	delete[] p1;
	delete[] p;
	return ret;
}

inline double BranchSiteSubstitutionProcess::StationaryLogProb(int state)	{

	return log(GetStationary()[state]);
}

inline int BranchSiteSubstitutionProcess::DrawStationary()	{

	const double* stat = GetStationary();
	double p = Random::Uniform();
	int k = -1;
	double tot = 0;
	do	{
		k++;
		tot += stat[k];
	}
	while ((k<GetNstate()) && (tot < p));
	if (k == GetNstate())	{
		cerr << "erroro in BranchSiteSubstitutionProcess::DrawOneStep\n";
		exit(1);
	}
	return k;
}


inline int BranchSiteSubstitutionProcess::DrawFiniteTime(int state)	{

	double* p = new double[GetNstate()];
	GetFiniteTimeTransitionProb(state,p);
	int newstate = Random::DrawFromDiscreteDistribution(p,GetNstate());
	delete[] p;
	return newstate;
}


inline bool BranchSiteSubstitutionProcess::CheckUniformizedSubstitutionNumber(int stateup, int statedown)	{

	double Z = GetFiniteTimeTransitionProb(stateup,statedown);
	double efflength = GetRate() * GetTime();
	double mu = GetSubMatrix()->GetUniformizationMu();
	double fact = exp(- efflength * mu);
	int m = 0;
	double total = (stateup==statedown) * fact;
	while ((m<SubMatrix::UniSubNmax)) 	{
		m++;
		fact *= mu * efflength / m;
		total += GetSubMatrix()->Power(m,stateup,statedown) * fact;
	}
	if (fabs(total-Z)>1e-12)	{
		cerr << "error in BranchSiteSubstitutionProcess::DrawUniformizedSubstitutionNumber: normalising constant\n";
		cerr << total << '\t' << Z << '\n';
		cerr << mu << '\n';
		cerr << m << '\n';
		cerr << stateup << '\t' << statedown << '\n';
		throw;
	}
	return true;
}

inline int BranchSiteSubstitutionProcess::DrawUniformizedSubstitutionNumber(int stateup, int statedown)	{

	double Z = GetFiniteTimeTransitionProb(stateup,statedown);
	double efflength = GetRate() * GetTime();
	double mu = GetSubMatrix()->GetUniformizationMu();
	double fact = exp(- efflength * mu);
	int m = 0;
	double total = (stateup==statedown) * fact;
	double q = Random::Uniform() * Z;

	while ((m<SubMatrix::UniSubNmax) && (total < q)) 	{
		m++;
		fact *= mu * efflength / m;
		total += GetSubMatrix()->Power(m,stateup,statedown) * fact;
		if ((total-Z)>1e-12)	{
			cerr << "error in BranchSiteSubstitutionProcess::DrawUniformizedSubstitutionNumber: normalising constant\n";
			cerr << total << '\t' << Z << '\t' << total - Z << '\n';
			/*
			cerr << mu << '\n';
			cerr << m << '\n';
			cerr << stateup << '\t' << statedown << '\n';

			GetSubMatrix()->ToStream(cerr);
			GetSubMatrix()->CheckReversibility();
			throw;
			*/
		}
	}

	/*
	if (m >= SubMatrix::UniSubNmax)	{
		cerr << "error in BranchSiteSubstitutionProcess::DrawUniformizedSubstitutionNumber: overflow\n";
		throw;
	}
	if (m > mParam->ObservedUniSubNmax)	{
		mParam->ObservedUniSubNmax = m;
	}
	*/
	return m;
}


inline int BranchSiteSubstitutionProcess::DrawUniformizedTransition(int state, int statedown, int n)	{

	double* p = new double[GetNstate()];
	double tot = 0;
	for (int l=0; l<GetNstate(); l++)	{
		tot += GetSubMatrix()->Power(1,state,l) * GetSubMatrix()->Power(n,l,statedown);
		p[l] = tot;
	}

	double s = tot * Random::Uniform();
	int k = 0;
	while ((k<GetNstate()) && (s > p[k]))	{
		k++;
	}
	if (k == GetNstate())	{
		cerr << "error in BranchSiteSubstitutionProcess::DrawUniformizedTransition: overflow\n";
		throw;
	}
	delete[] p;
	return k;
}

#endif // SUBPROCESS_H
