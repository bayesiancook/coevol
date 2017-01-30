
#ifndef EMPSUB_H
#define EMPSUB_H

#include "SubMatrix.h"
#include "BranchSitePath.h"
#include <vector>


// in charge of gathering statistics from a series of random branch site paths
class EmpiricalSubMatrix : public virtual SubMatrix	{

	public:

	EmpiricalSubMatrix(int Nstate, bool innormalise = false) : SubMatrix(Nstate, innormalise)	{

		paircounts = new int*[GetNstate()];
		for (int i=0; i<GetNstate(); i++)	{
			paircounts[i] = new int[GetNstate()];
		}
		statecounts = new int[GetNstate()];
	}

	virtual ~EmpiricalSubMatrix()	{
		delete[] statecounts;
		for (int i=0; i<GetNstate(); i++)	{
			delete[] paircounts[i];
		}
		delete[] paircounts;
	}

	void AddPath(BranchSitePath* inpath)	{
		pathlist.push_back(inpath);
	}

	void RefreshStatistics()	{
		for (int i=0; i<GetNstate(); i++)	{
			statecounts[i] = 0;
			for (int j=0; j<GetNstate(); j++)	{
				paircounts[i][j] = 0;
			}
		}
		for (unsigned int i=0; i<pathlist.size(); i++)	{
			pathlist[i]->AddCounts(paircounts,statecounts);
		}
		CorruptMatrix();
	}

	protected:

	int** paircounts;
	int* statecounts;

	vector<BranchSitePath*> pathlist;
	// a list of BranchSitePaths
};

class EmpiricalGTRSubMatrix : public EmpiricalSubMatrix	{

	public :

	EmpiricalGTRSubMatrix(int Nstate, bool innormalise = false) : SubMatrix(Nstate,innormalise), EmpiricalSubMatrix(Nstate,innormalise)	{
		rr = new double[GetNrr()];
		stat = new double[GetNstate()];
	}

	~EmpiricalGTRSubMatrix()	{
		delete[] rr;
		delete[] stat;
	}

	void ComputeArray(int i)	{
		if (! statflag)	{
			UpdateStationary();
			UpdateRelRates();
		}
		double total = 0;
		for (int j=0; j<GetNstate(); j++)	{
			if (i!=j)	{
				Q[i][j] = rr[rrindex(i,j,GetNstate())] * stat[j];
				total += Q[i][j];
			}
		}

		// should always ensure that the diagonal entry of the matrix Q[i][i] is such that
		// the sum over all entries of the row is equal to 0
		Q[i][i] = - total;
	}

	void ComputeStationary()	{
		double total = 0;
		for (int i=0; i<GetNstate(); i++)	{
			stat[i] = statecounts[i];
			total += stat[i];
		}
		for (int i=0; i<GetNstate(); i++)	{
			stat[i] /= total;
		}
	}

	void UpdateRelRates()	{
		for (int i=0; i<GetNrr(); i++)	{
			rr[i] = 0;
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				if (i != j)	{
					rr[rrindex(i,j,GetNstate())] += paircounts[i][j];
				}
			}
		}
		double total = 0;
		for (int i=0; i<GetNrr(); i++)	{
			total += rr[i];
		}
		for (int i=0; i<GetNrr(); i++)	{
			 rr[i] /= total;
		}
	}

	protected:

	int GetNrr()	{
		return Nstate * (Nstate-1) / 2;
	}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	double* rr;
	double* stat;

};

/*
class EmpiricalMGOmegaCodonSubMatrix : public EmpiricalSubMatrix, public MGOmegaCodonSubMatrix	{

};
*/


#endif
