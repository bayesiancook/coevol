
#ifndef HKYSUBMATRIX_H
#define HKYSUBMATRIX_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"

class HKYSubMatrix : public virtual SubMatrix	{

	public:

				HKYSubMatrix(const double kappa, const double* stat, const double* tsrr, const double* tvrr, bool innormalise) : SubMatrix(Nnuc,innormalise) {
					SetKappa(kappa);
					SetTsRR(tsrr);
					SetTvRR(tvrr);
					CopyStationary(stat);
				}

				~HKYSubMatrix() {};

	protected:

	void  CopyStationary(const double* instat)	{
		for (int k=0; k<Nstate; k++)	{
			mStationary[k] = instat[k];
		}
	}

	void			SetKappa(double inkappa) {kappa = inkappa;}
	void			SetTsRR(const double* intsrr) {tsrr = intsrr;}
	void			SetTvRR(const double* intvrr) {tvrr = intvrr;}

	private:

	void 			ComputeArray(int state);
	void 			ComputeStationary() {}

	double kappa;
	const double* tsrr;
	const double* tvrr;
};

class HKYRandomSubMatrix : public RandomSubMatrix, public HKYSubMatrix  {

	public:
	HKYRandomSubMatrix(Var<PosReal>* inKappa, Var<Profile>* instat, Var<Profile>* inTsRR, Var<Profile>* inTvRR, bool innormalise) :
		SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , HKYSubMatrix(1.0,instat->val().GetArray(), 0,0,innormalise), Kappa(inKappa), stat(instat), TsRR(inTsRR), TvRR(inTvRR)	{
		if (stat->GetDim() != Nnuc)	{
			cerr << "error in HKY : equilibrium frequency profile should be of dimension 4, and not " << stat->GetDim() << '\n';
			throw;
		}
		if (Kappa)	{
			Register(Kappa);
		}
		if (TsRR)	{
			if (TsRR->GetDim() != 2)	{
				cerr << "error in hky: transition rates should be of dimension 2\n";
				exit(1);
			}
			Register(TsRR);
		}
		if (TvRR)	{
			if (TvRR->GetDim() != 4)	{
				cerr << "error in hky: transversion rates should be of dimension 4\n";
				exit(1);
			}
			Register(TvRR);
		}
		Register(stat);
		specialUpdate();
	}

	Var<Profile>* GetRandomStationaries() {return stat;}

	protected:

	void SetParameters()	{
		if (Kappa)	{
			SetKappa(Kappa->val());
		}
		else	{
			SetKappa(2);
		}
		if (TsRR)	{
			SetTsRR(TsRR->val().GetArray());
		}
		if (TvRR)	{
			SetTvRR(TvRR->val().GetArray());
		}
		CopyStationary(stat->val().GetArray());
	}

	private:
	Var<PosReal>* Kappa;
	Var<Profile>* stat;
	Var<Profile>* TsRR;
	Var<Profile>* TvRR;
};


// ---------------------------------------------------------------------------
//		 ComputeArray
// ---------------------------------------------------------------------------

void	HKYSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<Nstate; j++)	{
		if (i!=j)	{
			Q[i][j] = mStationary[j];
			if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
				Q[i][j] *= kappa;
				if (tsrr)	{
					if (((i==0) && (j==2)) || ((i==2) && (j==0)))	{
						Q[i][j] *= tsrr[0];
					}
					else	{
						Q[i][j] *= tsrr[1];
					}
				}
			}
			else	{
				if (tvrr)	{
					if (((i==0) && (j==1)) || ((i==1) && (j==0)))	{
						Q[i][j] *= tvrr[0];
					}
					else if (((i==0) && (j==3)) || ((i==3) && (j==0)))	{
						Q[i][j] *= tvrr[1];
					}
					else if(((i==2) && (j==1)) || ((i==1) && (j==2)))	{
						Q[i][j] *= tvrr[2];
					}
					else	{
						Q[i][j] *= tvrr[3];
					}
				}
			}
			total += Q[i][j];
		}
	}

	// should always ensure that the diagonal entry of the matrix Q[i][i] is such that
	// the sum over all entries of the row is equal to 0
	Q[i][i] = - total;
}

#endif
