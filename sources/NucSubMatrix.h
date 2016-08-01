
#ifndef NUCSUBMATRIX_H
#define NUCSUBMATRIX_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"

class Nuc2X2SubMatrix : public virtual SubMatrix	{

	public:

				Nuc2X2SubMatrix(double tstv, const double* stat, const double* rr, bool innormalise) : SubMatrix(Nnuc,innormalise) {
					SetTsTv(tstv);
					SetRR(rr);
					CopyStationary(stat);
				}

				~Nuc2X2SubMatrix() {};

	int			GetNRelativeRate() {return Nrr;}
	double 			RelativeRate(int i, int j) {return rr[rrindex(i,j,GetNstate())];}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	protected:

	void  CopyStationary(const double* instat)	{
		for (int k=0; k<Nstate; k++)	{
			mStationary[k] = instat[k];
		}
	}

	void			SetTsTv(double intstv) {tstv= intstv;}
	void			SetRR(const double* inrr) {rr = inrr;}


	virtual void 		ComputeArray(int state);
	void 			ComputeStationary() {}

	double tstv;
	const double* rr;
	static const int Nrr = 6;
};

class Nuc3X3SubMatrix : public virtual Nuc2X2SubMatrix {

	public:

				Nuc3X3SubMatrix(double tstv, double tvgc, const double* stat, const double* rr, bool innormalise) : SubMatrix(Nnuc,innormalise), Nuc2X2SubMatrix(tstv,stat,rr,innormalise) {
					SetTvGC(tvgc);
				}

				~Nuc3X3SubMatrix() {};

	protected:

	void			SetTvGC(double intvgc) {tvgc = intvgc;}
	void			SetRR(const double* inrr) {rr = inrr;}

	private:

	void 			ComputeArray(int state);

	double tvgc;
};

class RandomNuc2X2SubMatrix : public RandomSubMatrix, public Nuc2X2SubMatrix  {

	public:
	RandomNuc2X2SubMatrix(Var<PosReal>* inTsTv, Var<Profile>* instat, Var<Profile>* inRR, bool innormalise) :
		SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , Nuc2X2SubMatrix(1.0,instat->val().GetArray(),0,innormalise), TsTv(inTsTv), stat(instat), RR(inRR)	{
		if (stat->GetDim() != Nnuc)	{
			cerr << "error in Nuc2X2 : equilibrium frequency profile should be of dimension 4, and not " << stat->GetDim() << '\n';
			throw;
		}
		if (TsTv)	{
			Register(TsTv);
		}
		if (RR)	{
			Register(RR);
		}
		Register(stat);
		specialUpdate();
	}

	Var<Profile>* GetRandomStationaries() {return stat;}

	protected:

	void SetParameters()	{
		if (TsTv)	{
			SetTsTv(TsTv->val());
		}
		else	{
			SetTsTv(2);
		}
		if (RR)	{
			SetRR(RR->val().GetArray());
		}
		CopyStationary(stat->val().GetArray());
	}

	Var<PosReal>* TsTv;
	Var<Profile>* stat;
	Var<Profile>* RR;
};

class RandomNuc3X3SubMatrix : public RandomSubMatrix, public Nuc3X3SubMatrix  {

	public:
	RandomNuc3X3SubMatrix(Var<PosReal>* inTsTv, Var<PosReal>* inTvGC, Var<Profile>* instat, Var<Profile>* inRR, bool innormalise) :
		SubMatrix(Nnuc, innormalise),
		Nuc2X2SubMatrix(1.0,instat->val().GetArray(),0,innormalise),
		RandomSubMatrix(Nnuc, innormalise),
		Nuc3X3SubMatrix(1.0,1.0,instat->val().GetArray(),0,innormalise),
		TsTv(inTsTv),
		stat(instat),
		RR(inRR),
		TvGC(inTvGC)
		{
		if (TsTv)	{
			Register(TsTv);
		}
		if (TvGC)	{
			Register(TvGC);
		}
		if (RR)	{
			Register(RR);
		}
		Register(stat);
		specialUpdate();
	}

	Var<Profile>* GetRandomStationaries() {return stat;}

	protected:

	void SetParameters()	{
		if (TsTv)	{
			SetTsTv(TsTv->val());
		}
		if (TvGC)	{
			SetTvGC(TvGC->val());
		}
		else	{
			SetTvGC(1);
		}
		if (RR)	{
			SetRR(RR->val().GetArray());
		}
		CopyStationary(stat->val().GetArray());
	}

	private:
	Var<PosReal>* TsTv;
	Var<Profile>* stat;
	Var<Profile>* RR;
	Var<PosReal>* TvGC;
};


// ---------------------------------------------------------------------------
//		 ComputeArray
// ---------------------------------------------------------------------------

void	Nuc2X2SubMatrix::ComputeArray(int i)	{

	if (rr)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				Q[i][j] = RelativeRate(i,j) * mStationary[j];
				if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
					Q[i][j] *= tstv;
				}
				total += Q[i][j];
			}
		}
		Q[i][i] = - total;
	}
	else	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] = mStationary[j];
			if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
				Q[i][j] *= tstv;
			}
		}
		Q[i][i] -= 1;
	}
}

void	Nuc3X3SubMatrix::ComputeArray(int i)	{

	if (rr)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				Q[i][j] = RelativeRate(i,j) * mStationary[j];
				if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
					Q[i][j] *= tstv;
				}
				else if (((i==0) && (j==1)) || ((i==1) && (j==0)) || ((i==2) && (j==3)) || ((i==3) && (j==2)))	{
					Q[i][j] *= tvgc;
				}
				total += Q[i][j];
			}
		}
		Q[i][i] = - total;
	}
	else	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] = mStationary[j];
			if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
				Q[i][j] *= tstv;
			}
			else if (((i==0) && (j==1)) || ((i==1) && (j==0)) || ((i==2) && (j==3)) || ((i==3) && (j==2)))	{
				Q[i][j] *= tvgc;
			}
		}
		Q[i][i] -= 1;
	}
}

#endif
