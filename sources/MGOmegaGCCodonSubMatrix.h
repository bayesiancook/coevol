#ifndef MGOMEGAGCCODONSUBMATRIX_H
#define MGOMEGAGCCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
#include "IncompleteGamma.h"

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmegaGCCodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmegaGCCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double inpopsize, double inbgc, double inalpha, double inbeta, bool innormalise = false, int inN = 1000, double inmax = 100) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			popsize(inpopsize),
			bgc(inbgc),
			alpha(inalpha),
			beta(inbeta),
			N(inN),
			max(inmax),
			fixprobflag(false), quadflag(false)	{
				CreateQuadrature();
			}


	double			GetPopSize() {return popsize;}
	double			GetBGC() {return bgc;}
	double			GetAlpha() {return alpha;}
	double			GetBeta() {return beta;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void			CreateQuadrature();
	void			UpdateQuadrature();
	void			UpdateFixProbRatios();
	double			FixProbRatio(double x);
	double			AverageFixProbRatio(double s0);

	void 			ComputeArray(int state);

	void			SetPopSize(double inpopsize) {
					popsize = inpopsize;
					fixprobflag = false;
					quadflag = false;
				}

	void			SetBGC(double inbgc) {
					bgc= inbgc;
					fixprobflag = false;
				}

	void			SetAlpha(double inalpha) {
					alpha= inalpha;
					fixprobflag = false;
					quadflag = false;
				}

	void			SetBeta(double inbeta) {
					beta = inbeta;
					fixprobflag = false;
					quadflag = false;
				}

	// data members
	double popsize;
	double bgc;
	double alpha;
	double beta;

	// quadrature
	int N;
	double max;
	double* x;
	double* w;
	double* y;
	double* v;

	// fixation probabilities
	double PsGC;
	double PsAT;
	double PnGC;
	double PnAT;
	double Pn0;

	bool fixprobflag;
	bool quadflag;
};

class RandomMGOmegaGCCodonSubMatrix : public MGOmegaGCCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmegaGCCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inPopSize, Var<Real>* inBGC, Var<PosReal>* inAlpha, Var<PosReal>* inBeta, bool innormalise = false, int N = 1000, double max = 100) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			MGOmegaGCCodonSubMatrix(instatespace,inmatrix,inPopSize->val(), inBGC->val(), inAlpha->val(), inBeta->val(), innormalise, N, max) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			PopSize(inPopSize),
			BGC(inBGC),
			Alpha(inAlpha),
			Beta(inBeta)
			{


		Register(matrix);
		Register(PopSize);
		Register(BGC);
		Register(Alpha);
		Register(Beta);
		popsize = 1;
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetPopSize(PopSize->val());
		SetBGC(BGC->val());
		SetAlpha(Alpha->val());
		SetBeta(Beta->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* PopSize;
	Var<Real>* BGC;
	Var<PosReal>* Alpha;
	Var<PosReal>* Beta;
};

void MGOmegaGCCodonSubMatrix::CreateQuadrature()	{

	x = new double[N+1];
	w = new double[N+1];
	y = new double[N];
	v = new double[N];
}

void MGOmegaGCCodonSubMatrix::UpdateQuadrature()	{

	if (! quadflag)	{
		double lnga1 = Random::logGamma(alpha);
		double Beta = beta / popsize;
		x[0] = 0;
		w[0] = 0;
		for (int i=1; i<=N; i++) {
			x[i] = max * exp(log(((double) i) / N) / alpha);
			w[i] = IncompleteGamma(x[i]*Beta,alpha,lnga1);
		}
		for (int i=0; i<N; i++) {
			y[i] = (x[i+1] + x[i]) / 2;
			v[i] = w[i+1] - w[i];
		}
		quadflag = true;
	}
}

void MGOmegaGCCodonSubMatrix::UpdateFixProbRatios()	{

	if (! fixprobflag)	{
		PsGC = FixProbRatio(-popsize*bgc);
		PsAT = FixProbRatio(popsize*bgc);
		Pn0 = AverageFixProbRatio(0);
		PnGC = AverageFixProbRatio(-popsize*bgc);
		PnAT = AverageFixProbRatio(popsize*bgc);

		fixprobflag = true;
	}
}

inline double MGOmegaGCCodonSubMatrix::FixProbRatio(double x)	{

	if (fabs(x)<1e-10)	{
		return 1;
	}
	if (x > 100)	{
		return 0;
	}
	if (x < -100)	{
		return -x;
	}
	return x / (exp(x) -1);
}


double MGOmegaGCCodonSubMatrix::AverageFixProbRatio(double s0)	{

	double total = 0;
	for (int i=0; i<N; i++) {
		total += v[i] * FixProbRatio(y[i]+s0);
	}
	if (std::isinf(total))	{
		cerr << "inf prob in MGOmegaGCCodoSubMatrix\n";
		cerr << "rescaled bgc : " << s0 << '\n';
		cerr << '\n';
		for (int i=0; i<N; i++) {
			cerr << v[i] << '\t' << y[i] << '\t' << FixProbRatio(y[i]+s0) << '\n';
		}
		exit(1);
	}
	return total;
}

void MGOmegaGCCodonSubMatrix::ComputeArray(int i)	{

	// will perform actual updates only if flags have been set to false
	UpdateQuadrature();
	UpdateFixProbRatios();

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				Q[i][j] = (*NucMatrix)(a,b);

				bool initialgc = ((a == 1) || (a == 2));
				bool finalgc = ((b == 1) || (b == 2));
				// towards gc
				if ((!initialgc) && (finalgc))	{
					if (! Synonymous(i,j))	{
						Q[i][j] *= PnGC;
					}
					else	{
						Q[i][j] *= PsGC;
					}
				}
				// towards at
				else if ((initialgc) && (! finalgc))	{
					if (! Synonymous(i,j))	{
						Q[i][j] *= PnAT;
					}
					else	{
						Q[i][j] *= PsAT;
					}
				}
				// gc neutral
				else	{
					if (! Synonymous(i,j))	{
						Q[i][j] *= Pn0;
					}
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

#endif

/*
void MGOmegaGCCodonSubMatrix::ComputeStationary()	{

	// compute stationary probabilities
	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
		total += mStationary[i];
	}

	// renormalize stationary probabilities
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i]  /= total;
	}
}
*/
