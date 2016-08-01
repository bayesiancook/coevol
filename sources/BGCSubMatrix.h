
#ifndef BGCSUBMATRIX_H
#define BGCSUBMATRIX_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"
#include "NonRevSubMatrix.h"
#include "IncompleteGamma.h"

// #include <gsl/gsl_sf_zeta.h>

// const double bgcupperlimit = 0;
const double bgcupperlimit = 25.0;

class MutSubMatrix : public virtual SubMatrix	{

	public:

				MutSubMatrix(double inlambda, double inat2cg, double inat2gc, double inat2ta, double ingc2cg, bool innormalise) : SubMatrix(Nnuc,innormalise) {
					lambda = inlambda;
					at2cg = inat2cg;
					at2gc = inat2gc;
					at2ta = inat2ta;
					gc2cg = ingc2cg;
				}

				~MutSubMatrix() {};

	protected:

	void 			SetRelRates(double inat2cg, double inat2gc, double inat2ta, double ingc2cg) {
					at2cg = inat2cg;
					at2gc = inat2gc;
					at2ta = inat2ta;
					gc2cg = ingc2cg;
				}

	void 			SetLambda(double inlambda) {lambda = inlambda;}

	void ComputeArray(int i)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				double tmp = 1;

				if ((j == 0) || (j == 3))	{
					tmp *= lambda;
				}

				if ( ((i==0) && (j==1)) || ((i==1) && (j==0)) || ((i==3) && (j==2)) || ((i==2) && (j==3)) )	{
					tmp *= at2cg;
				}
				if ( ((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==3) && (j==1)) || ((i==1) && (j==3)) )	{
					tmp *= at2gc;
				}
				if ( ((i==0) && (j==3)) || ((i==3) && (j==0)) )	{
					tmp *= at2ta;
				}
				if ( ((i==1) && (j==2)) || ((i==2) && (j==1)) )	{
					tmp *= gc2cg;
				}

				Q[i][j] = tmp;
				total += tmp;
			}
		}
		Q[i][i] = - total;
	}

	void 			ComputeStationary()	{
		double z = 2 * (lambda + 1);
		mStationary[0] = lambda / z;
		mStationary[1] = 1.0 / z;
		mStationary[2] = 1.0 / z;
		mStationary[3] = lambda / z;
	}

	// data members
	double 			lambda;
	double			at2cg;
	double			at2gc;
	double			at2ta;
	double			gc2cg;
};

class RandomMutSubMatrix : public virtual RandomSubMatrix, public virtual MutSubMatrix  {

	public:
	RandomMutSubMatrix(Var<PosReal>* inlambda,Var<PosReal>* inAT2CG, Var<PosReal>* inAT2GC, Var<PosReal>* inAT2TA, Var<PosReal>* inGC2CG, bool innormalise, Var<PosReal>* inlambda2 = 0) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , MutSubMatrix(1.0,1.0,1.0,1.0,1.0,innormalise) {
		lambda = inlambda;
		AT2CG = inAT2CG;
		AT2GC = inAT2GC;
		AT2TA = inAT2TA;
		GC2CG = inGC2CG;
		lambda2= inlambda2;
		Register(lambda);
		Register(lambda2);
		Register(AT2CG);
		Register(AT2GC);
		Register(AT2TA);
		Register(GC2CG);
		specialUpdate();
	}

	protected:

	void SetParameters()	{
		if (lambda2)	{
			SetLambda(lambda->val() * lambda2->val());
		}
		else	{
			SetLambda(lambda->val());
		}
		SetRelRates(AT2CG->val(),AT2GC->val(),AT2TA->val(),GC2CG->val());
	}

	private:
	Var<PosReal>* lambda;
	Var<PosReal>* lambda2;
	Var<PosReal>* AT2CG;
	Var<PosReal>* AT2GC;
	Var<PosReal>* AT2TA;
	Var<PosReal>* GC2CG;
};

class NonRevMutSubMatrix : public virtual NonRevSubMatrix	{

	public:

				NonRevMutSubMatrix(double inat2cg, double inat2gc, double inat2ta, double incg2at, double incg2gc, double incg2ta, double* instat = 0, bool innormalise = false, int discn = 10) : SubMatrix(Nnuc,innormalise), NonRevSubMatrix(Nnuc,innormalise, discn) {
					at2cg = inat2cg;
					at2gc = inat2gc;
					at2ta = inat2ta;
					cg2at = incg2at;
					cg2gc = incg2gc;
					cg2ta = incg2ta;
				}

				~NonRevMutSubMatrix() {};

	protected:

	void 			SetRelRates(double inat2cg, double inat2gc, double inat2ta, double incg2at, double incg2gc, double incg2ta)	{
					at2cg = inat2cg;
					at2gc = inat2gc;
					at2ta = inat2ta;
					cg2at = incg2at;
					cg2gc = incg2gc;
					cg2ta = incg2ta;
				}

	void 			CopyStationary(const double* instat)	{
		if (instat)	{
			mStationary[0] = instat[0];
			mStationary[1] = instat[1];
			mStationary[2] = instat[2];
			mStationary[3] = instat[3];
		}
		else	{
			mStationary[0] = 0.25;
			mStationary[1] = 0.25;
			mStationary[2] = 0.25;
			mStationary[3] = 0.25;
		}
	}

	void ComputeArray(int i)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				double& tmp = Q[i][j];

				if ( ((i==0) && (j==1)) || ((i==3) && (j==2)) )	{
					tmp = at2cg;
				}
				if ( ((i==0) && (j==2)) || ((i==3) && (j==1)) ) {
					tmp = at2gc;
				}
				if ( ((i==0) && (j==3)) || ((i==3) && (j==0)) )	{
					tmp = at2ta;
				}
				if ( ((i==1) && (j==0)) || ((i==2) && (j==3)) )	{
					tmp = cg2at;
				}
				if ( ((i==1) && (j==2)) || ((i==2) && (j==1)) )	{
					tmp = cg2gc;
				}
				if ( ((i==1) && (j==3)) || ((i==2) && (j==0)) )	{
					tmp = cg2ta;
				}

				Q[i][j] = tmp;
				total += tmp;
			}
		}
		Q[i][i] = - total;
	}

	void 			ComputeStationary()	{}

	// data members
	double			at2cg;
	double			at2gc;
	double			at2ta;
	double			cg2at;
	double			cg2gc;
	double			cg2ta;
};

class RandomNonRevMutSubMatrix : public virtual RandomNonRevSubMatrix, public virtual NonRevMutSubMatrix  {

	public:
	RandomNonRevMutSubMatrix(Var<PosReal>* inAT2CG, Var<PosReal>* inAT2GC, Var<PosReal>* inAT2TA, Var<PosReal>* inCG2AT, Var<PosReal>* inCG2GC, Var<PosReal>* inCG2TA, bool innormalise = false, int indiscn = 10) :
		SubMatrix(Nnuc, innormalise),
		NonRevSubMatrix(Nnuc,innormalise,indiscn),
		RandomSubMatrix(Nnuc, innormalise) ,
		RandomNonRevSubMatrix(Nnuc, innormalise, indiscn),
		NonRevMutSubMatrix(1.0,1.0,1.0,1.0,1.0,1.0,0,innormalise,indiscn) {
		AT2CG = inAT2CG;
		AT2GC = inAT2GC;
		AT2TA = inAT2TA;
		CG2AT = inCG2AT;
		CG2GC = inCG2GC;
		CG2TA = inCG2TA;
		Register(AT2CG);
		Register(AT2GC);
		Register(AT2TA);
		Register(CG2AT);
		Register(CG2GC);
		Register(CG2TA);
		specialUpdate();
	}

	protected:

	void SetParameters()	{
		CopyStationary(0);
		SetRelRates(AT2CG->val(),AT2GC->val(),AT2TA->val(),CG2AT->val(),CG2GC->val(),CG2TA->val());
	}

	private:

	Var<PosReal>* AT2CG;
	Var<PosReal>* AT2GC;
	Var<PosReal>* AT2TA;
	Var<PosReal>* CG2AT;
	Var<PosReal>* CG2GC;
	Var<PosReal>* CG2TA;
};

class BGCMutSelSubMatrix : public virtual SubMatrix	{

	public:

				BGCMutSelSubMatrix(SubMatrix* inmutmatrix, double inmeanbgc, double inalpha, bool innormalise, bool inwithmean, bool inwithhotspot, int indiscgam) : SubMatrix(Nnuc,innormalise)	{
					mutmatrix = inmutmatrix;
					meanbgc = inmeanbgc;
					alpha = inalpha;
					withmean = inwithmean;
					withhotspot = inwithhotspot;
					discgam = indiscgam;
					if (discgam)	{
						CreateDiscreteGamma();
					}
				}

				~BGCMutSelSubMatrix() {
					if (discgam)	{
						DeleteDiscreteGamma();
					}
				}

	int			GetBGCOverflowCount() {return bgccount;}

	protected:

	void 			SetMeanBGC(double inbgc)	{meanbgc = inbgc;}
	void 			SetAlpha(double inalpha) {
					alpha = inalpha;
					if (discgam)	{
						RefreshDiscreteGamma();
					}
				}

	void			SetMutMatrix(SubMatrix* inmatrix) {mutmatrix = inmatrix;}


	void CreateDiscreteGamma()	{
		xx = new double[discgam];
		yy = new double[discgam];
		rr = new double[discgam];
	}

	void DeleteDiscreteGamma()	{
		delete[] xx;
		delete[] yy;
		delete[] rr;
	}

	void RefreshDiscreteGamma()	{
		double lg = Random::logGamma(alpha+1.0);
		for (int i=0; i<discgam; i++)	{
			xx[i] = PointGamma((i+1.0)/discgam,alpha,alpha);
		}
		for (int i=0; i<discgam; i++)	{
			yy[i] = IncompleteGamma(alpha*xx[i],alpha+1,lg);
		}
		yy[discgam-1] = 1.0;
		rr[0] = discgam * yy[0];
		for (int i=1; i<discgam; i++)	{
			rr[i] = discgam * (yy[i] - yy[i-1]);
		}
	}

	double DiscreteWS(double beta)	{
		double ret = 0;
		for (int i=0; i<discgam; i++)	{
			double z = yy[i] * alpha / beta;
			if (z > bgcupperlimit)	{
				z = bgcupperlimit;
			}
			if (fabs(z) > 1e-6)	{
				ret += z / (1 - exp(-z));
			}
			else	{
				ret += 1;
			}
		}
		ret /= discgam;
		return ret;
	}

	double DiscreteSW(double beta)	{
		double ret = 0;
		for (int i=0; i<discgam; i++)	{
			double z = yy[i] * alpha / beta;
			if (z > bgcupperlimit)	{
				z = bgcupperlimit;
			}
			if (fabs(z) > 1e-6)	{
				ret += z / (exp(z) - 1);
			}
			else	{
				ret += 1;
			}
		}
		ret /= discgam;
		return ret;
	}

	void ComputeArray(int i)	{
		if (withhotspot)	{
			double total = 0;
			for (int j=0; j<Nstate; j++)	{
				if (i!=j)	{
					double tmp = (*mutmatrix)(i,j);
					double b = meanbgc;
					if (b > bgcupperlimit)	{
						b = bgcupperlimit;
					}
					double a = alpha / (1 + alpha);
					// this is for root;
					if (! a)	{
						a = 1;
					}
					if (fabs(b) > 1e-6)	{
						if (((i == 1) || (i == 2)) && ((j == 0) || (j == 3)))	{
							tmp *= 1 -a + a * b / (exp(b) - 1);
						}
						else if (((i == 0) || (i == 3)) && ((j == 1) || (j == 2)))	{
							tmp *= 1 -a + a * b / (1 - exp(-b));
						}
					}
					Q[i][j] = tmp;
					total += tmp;
				}
			}
			Q[i][i] = - total;
		}
		else if (alpha != 0)	{
			if (discgam)	{
				double total = 0;
				double beta = withmean ? (alpha / meanbgc) : (1.0 / meanbgc);
				double a = alpha;
				if (alpha < 0.001)	{
					alpha = 0.001;
				}

				double sw = DiscreteSW(beta);
				double ws = DiscreteWS(beta);

				for (int j=0; j<Nstate; j++)	{
					if (i!=j)	{
						double tmp = (*mutmatrix)(i,j);
						if (((i == 1) || (i == 2)) && ((j == 0) || (j == 3)))	{
							tmp *= sw;
						}
						else if (((i == 0) || (i == 3)) && ((j == 1) || (j == 2)))	{
							tmp *= ws;
						}
						Q[i][j] = tmp;
						total += tmp;
					}
				}
				Q[i][i] = - total;
				alpha = a;
			}
			else	{
				double total = 0;
				// double beta = withmean ? (alpha / meanbgc) : (1.0 / meanbgc);
				double a = alpha;
				if (a < 0.001)	{
					a = 0.001;
				}

				cerr << "gsl deactivated\n";
				exit(1);
				double sw = 1; // a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta+1);
				double ws = 1; //  a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta);

				for (int j=0; j<Nstate; j++)	{
					if (i!=j)	{
						double tmp = (*mutmatrix)(i,j);
						if (((i == 1) || (i == 2)) && ((j == 0) || (j == 3)))	{
							tmp *= sw;
						}
						else if (((i == 0) || (i == 3)) && ((j == 1) || (j == 2)))	{
							tmp *= ws;
						}
						Q[i][j] = tmp;
						total += tmp;
					}
				}
				Q[i][i] = - total;
			}
		}
		else	{
			double total = 0;
			for (int j=0; j<Nstate; j++)	{
				if (i!=j)	{
					double tmp = (*mutmatrix)(i,j);
					double b = meanbgc;
					if (b > bgcupperlimit)	{
						b = bgcupperlimit;
					}
					if (fabs(b) > 1e-6)	{
						if (((i == 1) || (i == 2)) && ((j == 0) || (j == 3)))	{
							tmp *= b / (exp(b) - 1);
						}
						else if (((i == 0) || (i == 3)) && ((j == 1) || (j == 2)))	{
							tmp *= b / (1 - exp(-b));
						}
					}
					Q[i][j] = tmp;
					total += tmp;
				}
			}
			Q[i][i] = - total;
		}
	}

	/*
	void ComputeStationary()	{
		// const double* mutstat = mutmatrix->GetStationary();
		double b = meanbgc;
		double total = 0;
		for (int i=0; i<Nnuc; i++)	{
			if ((i==1) || (i==2))	{
				mStationary[i] = b;
				// mStationary[i] = mutstat[i] * b;
			}
			else	{
				mStationary[i] = 1.0;
				// mStationary[i] = mutstat[i];
			}
			total += mStationary[i];
		}
		for (int i=0; i<Nnuc; i++)	{
			mStationary[i] /= total;
		}
	}
	*/

	void ComputeStationary() {
		const double* mutstat = mutmatrix->GetStationary();
		if (withhotspot)	{
			double b = meanbgc;
			if (b > bgcupperlimit)	{
				b = bgcupperlimit;
				bgccount++;
			}
			// double a = alpha;
			double a = alpha / (1 + alpha);
			// this is for root;
			if (! a)	{
				a = 1;
			}
			if (fabs(b) > 1e-6)	{
				mStationary[0] = mutstat[0] * (1-a + a * b / (exp(b) - 1) );
				mStationary[1] = mutstat[1] * (1-a + a * b / (1 - exp(-b)));
				mStationary[2] = mutstat[2] * (1-a + a * b / (1 - exp(-b)));
				mStationary[3] = mutstat[3] * (1-a + a * b / (exp(b) - 1) );
			}
			else	{
				mStationary[0] = mutstat[0];
				mStationary[1] = mutstat[1];
				mStationary[2] = mutstat[2];
				mStationary[3] = mutstat[3];
			}
		}
		else if (alpha != 0)	{
			if (discgam)	{
				double beta = withmean ? (alpha / meanbgc) : (1.0 / meanbgc);
				double a = alpha;
				if (alpha < 0.001)	{
					alpha = 0.001;
				}
				double sw = DiscreteSW(beta);
				double ws = DiscreteWS(beta);

				mStationary[0] = mutstat[0] * sw;
				mStationary[1] = mutstat[1] * ws;
				mStationary[2] = mutstat[2] * ws;
				mStationary[3] = mutstat[3] * sw;
				alpha = a;
			}
			else	{
			// 	double beta = withmean ? (alpha / meanbgc) : (1.0 / meanbgc);
				double a = alpha;
				if (a < 0.001)	{
					a = 0.001;
				}
				cerr << "gsl deactivated\n";
				exit(1);
				double sw = 1; // a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta+1);
				double ws = 1; //  a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta);

				mStationary[0] = mutstat[0] * sw;
				mStationary[1] = mutstat[1] * ws;
				mStationary[2] = mutstat[2] * ws;
				mStationary[3] = mutstat[3] * sw;
			}
		}
		else	{
			double b = meanbgc;
			if (b > bgcupperlimit)	{
				b = bgcupperlimit;
				bgccount++;
			}
			mStationary[0] = mutstat[0];
			mStationary[1] = mutstat[1] * exp(b);
			mStationary[2] = mutstat[2] * exp(b);
			mStationary[3] = mutstat[3];
		}
		double total = 0;
		for (int i=0; i<Nnuc; i++)	{
			total += mStationary[i];
		}
		for (int i=0; i<Nnuc; i++)	{
			mStationary[i] /= total;
		}
	}

	// data members
	double 			meanbgc;
	double			alpha;
	SubMatrix*		mutmatrix;

	static int		bgccount;
	bool 			withmean;
	bool			withhotspot;

	double* xx;
	double* yy;
	double* rr;
	int discgam;

};

class RandomBGCMutSelSubMatrix : public virtual RandomSubMatrix, public virtual BGCMutSelSubMatrix  {

	public:

	RandomBGCMutSelSubMatrix(RandomSubMatrix* inmutmatrix, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, bool innormalise, int indiscgam) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , BGCMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,true,false,indiscgam)	{
		mutmatrix = inmutmatrix;
		meanbgc = inmeanbgc;
		alpha = inalpha;
		Register(mutmatrix);
		Register(meanbgc);
		if (alpha)	{
			Register(alpha);
		}
		specialUpdate();
	}

	Var<PosReal>* GetMeanBGC() {return meanbgc;}
	Var<PosReal>* GetAlpha() {return alpha;}

	RandomSubMatrix* GetMutMatrix() {return mutmatrix;}

	protected:

	virtual void SetParameters()	{
		SetMutMatrix(mutmatrix);
		SetMeanBGC(meanbgc->val());
		double tmp = 0;
		if (alpha)	{
			tmp = alpha->val();
		}
		SetAlpha(tmp);
	}

	private:
	RandomSubMatrix* mutmatrix;
	Var<PosReal>* meanbgc;
	Var<PosReal>* alpha;
};


class BGCNonRevMutSelSubMatrix : public virtual NonRevSubMatrix, public virtual BGCMutSelSubMatrix {

	public:

				BGCNonRevMutSelSubMatrix(NonRevSubMatrix* inmutmatrix, double inmeanbgc, double inalpha, bool innormalise, bool inwithmean, bool inwithhotspot, int discgam, int discn) :
					SubMatrix(Nnuc,innormalise),
					NonRevSubMatrix(Nnuc,innormalise, discn),
					BGCMutSelSubMatrix(inmutmatrix,inmeanbgc,inalpha,innormalise,inwithmean,inwithhotspot,discgam)	{
					}

				~BGCNonRevMutSelSubMatrix() {};

	/*
	virtual void ComputeExponential(double l)	{
		cerr << "bgc expo : " << l << '\n';
		cerr << "expo flag : " << expoflag << '\n';
		for (int i=0; i<GetNstate(); i++)	{
			cerr << flagarray[i] << '\t';
		}
		cerr << '\n';

		NonRevSubMatrix::ComputeExponential(l);
		cerr << "bgc expo ok\n";

		cerr << "mut : " << '\n';
		mutmatrix->ToStream(cerr);
		cerr << '\n';
		exit(1);
	}
	*/

};

class RandomBGCNonRevMutSelSubMatrix : public virtual RandomNonRevSubMatrix, public virtual BGCNonRevMutSelSubMatrix, public virtual RandomBGCMutSelSubMatrix  {

	public:
	RandomBGCNonRevMutSelSubMatrix(RandomNonRevSubMatrix* inmutmatrix, Var<PosReal>* inlength, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, bool innormalise, int discgam, int discn) :
			SubMatrix(Nnuc,innormalise),
			NonRevSubMatrix(Nnuc,innormalise,discn),
			RandomSubMatrix(Nnuc,innormalise),
			RandomNonRevSubMatrix(Nnuc,innormalise,discn) ,
			BGCMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,true,false,discgam),
			BGCNonRevMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,true,false,discgam,discn),
			RandomBGCMutSelSubMatrix(inmutmatrix,inmeanbgc,inalpha,innormalise,discgam)	{

				Length = inlength;
				Register(Length);
				// SetParameters();
				specialUpdate();
	}

	protected:

	void SetParameters()	{
		if (Length)	{
			SetLength(Length->val());
		}
		RandomBGCMutSelSubMatrix::SetParameters();
	}

	private:

	Var<PosReal>* Length;

};
/*
class RandomNonRevBGCMutSelSubMatrix : public RandomNonRevSubMatrix, public BGCMutSelSubMatrix {

	public:
	RandomNonRevBGCMutSelSubMatrix(RandomNonRevSubMatrix* inmutmatrix, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, Var<Profile>* inStat, Var<PosReal>* inLength, bool innormalise = false) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , BGCMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,true)	{
		mutmatrix = inmutmatrix;
		meanbgc = inmeanbgc;
		alpha = inalpha;
		Register(mutmatrix);
		Register(meanbgc);
		if (alpha)	{
			Register(alpha);
		}
		specialUpdate();
	}

	Var<PosReal>* GetMeanBGC() {return meanbgc;}
	Var<PosReal>* GetAlpha() {return alpha;}

	RandomSubMatrix* GetMutMatrix() {return mutmatrix;}

	protected:

	void SetParameters()	{
		SetMutMatrix(mutmatrix);
		SetMeanBGC(meanbgc->val());
		double tmp = 0;
		if (alpha)	{
			tmp = alpha->val();
		}
		SetAlpha(tmp);
	}

	private:
	RandomSubMatrix* mutmatrix;
	Var<PosReal>* meanbgc;
	Var<PosReal>* alpha;
};
*/

class RandomAlphaBGCMutSelSubMatrix : public RandomSubMatrix, public BGCMutSelSubMatrix  {

	public:
	RandomAlphaBGCMutSelSubMatrix(RandomSubMatrix* inmutmatrix, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, Var<PosReal>* inalphaoffset, bool innormalise = false) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , BGCMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,false,false,0)	{
		mutmatrix = inmutmatrix;
		meanbgc = inmeanbgc;
		alpha = inalpha;
		alphaoffset = inalphaoffset;
		Register(mutmatrix);
		Register(meanbgc);
		if (alpha)	{
			Register(alpha);
		}
		if (alphaoffset)	{
			Register(alphaoffset);
		}
		specialUpdate();
	}

	Var<PosReal>* GetMeanBGC() {return meanbgc;}
	Var<PosReal>* GetAlpha() {return alpha;}
	Var<PosReal>* GetAlphaOffset() {return alphaoffset;}

	RandomSubMatrix* GetMutMatrix() {return mutmatrix;}

	protected:

	void SetParameters()	{
		SetMutMatrix(mutmatrix);
		SetMeanBGC(meanbgc->val());
		double tmp = 0;
		if (alpha)	{
			tmp = alpha->val() * alphaoffset->val();
		}
		if (tmp > bgcupperlimit)	{
			tmp = bgcupperlimit;
		}
		SetAlpha(tmp);
	}

	private:
	RandomSubMatrix* mutmatrix;
	Var<PosReal>* meanbgc;
	Var<PosReal>* alpha;
	Var<PosReal>* alphaoffset;
};

class RandomHotSpotBGCMutSelSubMatrix : public RandomSubMatrix, public BGCMutSelSubMatrix  {

	public:
	RandomHotSpotBGCMutSelSubMatrix(RandomSubMatrix* inmutmatrix, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, Var<PosReal>* inalphaoffset, bool innormalise = false) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , BGCMutSelSubMatrix(inmutmatrix,1.0,0,innormalise,true,true,0)	{
		mutmatrix = inmutmatrix;
		meanbgc = inmeanbgc;
		alpha = inalpha;
		alphaoffset = inalphaoffset;
		Register(mutmatrix);
		Register(meanbgc);
		if (alpha)	{
			Register(alpha);
		}
		if (alphaoffset)	{
			Register(alphaoffset);
		}
		specialUpdate();
	}

	Var<PosReal>* GetMeanBGC() {return meanbgc;}
	Var<PosReal>* GetAlpha() {return alpha;}
	Var<PosReal>* GetAlphaOffset() {return alphaoffset;}

	RandomSubMatrix* GetMutMatrix() {return mutmatrix;}

	protected:

	void SetParameters()	{
		SetMutMatrix(mutmatrix);
		SetMeanBGC(meanbgc->val());
		double tmp = 0;
		if (alpha)	{
			tmp = alpha->val() * alphaoffset->val();
		}
		if (tmp > bgcupperlimit)	{
			tmp = bgcupperlimit;
		}
		SetAlpha(tmp);
	}

	private:
	RandomSubMatrix* mutmatrix;
	Var<PosReal>* meanbgc;
	Var<PosReal>* alpha;
	Var<PosReal>* alphaoffset;
};

/*
class BGCMutSelSubMatrix : public virtual SubMatrix	{

	public:

				BGCMutSelSubMatrix(SubMatrix* inmutmatrix, double inmeanbgc, bool innormalise = false) : SubMatrix(Nnuc,innormalise)	{
					mutmatrix = inmutmatrix;
					meanbgc = inmeanbgc;
				}

				~BGCMutSelSubMatrix() {};

	int			GetBGCOverflowCount() {return bgccount;}

	protected:

	void 			SetMeanBGC(double inbgc)	{meanbgc = inbgc;}
	void			SetMutMatrix(SubMatrix* inmatrix) {mutmatrix = inmatrix;}


	void ComputeArray(int i)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			if (i!=j)	{
				double tmp = (*mutmatrix)(i,j);
				double b = meanbgc;
				if (b > bgcupperlimit)	{
					b = bgcupperlimit;
				}
				if (fabs(b) > 1e-6)	{
					if (((i == 1) || (i == 2)) && ((j == 0) || (j == 3)))	{
						tmp *= b / (exp(b) - 1);
					}
					else if (((i == 0) || (i == 3)) && ((j == 1) || (j == 2)))	{
						tmp *= b / (1 - exp(-b));
					}
				}
				Q[i][j] = tmp;
				total += tmp;
			}
		}
		Q[i][i] = - total;
	}

	void ComputeStationary() {
		const double* mutstat = mutmatrix->GetStationary();
		double b = meanbgc;
		if (b > bgcupperlimit)	{
			// cerr << b << '\n';
			b = bgcupperlimit;
			bgccount++;
		}
		mStationary[0] = mutstat[0];
		mStationary[1] = mutstat[1] * exp(b);
		mStationary[2] = mutstat[2] * exp(b);
		mStationary[3] = mutstat[3];

		double total = 0;
		for (int i=0; i<Nnuc; i++)	{
			total += mStationary[i];
		}
		for (int i=0; i<Nnuc; i++)	{
			mStationary[i] /= total;
		}
	}

	// data members
	double 			meanbgc;
	SubMatrix*		mutmatrix;

	static int		bgccount;
};

class RandomBGCMutSelSubMatrix : public RandomSubMatrix, public BGCMutSelSubMatrix  {

	public:
	RandomBGCMutSelSubMatrix(RandomSubMatrix* inmutmatrix, Var<PosReal>* inmeanbgc, bool innormalise = false) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , BGCMutSelSubMatrix(inmutmatrix,1.0,innormalise)	{
		mutmatrix = inmutmatrix;
		meanbgc = inmeanbgc;
		Register(mutmatrix);
		Register(meanbgc);
		specialUpdate();
	}

	Var<PosReal>* GetMeanBGC() {return meanbgc;}

	RandomSubMatrix* GetMutMatrix() {return mutmatrix;}

	protected:

	void SetParameters()	{
		SetMutMatrix(mutmatrix);
		SetMeanBGC(meanbgc->val());
	}

	private:
	RandomSubMatrix* mutmatrix;
	Var<PosReal>* meanbgc;
};
*/

#endif
