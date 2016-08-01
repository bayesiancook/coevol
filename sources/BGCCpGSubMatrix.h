#ifndef BGCCPG_H
#define BGCCPG_H

#include "BGCSubMatrix.h"


class BGCCpGMutSelSubMatrix : public virtual SubMatrix	{

	public:

				BGCCpGMutSelSubMatrix(SubMatrix* inmutmatrix, SubMatrix* inctmutmatrix, double incpgrate, double inmeanbgc, double inalpha, bool innormalise, bool inwithmean, bool inwithhotspot, int indiscgam) : SubMatrix(64,innormalise)	{
					mutmatrix = inmutmatrix;
					ctmutmatrix = inctmutmatrix;
					cpgrate = incpgrate;
					meanbgc = inmeanbgc;
					alpha = inalpha;
					withmean = inwithmean;
					withhotspot = inwithhotspot;
					discgam = indiscgam;
					if (discgam)	{
						CreateDiscreteGamma();
					}
				}

				~BGCCpGMutSelSubMatrix() {
					if (discgam)	{
						DeleteDiscreteGamma();
					}
				}

	int			GetBGCOverflowCount() {return bgccount;}

	protected:

	void			SetCpGRate(double incpgrate)	{cpgrate = incpgrate;}
	void 			SetMeanBGC(double inbgc)	{meanbgc = inbgc;}
	void 			SetAlpha(double inalpha) {
					alpha = inalpha;
					if (discgam)	{
						RefreshDiscreteGamma();
					}
				}

	void			SetMutMatrix(SubMatrix* inmatrix) {mutmatrix = inmatrix;}
	void			SetContextMutMatrix(SubMatrix* inmatrix) {ctmutmatrix = inmatrix;}


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
		double b = meanbgc;
		if (b > bgcupperlimit)	{
			b = bgcupperlimit;
		}
		double fixpws = 1;
		double fixpsw = 1;
		if (fabs(b) > 1e-6)	{
			fixpsw = b / (exp(b) - 1);
			fixpws = b / (1 - exp(-b));
		}
		if (withhotspot)	{
			cerr << "hot spot + CpG not implemented yet\n";
			exit(1);
		}
		else if (alpha != 0)	{
			cerr << "alpha + CpG not implemented yet\n";
			/*
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
				double beta = withmean ? (alpha / meanbgc) : (1.0 / meanbgc);
				double a = alpha;
				if (a < 0.001)	{
					a = 0.001;
				}

				double sw = a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta+1);
				double ws =  a * exp(a * log(beta)) * gsl_sf_hzeta(a+1,beta);

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
			*/
		}
		else	{
			double* q = Q[i];
			for (int j = 0; j<Nstate; j++)	{
				(*q++) = 0;
				// Q[i][j] = 0;
			}


			int i1 = i / 16;
			int i2 = (i - 16*i1) / 4;
			int i3 = i - 16*i1 - 4*i2;

			for (int j1 = 0; j1<Nnuc; j1++)	{
				if (j1 != i1)	{
					double tmp = (*ctmutmatrix)(i1,j1);
					int j = 16*j1 + 4*i2 + i3;
					Q[i][j] = tmp;
				}
			}

			for (int j3 = 0; j3<Nnuc; j3++)	{
				if (j3 != i3)	{
					double tmp = (*ctmutmatrix)(i3,j3);
					int j = 16*i1 + 4*i2 + j3;
					Q[i][j] = tmp;
				}
			}

			for (int j2 = 0; j2<Nnuc; j2++)	{
				if (j2 != i2)	{
					int j = 16*i1 + 4*j2 + i3;
					double tmp = (*mutmatrix)(i2,j2);
					if (((i1 == 1) && (i2 == 2) && (j2 == 0)) || ((i3 == 2) && (i2 == 1) && (j2 == 3)))	{
						tmp += cpgrate;
					}
					/*
					double b = meanbgc;
					if (b > bgcupperlimit)	{
						b = bgcupperlimit;
					}
					if (fabs(b) > 1e-6)	{
					*/
					if (((i2 == 1) || (i2 == 2)) && ((j2 == 0) || (j2 == 3)))	{
						// tmp *= b / (exp(b) - 1);
						tmp *= fixpsw;
					}
					else if (((i2 == 0) || (i2 == 3)) && ((j2 == 1) || (j2 == 2)))	{
						// tmp *= b / (1 - exp(-b));
						tmp *= fixpws;
					}
					// }
					Q[i][j] = tmp;
				}
			}

			double total = 0;
			q = Q[i];
			for (int j = 0; j<Nstate; j++)	{
				total += (*q++);
				// total += Q[i][j];
			}
			Q[i][i] = - total;
		}
	}

	void ComputeStationary()	{
		/*
		const double* mutstat1 = ctmutmatrix->GetStationary();
		const double* mutstat2 = mutmatrix->GetStationary();
		const double* mutstat3 = ctmutmatrix->GetStationary();
		*/
		double b = meanbgc;
		double total = 0;
		for (int i1=0; i1<Nnuc; i1++)	{
			for (int i2 = 0; i2 <Nnuc; i2++)	{
				for (int i3 = 0; i3 <Nnuc; i3++)	{
					int i = 16*i1 + 4*i2 + i3;
					if ((i2 == 1) || (i2 == 2))	{
						mStationary[i] = b;
						// mStationary[i] = mutstat1[i1] * mutstat2[i2] * mutstat3[i3] * b;
					}
					else	{
						mStationary[i] = 1.0;
						// mStationary[i] = mutstat1[i1] * mutstat2[i2] * mutstat3[i3];
					}
					total += mStationary[i];
				}
			}
		}
		for (int i=0; i<Nstate; i++)	{
			mStationary[i] /= total;
		}
	}

	// data members
	double 			meanbgc;
	double			alpha;
	double			cpgrate;
	SubMatrix*		mutmatrix;
	SubMatrix*		ctmutmatrix;

	static int		bgccount;
	bool 			withmean;
	bool			withhotspot;

	double* xx;
	double* yy;
	double* rr;
	int discgam;

};

class RandomBGCCpGMutSelSubMatrix : public virtual RandomSubMatrix, public virtual BGCCpGMutSelSubMatrix  {

	public:

	RandomBGCCpGMutSelSubMatrix(RandomSubMatrix* inmutmatrix, RandomSubMatrix* inctmutmatrix, Var<PosReal>* incpgrate, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, bool innormalise, int indiscgam) : SubMatrix(64, innormalise), RandomSubMatrix(64, innormalise) , BGCCpGMutSelSubMatrix(inmutmatrix,inctmutmatrix,1.0,1.0,0,innormalise,true,false,indiscgam)	{
		mutmatrix = inmutmatrix;
		ctmutmatrix = inctmutmatrix;
		cpgrate = incpgrate;
		meanbgc = inmeanbgc;
		alpha = inalpha;
		Register(mutmatrix);
		Register(ctmutmatrix);
		Register(cpgrate);
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
		SetContextMutMatrix(ctmutmatrix);
		SetCpGRate(cpgrate->val());
		SetMeanBGC(meanbgc->val());
		double tmp = 0;
		if (alpha)	{
			tmp = alpha->val();
		}
		SetAlpha(tmp);
	}

	private:
	RandomSubMatrix* mutmatrix;
	RandomSubMatrix* ctmutmatrix;
	Var<PosReal>* cpgrate;
	Var<PosReal>* meanbgc;
	Var<PosReal>* alpha;
};


class BGCCpGNonRevMutSelSubMatrix : public virtual NonRevSubMatrix, public virtual BGCCpGMutSelSubMatrix {

	public:

				BGCCpGNonRevMutSelSubMatrix(NonRevSubMatrix* inmutmatrix, NonRevSubMatrix* inctmutmatrix, double incpgrate, double inmeanbgc, double inalpha, bool innormalise, bool inwithmean, bool inwithhotspot, int discgam, int discn) :
					SubMatrix(64,innormalise),
					NonRevSubMatrix(64,innormalise, discn),
					BGCCpGMutSelSubMatrix(inmutmatrix,inctmutmatrix,incpgrate,inmeanbgc,inalpha,innormalise,inwithmean,inwithhotspot,discgam)	{
					}

				~BGCCpGNonRevMutSelSubMatrix() {};

};

class RandomBGCCpGNonRevMutSelSubMatrix : public virtual RandomNonRevSubMatrix, public virtual BGCCpGNonRevMutSelSubMatrix, public virtual RandomBGCCpGMutSelSubMatrix  {

	public:
	RandomBGCCpGNonRevMutSelSubMatrix(RandomNonRevSubMatrix* inmutmatrix, RandomNonRevSubMatrix* inctmutmatrix, Var<PosReal>* incpgrate, Var<PosReal>* inlength, Var<PosReal>* inmeanbgc, Var<PosReal>* inalpha, bool innormalise, int discgam, int discn) :
			SubMatrix(64,innormalise),
			NonRevSubMatrix(64,innormalise,discn),
			RandomSubMatrix(64,innormalise),
			RandomNonRevSubMatrix(64,innormalise,discn) ,
			BGCCpGMutSelSubMatrix(inmutmatrix,inctmutmatrix,1.0,1.0,0,innormalise,true,false,discgam),
			BGCCpGNonRevMutSelSubMatrix(inmutmatrix,inctmutmatrix,1.0,1.0,0,innormalise,true,false,discgam,discn),
			RandomBGCCpGMutSelSubMatrix(inmutmatrix,inctmutmatrix,incpgrate,inmeanbgc,inalpha,innormalise,discgam)	{

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
		RandomBGCCpGMutSelSubMatrix::SetParameters();
	}

	private:

	Var<PosReal>* Length;

};


#endif

