
#include "IID.h"
#include "BaseType.h"
#include "RandomTypes.h"
#include "Move.h"
#include "Normal.h"


class DirichletIIDArray : public IIDArray<Profile>	{

	public:

	DirichletIIDArray(int insize, Var<Profile>* incenter, Var<PosReal>* inconcentration) : IIDArray<Profile>(insize)	{
		center = incenter;
		concentration = inconcentration;
		Create();
	}

	Dirichlet* GetDirichletVal(int site)	{
		return dynamic_cast<Dirichlet*>(GetVal(site));
	}

	double GetMeanEntropy()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetDirichletVal(i)->GetEntropy();
		}
		mean /= GetSize();
		return mean;
	}

	void SetUniform()	{
		for (int i=0; i<GetSize(); i++)	{
			GetDirichletVal(i)->setuniform();
		}
	}

	void SetSemiUniform()	{
		for (int i=0; i<GetSize(); i++)	{
			GetDirichletVal(i)->setsemiuniform();
		}
	}




	protected:

	Rvar<Profile>* CreateVal(int site)	{
		return new Dirichlet(center, concentration);
	}

	Var<Profile>* center;
	Var<PosReal>* concentration;
};

class DirichletIIDArrayMove : public MCUpdate	{

	public: 

	DirichletIIDArrayMove(DirichletIIDArray* inselectarray, double intuning, int inm) : selectarray(inselectarray), tuning(intuning), m(inm)	{}

	double Move(double tuning_modulator)	{
		double total = 0;

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif

			for (int i=0; i<selectarray->GetSize(); i++)	{
				total += selectarray->GetDirichletVal(i)->Move(tuning_modulator * tuning,m);
			}


		return total / selectarray->GetSize();
	}	

	private:

	DirichletIIDArray* selectarray;
	double tuning;
	int m;
};

class GammaIIDArray : public IIDArray<PosRealVector>	{
	
	public:

	GammaIIDArray(int insize, int indim, Var<PosReal>* inalpha, Var<PosReal>* inbeta ) : IIDArray<PosRealVector>(insize)	{

		dim = indim;
		alpha = inalpha;
		beta = inbeta;
		Create();
	}


	// make a MoveG(double tuning, int m)  calls over all sites
	virtual double MoveG(double tuning, int m)	{

		double tot = 0;

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetGammaVal(i)->Move(tuning,m);
		}
		return tot / GetSize();
	}

	IIDGamma* GetGammaVal(int site)	{
		return dynamic_cast<IIDGamma*>(GetVal(site));
	}

	double GetMeanVar()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetGammaVal(i)->GetVar();
		}
		mean /= GetSize();
		return mean;
	}

	protected:

	Rvar<PosRealVector>* CreateVal(int site)	{
		IIDGamma* tmp = new IIDGamma(dim,alpha,beta);
		return tmp;
	}

	int dim;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};


class GammaIIDArrayMove : public MCUpdate	{

	public:

	GammaIIDArrayMove(GammaIIDArray* ingamarray, double intuning, int inm)	{
		gamarray = ingamarray;
		tuning = intuning;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
	    return gamarray->MoveG(tuning_modulator * tuning, m);
	}

	private:

	GammaIIDArray* gamarray;
	double tuning;
	int m;
};


/*
class ExpIIDArray : public IIDArray<PosRealVector>	{
	
	public:

	ExpIIDArray(int insize, int indim, Var<PosReal>* inalpha ) : IIDArray<PosRealVector>(insize)	{

		dim = indim;
		alpha = inalpha;
		Create();
	}


	// make a MoveE(double tuning, int m)  calls over all sites
	virtual double MoveE(double tuning, int m)	{

		double tot = 0;

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetExpVal(i)->Move(tuning,m);
		}
		return tot / GetSize();
	}

	IIDExp* GetExpVal(int site)	{
		return dynamic_cast<IIDExp*>(GetVal(site));
	}

	double GetMeanVar()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetExpVal(i)->GetVar();
		}
		mean /= GetSize();
		return mean;
	}

	protected:

	Rvar<PosRealVector>* CreateVal(int site)	{
		IIDExp* tmp = new IIDExp(dim,alpha);
		return tmp;
	}

	int dim;
	Var<PosReal>* alpha;

};


class ExpIIDArrayMove : public MCUpdate	{

	public:

	ExpIIDArrayMove(ExpIIDArray* inexparray, double intuning, int inm)	{
		exparray = inexparray;
		tuning = intuning;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
	    return exparray->MoveE(tuning_modulator * tuning, m);
	}

	private:

	ExpIIDArray* exparray;
	double tuning;
	int m;
};

*/


class NormalIIDArray : public IIDArray<RealVector>	{

	public:

	NormalIIDArray(int insize,int indim, Var<Real>* inmean, Var<PosReal>*  invariance) : IIDArray<RealVector>(insize)	{

		mean = inmean;
		variance = invariance;
		dim = indim;
		Create();
	}

	void ClampAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetNormalVal(i)->ClampAtZero();
		}
	}

	// make a MoveN(double tuning, int m)  calls over all sites
	virtual double MoveN(double tuning, int m)	{

		double tot = 0;
//		int total, rank;

		#ifdef _OPENMP
//		cerr << "compiled by an openmp-compliant" << '\n';
//		cerr << omp_in_parallel() << '\n';;

//		cerr << "total = " << omp_get_num_threads() << '\n';;

		#pragma omp parallel for
		#endif

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetNormalVal(i)->Move(tuning,m);
		#ifdef _OPENMP
//		cerr << omp_in_parallel() << '\n';
//		cerr << "rank = " << omp_get_thread_num() << '\n';;

		#endif
		}

		return tot / GetSize();
	}




	IIDNormal* GetNormalVal(int site)	{
		return dynamic_cast<IIDNormal*>(GetVal(site));
	}




	double GetMeanVar()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetNormalVal(i)->GetVar();
		}
		mean /= GetSize();
		return mean;
	}


	protected:

	Rvar<RealVector>* CreateVal(int site)	{
		return new IIDNormal(dim, mean, variance);
	}

	Var<Real>* mean;
	Var<PosReal>* variance;
	int dim;
};


class NormalIIDArrayMove : public MCUpdate	{

	public:

	NormalIIDArrayMove(NormalIIDArray* inselectarray, double intuning, int inm)	{
		selectarray = inselectarray;
		tuning = intuning;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
	    return selectarray->MoveN(tuning_modulator * tuning, m);
	}

	private:

	NormalIIDArray* selectarray;
	double tuning;
	int m;
};




class MultinomialIIDArray : public IIDArray<IntVector>	{

	public:

	MultinomialIIDArray(int insize,int indim, Var<Profile>* inweight) : IIDArray<IntVector>(insize)	{

		dim = indim;
		weight = inweight;
		Create();
//		ComputeMeanCount()
	}


	// make a MoveM(double tuning, int m)  calls over all sites
	virtual double MoveM(int m)	{

		double tot = 0;
//		int total, rank;

		#ifdef _OPENMP
//		cerr << "compiled by an openmp-compliant" << '\n';
//		cerr << omp_in_parallel() << '\n';;

//		cerr << "total = " << omp_get_num_threads() << '\n';;

		#pragma omp parallel for
		#endif

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetMultinomialVal(i)->Move(m);
		#ifdef _OPENMP
//		cerr << omp_in_parallel() << '\n';
//		cerr << "rank = " << omp_get_thread_num() << '\n';

		#endif
		}

		return tot/GetSize();
	}



	IIDMultinomial* GetMultinomialVal(int site)	{
		return dynamic_cast<IIDMultinomial*>(GetVal(site));
	}


//counts the # of zero, poitive and negative mixture categories
	void ComputeMeanCount()	{

		totzero = 0;
		totpos = 0;
		totneg = 0;

		for (int i=0; i<GetSize(); i++)	{
			zero = 0;
			pos = 0;
			neg = 0;

			for (int j=0; j<dim; j++)	{

				if((*GetMultinomialVal(i))[j] == 0)	{  
					zero++;
				}
				if((*GetMultinomialVal(i))[j] == 1)	{
					pos++;
				}
				if((*GetMultinomialVal(i))[j] == 2)	{
					neg++;
				}
			}

			zero /= dim;
			pos /= dim;
			neg /= dim;

			totzero += zero;
			totpos += pos;
			totneg += neg;
		}

		totzero /= GetSize();
		totpos /= GetSize();
		totneg /= GetSize();

	}


	double GetZeroCount()	{
		return totzero;
	}


	double GetPosCount()	{
		return totpos;
	}


	double GetNegCount()	{
		return totneg;
	}




	protected:

	Rvar<IntVector>* CreateVal(int site)	{
		return new IIDMultinomial(weight, dim);
	}

	Var<Profile>* weight;
	int dim;
	double zero;
	double pos;
	double neg;


	public:

	double totzero;
	double totpos;
	double totneg;

};



class MultinomialIIDArrayMove : public MCUpdate	{

	public:

	MultinomialIIDArrayMove(MultinomialIIDArray* inmixcat, int inm)	{
		mixcat = inmixcat;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
	    return mixcat->MoveM(m);
	}

	private:

	MultinomialIIDArray* mixcat;
	int m;
};






class IIDSelector : public virtual Rvar<IntVector> {

	public:
	
	IIDSelector(Var<Profile>* inprobarray, int N)	{
		setval(IntVector(N));
		bkvalue = IntVector(N);
		probarray = inprobarray;
		Register(probarray);
		Sample();
	}

	virtual ~IIDSelector() {}
	
	void drawSample()	{
		for (int i=0; i<GetDim(); i++)	{
			int k = (int) (GetRange() * rnd::GetRandom().Uniform());
			(*this)[i] = k;	
		}
	}

	int GetRange()	{return probarray->GetDim();}

	virtual double	ProposeMove(double epsilon)	{
		for (int i=0; i<GetDim(); i++)	{
			if (rnd::GetRandom().Uniform() < epsilon)	{
				int k = (int) (GetRange() * rnd::GetRandom().Uniform());
				(*this)[i] = k;	
			}
		}
		return 0;
	}

	virtual double logProb()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += log((*probarray)[(*this)[i]]);
		}
		return total;
	}

	private:

	Var<Profile>* probarray;

};

class SelectorIIDArray : public IIDArray<IntVector>	{

	public:

	SelectorIIDArray(int insize,int indim, Var<Profile>* inweight) : IIDArray<IntVector>(insize)	{

		dim = indim;
		weight = inweight;
		Create();
	}

	IIDSelector* GetSelector(int site)	{
		return dynamic_cast<IIDSelector*>(GetVal(site));
	}

	//counts the # of zero, poitive and negative mixture categories
	double GetMeanCount(int k)	{

		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			for (int j=0; j<dim; j++)	{
				if((*GetSelector(i))[j] == k)	{  
					total++;
				}
			}
		}
		return total / GetSize() / dim;
	}

	virtual double MoveAll(double tuning)	{

		double tot = 0;

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetSelector(i)->Move(tuning);
		}
		return tot / GetSize();
	}

	protected:

	Rvar<IntVector>* CreateVal(int site)	{
		return new IIDSelector(weight, dim);
	}

	Var<Profile>* weight;
	int dim;

	public:


};


class SelectorIIDArrayMove : public MCUpdate	{

	public:

	SelectorIIDArrayMove(SelectorIIDArray* inmixcat, double intuning)	{
		mixcat = inmixcat;
		tuning = intuning;
	}

	double Move(double tuning_modulator = 1)	{
	    return mixcat->MoveAll(tuning);
	}

	private:

	SelectorIIDArray* mixcat;
	double tuning;
};





class GammaIIDMix: public Dvar<RealVector>  	{

	public:
	GammaIIDMix(int indim, IIDMultinomial* iniidmixcat, Var<PosRealVector>* iniidgammaPlus, Var<PosRealVector>* iniidgammaMinus) 	{
//	GammaIIDMix(int indim, IIDSelector* iniidmixcat, Var<PosRealVector>* iniidgammaPlus, Var<PosRealVector>* iniidgammaMinus) 	{

		setval(PosRealVector(iniidgammaPlus->GetDim()));
		bkvalue = PosRealVector(iniidgammaPlus->GetDim());
		dim = indim;
		iidmixcat = iniidmixcat;
		iidgammaPlus = iniidgammaPlus;
		iidgammaMinus = iniidgammaMinus;
		Register(iidmixcat);
		Register(iidgammaPlus);
		Register(iidgammaMinus);
		specialUpdate();

	}

	protected:
	
	void specialUpdate()	{

		for (int i=0; i<GetDim(); i++) {
			if ((*iidmixcat)[i] == 0)	{
				(*this)[i] = 0;
			}
			if ((*iidmixcat)[i] == 1)	{
				(*this)[i] = (*iidgammaPlus)[i];
			}
			if ((*iidmixcat)[i] == 2)	{
				(*this)[i] = -((*iidgammaMinus)[i]);
			}
		}
	}

	IIDMultinomial* iidmixcat;
//	IIDSelector* iidmixcat;

	Var<PosRealVector>* iidgammaPlus;
	Var<PosRealVector>* iidgammaMinus;
	int dim;

};


class GammaIIDMixArray: public ValPtrArray<Dvar<RealVector> > 	{

	public:
	GammaIIDMixArray(int insize, int indim, MultinomialIIDArray* inmixcat, GammaIIDArray* ingammaPlus, GammaIIDArray* ingammaMinus) : ValPtrArray<Dvar<RealVector> > (insize)	{
//	GammaIIDMixArray(int insize, int indim, SelectorIIDArray* inmixcat, GammaIIDArray* ingammaPlus, GammaIIDArray* ingammaMinus) : ValPtrArray<Dvar<RealVector> > (insize)	{


		dim = indim;
		mixcat = inmixcat;
		gammaPlus = ingammaPlus;
		gammaMinus = ingammaMinus;
		Create();
	}

	GammaIIDMix* GetGammaIIDMixVal(int site)	{
		return dynamic_cast<GammaIIDMix*>(GetVal(site));
	}

	double GetMeanVar()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetGammaIIDMixVal(i)->GetVar();
		}
		mean /= GetSize();
		return mean;
	}	

	protected:

	Dvar<RealVector>* CreateVal(int site)	{
		return new GammaIIDMix(dim, mixcat->GetMultinomialVal(site), gammaPlus->GetGammaVal(site), gammaMinus->GetGammaVal(site));
//		return new GammaIIDMix(dim, mixcat->GetSelector(site), gammaPlus->GetGammaVal(site), gammaMinus->GetGammaVal(site));
	}

	int dim;
	MultinomialIIDArray* mixcat;
//	SelectorIIDArray* mixcat;

	GammaIIDArray* gammaPlus;
	GammaIIDArray* gammaMinus;


};
