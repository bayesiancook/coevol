#ifndef IIDGAMMA_H
#define IIDGAMMA_H

#include "IID.h"

class GammaArray : public IIDArray<PosReal>	{

	public:

	GammaArray(int insize, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDArray<PosReal>(insize)	{
		alpha = inalpha;
		beta = inbeta;
		Create();
	}


	double GetMean()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetVal(i)->val();
		}
		mean /= GetSize();
		return mean;
	}

	Gamma* GetGammaVal(int site)	{
		return dynamic_cast<Gamma*>(GetVal(site));
	}

/*
	virtual double MoveG(double tuning, int m)	{

		double tot = 0;

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetVal(i)->Move(tuning,m);
		}
		return tot / GetSize();
	}
*/

	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new Gamma(alpha, beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
};


class GammaArrayMove : public MCUpdate	{
  
	public:
	
	GammaArrayMove(GammaArray* inarray, double intuning) {
		array = inarray;
		tuning = intuning;
	}
	
	double Move(double tuning_modulator = 1)	{
	  
		double total = 0;
	  
		for (int i=0; i<array->GetSize(); i++)	{
			total += array->GetGammaVal(i)->Move(tuning_modulator * tuning);
		}
		return total/array->GetSize();
	}

	private:

	GammaArray* array;
	double tuning;
};  

class ExpArray : public IIDArray<PosReal>	{
	
	public:
	  
//	enum ParentType {MEAN,RATE};

	ExpArray(int insize, Var<PosReal>* inalpha, Exponential::ParentType intype) : IIDArray<PosReal>(insize)	{

		alpha = inalpha;
		type = intype;
		Create();
	}

	Exponential* GetExpVal(int site)	{
		return dynamic_cast<Exponential*>(GetVal(site));
	}

	
	double GetMean()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetVal(i)->val();
		}
		mean /= GetSize();
		return mean;
	}


/*
	virtual double MoveE(double tuning, int m)	{

		double tot = 0;

		for (int i=0; i<GetSize(); i++) {
			tot += this->GetVal(i)->Move(tuning,m);
		}
		return tot / GetSize();
	}
*/

	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new Exponential(alpha,type);
	}

	int dim;
	Var<PosReal>* alpha;
	Exponential::ParentType type;

};

class ExpArrayMove : public MCUpdate	{
  
	public:
	
	ExpArrayMove(ExpArray* inarray, double intuning) {
		array = inarray;
		tuning = intuning;
	}

	double Move(double tuning_modulator = 1)	{
	  
		double total = 0;
  
		for (int i=0; i<array->GetSize(); i++)	{
			total += array->GetExpVal(i)->Move(tuning_modulator * tuning);
		}
		return total/array->GetSize();
	}

	private:

	ExpArray* array;
	double tuning;
};  

#endif


