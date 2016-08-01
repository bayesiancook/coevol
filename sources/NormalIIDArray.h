
#ifndef NORMALIID_H
#define NORMALIID_H

#include "Normal.h"
#include "IID.h"
#include "Move.h"


class NormalIIDArray : public IIDArray<RealVector>	{

	public:
	NormalIIDArray(int indim, int insize, Var<Real>* inmean, Var<PosReal>* invar) : IIDArray<RealVector>(insize)	{
		dim = indim;
		mean = inmean;
		meanvector = 0;
		var = invar;
		Create();
	}

	NormalIIDArray(int insize, Var<RealVector>* inmeanvector, Var<PosReal>* invar) : IIDArray<RealVector>(insize)	{
		dim = inmeanvector->GetDim();
		meanvector = inmeanvector;
		mean = 0;
		var = invar;
		Create();
	}

	IIDNormal* operator[](int site)	{
		return dynamic_cast<IIDNormal*>(array[site]);
	}

	double GetGrandMean()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += (*this)[i]->GetMean();
		}
		mean /= GetSize();
		return mean;
	}

	double GetGrandVar()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetSize(); i++)	{
			double tmp = (*this)[i]->GetMean();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetSize();
		var /= GetSize();
		var -= mean * mean;
		return var;
	}

	double Move(double tuning, int n)	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += (*this)[i]->Move(tuning,n);
		}
		return total / this->GetSize();
	}

	void SetAtZero()	{
		for (int i=0; i<this->GetSize(); i++)	{
			(*this)[i]->SetAtZero();
		}
	}

	double PiecewiseTranslationMove(double tuning, int index, int n)	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += (*this)[i]->PiecewiseTranslationMove(tuning,index,n);
		}
		return total / this->GetSize();
	}

	protected:

	Rvar<RealVector>* CreateVal(int site)	{
		if (mean)	{
			return new IIDNormal(dim,mean,var);
		}
		return new IIDNormal(meanvector,var);
	}

	int dim;
	Var<Real>* mean;
	Var<RealVector>* meanvector;
	Var<PosReal>* var;
};


class NormalIIDArrayMove : public MCUpdate	{

	public:

	NormalIIDArrayMove(NormalIIDArray* invararray, double intuning, int inm) {
		vararray = invararray;
		tuning = intuning;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
		return vararray->Move(tuning * tuning_modulator,m);
	}

	private:

	NormalIIDArray* vararray;
	double tuning;
	int m;
};


class NormalIIDArrayPiecewiseTranslationMove : public MCUpdate	{

	public:

	NormalIIDArrayPiecewiseTranslationMove(NormalIIDArray* invararray, double intuning, int inindex, int inm)	{
		vararray = invararray;
		tuning = intuning;
		index = inindex;
		m = inm;
	}

	double Move(double tuning_modulator = 1)	{
		return vararray->PiecewiseTranslationMove(tuning * tuning_modulator,index,m);
	}

	private:

	NormalIIDArray* vararray;
	double tuning;
	int index;
	int m;
};

#endif

