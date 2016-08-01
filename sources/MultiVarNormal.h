
#ifndef MULTIVAR_H
#define MULTIVAR_H

#include "Normal.h"
#include "InverseWishartMatrix.h"

class MultiVarNormal : public IIDNormal	{


	public:

	MultiVarNormal(int dim, Var<Real>* inmean, Var<PosReal>* invar)	: IIDNormal(dim,inmean,invar) {
		Sample();
	}

	MultiVarNormal(Var<Real>* inmean, VarArray<PosReal>*  invararray) : IIDNormal(invararray->GetSize()) {
		// setval(RealVector(dim));
		// bkvalue = RealVector(dim);
		mean = inmean;
		meanvector = 0;
		variance = 0;
		vararray = invararray;
		Register(mean);
		for (int i=0; i<GetDim(); i++)	{
			Register(vararray->GetVal(i));
		}
		Sample();
	}

	~MultiVarNormal()	{
	}

	double logProb()	{

		if (variance)	{
			return IIDNormal::logProb();
		}
		double total = 0;
		if (mean)	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - mean->val();
				total -= 0.5 * log(2*Pi*vararray->GetVal(i)->val()) + 0.5 * tmp * tmp / vararray->GetVal(i)->val();
			}
		}
		else	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - (*meanvector)[i];
				total -= 0.5 * log(2*Pi*vararray->GetVal(i)->val()) + 0.5 * tmp * tmp / vararray->GetVal(i)->val();
			}
		}
		// cerr << "multivar normal log prob: " << total << '\n';
		return total;
	}


	protected:

	void drawSample()	{
		if (variance)	{
			IIDNormal::drawSample();
		}
		else	{
			if (mean)	{
				for (int i=0; i<GetDim(); i++)	{
					(*this)[i] = Random::sNormal() * sqrt(vararray->GetVal(i)->val()) + mean->val();
				}
			}
			else	{
				for (int i=0; i<GetDim(); i++)	{
					(*this)[i] = Random::sNormal() * sqrt(vararray->GetVal(i)->val()) + (*meanvector)[i];
				}
			}
		}
	}

	private:

	VarArray<PosReal>* vararray;

};

class MultiVarArray : public IIDArray<RealVector>	{

	public:

	MultiVarArray(int insize, Var<Real>* inmean, VarArray<PosReal>* invar, VarArray<PosReal>* invar0 = 0) : IIDArray<RealVector>(insize)	{
		mean = inmean;
		var = invar;
		if (invar0)	{
			var0 = invar0;
		}
		else	{
			var0 = var;
		}
		Create();
	}

	MultiVarNormal* GetMultiVarNormal(int i)	{
		return dynamic_cast<MultiVarNormal*>(GetVal(i));
	}

	void ClampAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetMultiVarNormal(i)->ClampAtZero();
		}
	}

	protected:

	Rvar<RealVector>* CreateVal(int mat)	{
		if (!mat)	{
			return new MultiVarNormal(mean,var);
		}
		return new MultiVarNormal(mean,var0);
	}

	private:

	Var<Real>* mean;
	VarArray<PosReal>* var;
	VarArray<PosReal>* var0;
};

class MultiUniArray : public IIDArray<RealVector>	{

	public:

	MultiUniArray(int insize, int indim, Var<Real>* inroot, double inmax) : IIDArray<RealVector>(insize)	{
		root = inroot;
		max = inmax;
		dim = indim;
		Create();
	}

	IIDUniform* GetIIDUniform(int i)	{
		return dynamic_cast<IIDUniform*>(GetVal(i));
	}

	void ClampAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetIIDUniform(i)->ClampAtZero();
		}
	}

	protected:

	Rvar<RealVector>* CreateVal(int mat)	{
		return new IIDUniform(root,dim,max);
	}

	private:

	double max;
	int dim;
	Var<Real>* root;
};

#endif

