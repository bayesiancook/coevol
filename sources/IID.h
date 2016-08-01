#ifndef IID_H
#define IID_H

#include <iostream>
#include <cstdlib>
using namespace std;

#include "RandomTypes.h"
#include "ValArray.h"

template <class V> class IIDArray : public ValPtrArray< Rvar<V> >	{

	public:

	IIDArray(int insize) : ValPtrArray<Rvar<V> >(insize) {}

	double Move(double tuning)	{
		double total = 0;
		int n = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			if (! this->GetVal(i)->isClamped())	{
				double tmp = this->GetVal(i)->Move(tuning);
				total += tmp;
				n++;
			}
		}
		double acc = n ? total / n : 0;
		return acc;
	}

	void drawSample()	{
		for (int i=0; i<this->GetSize(); i++)	{
			this->GetVal(i)->Sample();
		}
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += this->GetVal(i)->GetLogProb();
		}
		return total ;
	}

	void Register(DAGnode* in)	{
		for (int i=0; i<this->GetSize(); i++)	{
			this->GetVal(i)->Register(in);
		}
	}
};

class BetaIIDArray : public IIDArray<UnitReal>	{

	public:
	BetaIIDArray(int insize, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDArray<UnitReal>(insize)	{
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

	double GetVar()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetSize(); i++)	{
			double tmp = GetVal(i)->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetSize();
		var /= GetSize();
		var -= mean * mean;
		return var;
	}

	protected:

	Rvar<UnitReal>* CreateVal(int site)	{
		return new Beta(alpha, beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
};

class PosUniIIDArray : public IIDArray<PosReal>	{

	public:
	PosUniIIDArray(int insize, Var<PosReal>* inroot, double inmax) : IIDArray<PosReal>(insize)	{
		root = inroot;
		max = inmax;
		Create();
	}

	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new PosUniform(root, max);
	}

	Var<PosReal>* root;
	double max;
};

class GammaIIDArray : public IIDArray<PosReal>	{

	public:
	GammaIIDArray(int insize, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDArray<PosReal>(insize)	{
		alpha = inalpha;
		beta = inbeta;
		Create();
	}

	double* SetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->setval(*ptr++);
		}
		return ptr;
	}

	double* GetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			(*ptr++) = GetVal(i)->val();
		}
		return ptr;
	}

	Var<PosReal>* GetAlpha()	{
		return alpha;
	}

	Var<PosReal>* GetBeta()	{
		return beta;
	}

	double GetMean()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetVal(i)->val();
		}
		mean /= GetSize();
		return mean;
	}

	double GetVar()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetSize(); i++)	{
			double tmp = GetVal(i)->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetSize();
		var /= GetSize();
		var -= mean * mean;
		return var;
	}

	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new Gamma(alpha, beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
};

class DirichletIIDArray : public IIDArray<Profile>	{

	public:
	DirichletIIDArray(int insize, Var<Profile>* incenter, Var<PosReal>* inconcentration) : IIDArray<Profile>(insize)	{
		center = incenter;
		concentration = inconcentration;
		Create();
	}

	double* SetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			for (int k=0; k<GetDim(); k++)	{
				(*GetVal(i))[k] = (*ptr++);
			}
		}
		return ptr;
	}

	double* GetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			for (int k=0; k<GetDim(); k++)	{
				(*ptr++) = (*GetVal(i))[k];
			}
		}
		return ptr;
	}

	int GetDim()	{
		return GetVal(0)->GetDim();
	}

	protected:

	Rvar<Profile>* CreateVal(int site)	{
		return new Dirichlet(center,concentration);
	}

	Var<Profile>* center;
	Var<PosReal>* concentration;
};

class NormalIIDArray : public IIDArray<Real>	{

	public:
	NormalIIDArray(int insize, Var<Real>* inmean, Var<PosReal>* invar) : IIDArray<Real>(insize)	{
		mean = inmean;
		var = invar;
		Create();
	}

	double* SetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->setval(*ptr++);
		}
		return ptr;
	}

	double* GetVals(double* ptr)	{
		for (int i=0; i<GetSize(); i++)	{
			(*ptr++) = GetVal(i)->val();
		}
		return ptr;
	}

	void Reset()	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->setval(0);
		}
	}

	Normal* GetNormal(int site)	{
		Normal* n = dynamic_cast<Normal*>(GetVal(site));
		if (! n)	{
			cerr << "error in NormalIIDArray::GetNormal\n";
			exit(1);
		}
		return n;
	}

	protected:

	Rvar<Real>* CreateVal(int site)	{
		return new Normal(mean,var);
	}

	Var<Real>* mean;
	Var<PosReal>* var;
};


#endif // IID_H
