#ifndef JEFF_H
#define JEFF_H

#include "IID.h"

#include "ValArray.h"

class Jeffreys: public Rvar<PosReal>	{

	public:

	Jeffreys(double inmin, double inmax, DAGnode* inroot = 0)	{
		min = inmin;
		max = inmax;
		root = inroot;
		if (root)	{
			Register(root);
		}
		Sample();
	}

	double GetMin() {return min;}
	double GetMax() {return max;}

	double logProb()	{
		if (val() < min)	{
			return -1000;
		}
		if (val() > max)	{
			return -1000;
		}
		return -log(val());
	}

	void drawSample()	{
		setval(min * exp(Random::Uniform() * (log(max) - log(min))));
	}

	protected:

	double min;
	double max;
	DAGnode* root;
};

class Uniform: public Rvar<Real>	{

	public:

	Uniform(double inmin, double inmax, DAGnode* inroot = 0)	{
		min = inmin;
		max = inmax;
		root = inroot;
		if (root)	{
			Register(root);
		}
		Sample();
	}

	double logProb()	{
		if (val() < min)	{
			return -1000;
		}
		if (val() > max)	{
			return -1000;
		}
		return 0;
	}

	void drawSample()	{
		setval(min +Random::Uniform() * (max - min));
	}

	protected:

	double min;
	double max;
	DAGnode* root;
};

class JeffreysIIDArray : public IIDArray<PosReal>	{

	public:
	JeffreysIIDArray(int insize, double inmin, double inmax, DAGnode* inroot) : IIDArray<PosReal>(insize)	{
		min = inmin;
		max = inmax;
		root = inroot;
		Create();
	}

	Jeffreys* GetJeffreys(int i)	{
		Jeffreys* tmp = dynamic_cast<Jeffreys*>(GetVal(i));
		return tmp;
	}

	virtual void ClampAt(double in)	{
		for (int i=0; i<this->GetSize(); i++)	{
			this->GetVal(i)->ClampAt(in);
		}
	}

	virtual void ClampAt(double in, int i)	{
		this->GetVal(i)->ClampAt(in);
	}

	virtual void setval(double in)	{
		for (int i=0; i<this->GetSize(); i++)	{
			this->GetVal(i)->setval(in);
		}
	}
	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new Jeffreys(min,max,root);
	}

	double min;
	double max;
	DAGnode* root;
};

#endif
