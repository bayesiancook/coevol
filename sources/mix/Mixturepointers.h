#ifndef MIXTURE_H
#define MIXTURE_H

#include <iostream>
#include <cstdlib>
using namespace std;

#include "RandomTypes.h"
#include "ValArray.h"

template <class V> class PointerKeeper : public Dvar<V>	{

	public:

	PointerKeeper(Rnode* incomponent = 0) : component(incomponent) {
		this->Register(component);
	}

	void DetachFrom()	{
		if (component)	{
			for (clit i=this->down.begin(); i!=this->down.end(); i++)	{ 
				(*i)->DeregisterFrom(component);
			}
			this->DeregisterFrom(component);
		}
	}

	void AttachTo(Rnode* newcomponent)	{
		component = newcomponent;
		this->Register(component);
		for (clit i=this->down.begin(); i!=this->down.end(); i++)	{ 
			(*i)->Register(component);
		}
	}

	void specialUpdate() {}

	private:

	Rnode* component;
	
};

/*
template <class T> class MixtureAware : public T	{

	MixtureAware(const T& from) : T(from) , p(0) {

	}

	void SetCurrentPointerKeeper(PointerKeeper* inp) {p = inp;}

	RegisterChild(DAGnode* child)	{
		child->Register(p);
		child->Register(this);
	}
		
	private:

	PointerKeeper* p;

};
*/

template <class U, class V = U> class FiniteMixture : public RandomVarArray<V> {


	public:

	FiniteMixture(int insize, int incomponentnumber)	{
		size = insize;
		Ncomponent = incomponentnumber;
	}

	int GetComponentNumber()	{
		return Ncomponent;
	}

	int GetSize()	{
		return size;
	}

	Var<V>* GetVal(int site)	{
		return GetComponent(GetAllocation(site)); 
	}
	
	Var<V>* operator[](int site)	{
		return GetComponent(GetAllocation(site)); 
	}

	
	void Create()	{
		weight = new Dirichlet(Ncomponent);
		allocation = new FiniteDiscrete*[size];
		for (int i=0; i<this->GetSize(); i++)	{
			allocation[i] = new FiniteDiscrete(weight);
		}
		component = new Rvar<U>*[GetComponentNumber()];
		for (int k=0; k<GetComponentNumber(); k++)	{
			component[k] = CreateComponent(k);
		}
		grid = new PointerKeeper<V>*[this->GetSize()];
		for (int i=0; i<this->GetSize(); i++)	{
			grid[i] = new PointerKeeper<V>();
		}
		// this->Sample();
	}

	int GetAllocation(int i)	{
		return allocation[i]->val();
	}

	double GetAllocationLogProb(int i)	{
		return allocation[i]->logProb();
	}
	

	double GetWeight(int k)	{
		return (*weight)[k];
	}

	Var<Profile>* GetWeightVector()	{
		return weight;
	}

	double GetWeightLogProb(int k)	{
		return weight->GetLogProb();
	}

	Rvar<U>* GetComponent(int k)	{
		if ((k<0) || (k>=GetComponentNumber()))	{
			cerr << "error in FiniteMixture::SetAllocation : component " << k << " out of " << GetComponentNumber() << '\n';
			throw;
		}
		return component[k];
	}

	void SampleWeights()	{
		weight->Sample();
	}

	void SampleAllocation()	{
		for (int i=0; i<this->GetSize(); i++)	{
			allocation[i]->Sample();
			SetAllocation(i,allocation[i]->val());
		}
	}

	void SampleValues()	{
		for (int k=0; k<GetComponentNumber(); k++)	{
			GetComponent(k)->Sample();
		}
	}

	double MoveAllocation(int nrep = 1)	{
		double* logp = new double[GetComponentNumber()];
		double* p = new double[GetComponentNumber()];

		int nsuccess = 0;

		for (int rep=0; rep<nrep; rep++)	{
			for (int i=0; i<this->GetSize(); i++)	{

				int bk = GetAllocation(i);

				// compute probs of each allocation
				double max = 0;
				/*
				double prev = 0;
				for (int k=0; k<GetComponentNumber(); k++)	{
					logp[k] = prev + SetAllocation(i,k);
					prev = logp[k];
					if ((!k) || (max < logp[k]))	{
						max = logp[k];
					}
				}
				*/
				for (int k=0; k<GetComponentNumber(); k++)	{
					SetAllocation(i,bk);
					logp[k] = SetAllocation(i,k);
					if ((!k) || (max < logp[k]))	{
						max = logp[k];
					}
				}
				for (int k=0; k<GetComponentNumber(); k++)	{
					p[k] = exp(logp[k] - max);
				}

				// gibbs choice
				int choose = Random::FiniteDiscrete(GetComponentNumber(),p);
				if (choose == bk)	{
					nsuccess++;
				}
				SetAllocation(i,choose);
			}
		}

		delete[] p;
		delete[] logp;
		return ((double) nsuccess) / this->GetSize() / nrep;
	}

	double MoveValues(double tuning)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetComponent(k)->Move(tuning);
		}
		return total / GetComponentNumber();
	};

	double Move(double tuning)	{
		double tot = 0;
		MoveAllocation();
		MoveValues(tuning);
		tot += weight->Move(tuning);
		tot += weight->Move(tuning/10);
		tot += weight->Move(tuning/100);
		return tot / 3;
	}

	double GetLogProb()	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetComponent(k)->GetLogProb();
		}
		for (int i=0; i<this->GetSize(); i++)	{
			total += allocation[i]->GetLogProb();
		}
		total += weight->GetLogProb();
		return total;
	}

	void ToStream(ostream& os)	{
		os << *weight << '\n';
		for (int k=0; k<GetComponentNumber(); k++)	{
			os << *(component[k]) << '\t';
		}
		os << '\n';
		for (int i=0; i<GetSize(); i++)	{
			os << *allocation[i] << '\t';
		}
	}

	void FromStream(istream& is)	{
		is >> *weight;
		if (weight->GetDim() != GetComponentNumber())	{
			cerr << "error in FiniteMixture::FromStream: wrong number of components : " << weight->GetDim() << " instead of " << GetComponentNumber() << '\n';
			throw;
		}
		for (int k=0; k<GetComponentNumber(); k++)	{
			is >> *(component[k]);
		}
		for (int i=0; i<GetSize(); i++)	{
			is >> *allocation[i];
		}
	}

	friend ostream& operator<<(ostream& os, FiniteMixture& mix)	{
		mix.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, FiniteMixture& mix)	{
		mix.FromStream(is);
		return is;
	}
		
	protected:

	void Corrupt(int i, bool bk=false)	{
		grid[i]->Corrupt(bk);
		allocation[i]->Corrupt(bk);
	}

	double Update(int i)	{
		return grid[i]->Update() + allocation[i]->Update();
	}

	double SetAllocation(int i, int k)	{
		Corrupt(i);
		grid[i]->DetachFrom();
		allocation[i]->setval(k);
		grid[i]->AttachTo(GetComponent(k));
		return Update(i);
	}

	virtual Rvar<U>* CreateComponent(int k) = 0;

	void drawSample()	{
		SampleWeights();
		SampleAllocation();
		SampleValues();
	}

	int size;
	int Ncomponent;
	FiniteDiscrete** allocation;
	Rvar<U>** component;
	PointerKeeper<V>** grid;
	Rvar<Profile>* weight;

};

class RateFiniteMixture : public FiniteMixture<PosReal>	{

	public:

	RateFiniteMixture(int insize, int incomponentnumber) : FiniteMixture<PosReal>(insize, incomponentnumber) {}

	void Register(DAGnode* in)	{
		for (int k=0; k<GetComponentNumber(); k++)	{
			GetComponent(k)->Register(in);
		}
	}

	int ScalarMultiplication(double d)	{
		for (int k=0; k<GetComponentNumber(); k++)	{
			GetComponent(k)->ScalarMultiplication(d);
		}
		return GetComponentNumber();
	}

};

class GammaMixture : public RateFiniteMixture	{

	public:

	GammaMixture(int insize, int incomponentnumber, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : RateFiniteMixture(insize, incomponentnumber)	{
		alpha = inalpha;
		beta = inbeta;
		Create();
	}

	protected:

	Rvar<PosReal>* CreateComponent(int k)	{
		cerr << *alpha << '\t' << *beta << '\n';
		return new Gamma(alpha,beta);
	}

	private:
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};

#endif
