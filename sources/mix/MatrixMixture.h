#ifndef MATRIXMIXTURE_H
#define MATRIXMIXTURE_H

#include "Mixture.h"
#include "RandomSubMatrix.h"
#include "RandomBranchSitePath.h"
#include "PhyloProcess.h"
#include "Move.h"


class RootPointer : public Dnode	{
	
	public:

	void specialUpdate() {}

};

class PointerKeeper : public Dnode	{

	public:

	PointerKeeper(RootPointer* inroot) : component(0) , root(inroot) {
		Register(root);
		SetName("pointer keeper");
	}

	void DetachFrom()	{
		if (component)	{
			for (clit i=this->down.begin(); i!=this->down.end(); i++)	{ 
				(*i)->DeregisterFrom(component);
			}
		}
	}

	void AttachTo(RandomSubMatrix* newcomponent)	{
		component = newcomponent;
		for (clit i=this->down.begin(); i!=this->down.end(); i++)	{ 
			(*i)->Register(component);
			RandomBranchSitePath* path = dynamic_cast<RandomBranchSitePath*>(*i);
			if (path)	{
				path->SetMatrix(component);
			}
		}
	}

	void specialUpdate() {}

	private:

	RandomSubMatrix* component;
	RootPointer* root;
};

template <class U> class MixtureRandomMatrix {

	public:

	MixtureRandomMatrix() {}
	virtual ~MixtureRandomMatrix() {}

	virtual Rvar<U>* CreateRandomVariable() = 0;
	virtual RandomSubMatrix* CreateRandomSubMatrix() = 0;

	Rvar<U>* GetRandomVariable() {return rvar;}
	RandomSubMatrix* GetRandomSubMatrix() {return matrix;}

	protected:
	Rvar<U>* rvar;
	RandomSubMatrix* matrix;
};


template <class U> class MatrixFiniteMixture  : public MCMC {

	public:

	MatrixFiniteMixture(int insize, int incomponentnumber)	{
		size = insize;
		Ncomponent = incomponentnumber;
	}

	void Register(DAGnode* innode, int site)	{
		innode->Register(grid[site]);
	}

	RootPointer* GetRoot() {return root;}

	int GetComponentNumber()	{
		return Ncomponent;
	}

	int GetSize()	{
		return size;
	}

	RandomSubMatrix* GetMatrix(int site)	{
		return component[GetAllocation(site)]->GetRandomSubMatrix(); 
	}
	
	RandomSubMatrix* GetComponentMatrix(int k) {
		return component[k]->GetRandomSubMatrix(); 
	}
	
	Rvar<U>* GetRandomVariable(int site)	{
		return GetComponent(GetAllocation(site)); 
	}

	Rvar<U>* operator[](int site)	{
		return GetComponent(GetAllocation(site)); 
	}

	Rvar<U>* GetComponent(int k)	{
		if ((k<0) || (k>=GetComponentNumber()))	{
			cerr << "error in FiniteMixture::SetAllocation : component " << k << " out of " << GetComponentNumber() << '\n';
			throw;
		}
		return component[k]->GetRandomVariable();
	}

	void Create()	{
		root = new RootPointer();
		weight = new Dirichlet(Ncomponent);
		allocation = new FiniteDiscrete*[size];
		for (int i=0; i<this->GetSize(); i++)	{
			allocation[i] = new FiniteDiscrete(weight);
			allocation[i]->setval(0);
		}
		component = new MixtureRandomMatrix<U>*[GetComponentNumber()];
		for (int k=0; k<GetComponentNumber(); k++)	{
			component[k] = CreateComponent(k);
		}
		grid = new PointerKeeper*[this->GetSize()];
		for (int i=0; i<this->GetSize(); i++)	{
			grid[i] = new PointerKeeper(root);
			grid[i]->AttachTo(GetComponentMatrix(0));
		}
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

	void SampleWeights()	{
		weight->Sample();
	}

	void ResampleWeights()	{
		weight->Corrupt(false);
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			double tmp = Random::sGamma(1 + GetAllocationStatistic(k));
			total += tmp;
			(*weight)[k] = tmp;
		}
		for (int k=0; k<GetComponentNumber(); k++)	{
			(*weight)[k] /= total;
		}
		weight->Update();
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
				for (int k=0; k<GetComponentNumber(); k++)	{
					logp[k] = SetAllocation(i,k);
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
				if (choose != bk)	{
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

	int GetAllocationStatistic(int k)	{
		int tot = 0;
		for (int i=0; i<GetSize(); i++)	{
			if (GetAllocation(i) == k)	{
				tot++;
			}
		}
		return tot;
	}

	void Check()	{
		for (int i=0; i<GetSize(); i++)	{
			if (! grid[i]->flag)	{
				cerr << "grid : " << i << '\n';
				exit(1);
			}
		}
		for (int k=0; k<GetComponentNumber(); k++)	{
			cerr << GetAllocationStatistic(k) << '\t';
		}
		cerr << '\n';
	}

	double Move(double tuning)	{
		double tot = 0;
		tot += MoveAllocation();
		tot += MoveValues(tuning);
		ResampleWeights();
		/*
		tot += weight->Move(tuning);
		tot += weight->Move(tuning/10);
		tot += weight->Move(tuning/100);
		*/
		return tot / 2;
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

	double GetWeightEntropy()	{
		return weight->GetEntropy();
	}

	double GetEffSize()	{
		return exp(weight->GetEntropy());
	}

	void ToStream(ostream& os)	{
		os << *weight << '\n';
		cout << "::" << *weight << "::" << '\n';
		for (int k=0; k<GetComponentNumber(); k++)	{
			os << *(GetComponent(k)) << '\t';
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
			is >> *(GetComponent(k));
		}
		for (int i=0; i<GetSize(); i++)	{
			is >> *allocation[i];
		}
	}

	friend ostream& operator<<(ostream& os, MatrixFiniteMixture& mix)	{
		mix.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, MatrixFiniteMixture& mix)	{
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
		grid[i]->AttachTo(GetComponentMatrix(k));
		return Update(i);
	}

	virtual MixtureRandomMatrix<U>* CreateComponent(int k) = 0;

	void drawSample()	{
		SampleWeights();
		SampleAllocation();
		SampleValues();
	}

	RootPointer* root;
	int size;
	int Ncomponent;
	FiniteDiscrete** allocation;
	MixtureRandomMatrix<U>** component;
	PointerKeeper** grid;
	Rvar<Profile>* weight;

};

template<class U> class MatMixValMove : public MCUpdate	{

	public:

	MatMixValMove(MatrixFiniteMixture<U>* inmix, double intuning, int innrep) : mix(inmix), tuning(intuning), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator);
		}
		return total / nrep;
	}	
	
	protected:

	MatrixFiniteMixture<U>* mix;
	double tuning;
	int nrep;
};



template<class U> class MatMixWeightAllocMove : public MCUpdate	{

	public:

	MatMixWeightAllocMove(MatrixFiniteMixture<U>* inmix, int innrep) : mix(inmix), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double ret = mix->MoveAllocation(nrep);
		mix->ResampleWeights();
		return ret;
	}	
	
	protected:

	MatrixFiniteMixture<U>* mix;
	int nrep;
};


template<class U> class MatrixMixturePhyloProcess : public PhyloProcess	{
	
	protected:

	public:

	MatrixMixturePhyloProcess(LengthTree* intree, MatrixFiniteMixture<U>* inmatmix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matmix = inmatmix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	// 
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		RandomBranchSitePath* path = new  RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(site),0);
		matmix->Register(path,site);
		return path;
	}

	RandomSubMatrix* 	GetMatrix(int site) {return matmix->GetMatrix(site);}

	protected:
	MatrixFiniteMixture<U>* matmix;
};

#endif
