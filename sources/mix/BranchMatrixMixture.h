#ifndef BRANCHMATRIXMIXTURE_H
#define BRANCHMATRIXMIXTURE_H

#include "MatrixMixture.h"
#include "RandomSubMatrix.h"

// should maintain a pointer to a BranchValTree<RandomSubMatrix>*
//
class BranchPointerKeeper : public BranchValPtrTree<PointerKeeper>	{

	public:

	BranchPointerKeeper(RootPointer* inroot, BranchValPtrTree<RandomSubMatrix>* incomponent) : component(incomponent) {
		root = inroot;
		SetWithRoot(component->WithRoot());
		RecursiveCreate(GetRoot());
		AttachTo(component);
	}

	Tree* GetTree() {return component->GetTree();}

	void DetachFrom()	{
		if (component)	{
			RecursiveDetachFrom(GetRoot());
		}
	}

	void AttachTo(BranchValPtrTree<RandomSubMatrix>* newcomponent)	{
		component = newcomponent;
		RecursiveAttachTo(GetRoot(),newcomponent);
	}

	void Corrupt(bool bk=false)	{
		RecursiveCorrupt(GetRoot(),bk);
	}

	double Update()	{
		return RecursiveUpdate(GetRoot());
	}

	private:

	void RecursiveDetachFrom(const Link* from)	{
		if (WithRoot())	{
			GetBranchVal(from->GetBranch())->DetachFrom();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDetachFrom(link->Out());
			GetBranchVal(link->GetBranch())->DetachFrom();
		}
	}

	void RecursiveAttachTo(const Link* from, BranchValPtrTree<RandomSubMatrix>* matrixtree)	{
		if (WithRoot())	{
			GetBranchVal(from->GetBranch())->AttachTo(matrixtree->GetBranchVal(from->GetBranch()));
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAttachTo(link->Out(), matrixtree);
			GetBranchVal(link->GetBranch())->AttachTo(component->GetBranchVal(link->GetBranch()));
		}
	}

	void RecursiveCorrupt(const Link* from, bool bk)	{
		if (WithRoot())	{
			GetBranchVal(from->GetBranch())->Corrupt(bk);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCorrupt(link->Out(),bk);
			GetBranchVal(link->GetBranch())->Corrupt(bk);
		}
	}

	double RecursiveUpdate(const Link* from)	{
		double total = 0;
		if (WithRoot())	{
			total += GetBranchVal(from->GetBranch())->Update();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveUpdate(link->Out());
			total += GetBranchVal(link->GetBranch())->Update();
		}
		return total;
	}

	PointerKeeper* CreateBranchVal(const Link* link)	{
		return new PointerKeeper(root);
	}

	BranchValPtrTree<RandomSubMatrix>* component;
	RootPointer* root;
};


template <class U> class MixtureRandomMatrixTree : public BranchValPtrTree<RandomSubMatrix> {

	public:

	MixtureRandomMatrixTree()	{
		SetWithRoot(true);
	}
	Rvar<U>* GetRandomVariable() {return rvar;}
	RandomSubMatrix* GetRandomSubMatrix(const Branch* branch) {return GetBranchVal(branch);}

	protected:

	virtual Rvar<U>* CreateRandomVariable() = 0;
	// virtual RandomSubMatrix* CreateRandomSubMatrix(const Link* link) = 0;

	Rvar<U>* rvar;
};


template <class U> class BranchMatrixFiniteMixture  : public MCMC {

	public:

	BranchMatrixFiniteMixture(Tree* intree, int insize, int incomponentnumber)	{
		tree = intree;
		size = insize;
		Ncomponent = incomponentnumber;
	}

	DAGnode* GetRoot() {return root;}

	void Register(DAGnode* innode, int site, const Branch* branch)	{
		innode->Register(grid[site]->GetBranchVal(branch));
	}

	int GetComponentNumber()	{
		return Ncomponent;
	}

	int GetSize()	{
		return size;
	}

	RandomSubMatrix* GetMatrix(int site, const Branch* branch)	{
		return component[GetAllocation(site)]->GetRandomSubMatrix(branch); 
	}
	
	RandomSubMatrix* GetComponentMatrix(int k, const Branch* branch) {
		return component[k]->GetRandomSubMatrix(branch); 
	}
	
	BranchValPtrTree<RandomSubMatrix>* GetComponentMatrixTree(int k) {
		return component[k]; 
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
		component = new MixtureRandomMatrixTree<U>*[GetComponentNumber()];
		for (int k=0; k<GetComponentNumber(); k++)	{
			component[k] = CreateComponent(k);
		}
		grid = new BranchPointerKeeper*[this->GetSize()];
		for (int i=0; i<this->GetSize(); i++)	{
			grid[i] = new BranchPointerKeeper(root,GetComponentMatrixTree(0));
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
		// ??
	}

	double Move(double tuning)	{
		double tot = 0;
		tot += MoveAllocation();
		tot += MoveValues(tuning);
		ResampleWeights();
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

	friend ostream& operator<<(ostream& os, BranchMatrixFiniteMixture& mix)	{
		mix.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, BranchMatrixFiniteMixture& mix)	{
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
		grid[i]->AttachTo(GetComponentMatrixTree(k));
		return Update(i);
	}

	virtual MixtureRandomMatrixTree<U>* CreateComponent(int k) = 0;

	void drawSample()	{
		SampleWeights();
		SampleAllocation();
		SampleValues();
	}

	RootPointer* root;
	Tree* tree;
	int size;
	int Ncomponent;
	FiniteDiscrete** allocation;
	MixtureRandomMatrixTree<U>** component;
	BranchPointerKeeper** grid;
	Rvar<Profile>* weight;

};

template<class U> class BranchMatrixMixturePhyloProcess : public PhyloProcess	{
	
	protected:

	public:

	BranchMatrixMixturePhyloProcess(LengthTree* intree, BranchMatrixFiniteMixture<U>* inmatmix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matmix = inmatmix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	// 
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		RandomBranchSitePath* path = new  RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,matmix->GetMatrix(site,link->GetBranch()),0);
		matmix->Register(path,site,link->GetBranch());
		return path;
	}

	protected:
	BranchMatrixFiniteMixture<U>* matmix;
};

template<class U> class BranchMatMixValMove : public MCUpdate	{

	public:

	BranchMatMixValMove(BranchMatrixFiniteMixture<U>* inmix, double intuning, int innrep) : mix(inmix), tuning(intuning), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator);
		}
		return total / nrep;
	}	
	
	protected:

	BranchMatrixFiniteMixture<U>* mix;
	double tuning;
	int nrep;
};



template<class U> class BranchMatMixWeightAllocMove : public MCUpdate	{

	public:

	BranchMatMixWeightAllocMove(BranchMatrixFiniteMixture<U>* inmix, int innrep) : mix(inmix), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double ret = mix->MoveAllocation(nrep);
		mix->ResampleWeights();
		return ret;
	}	
	
	protected:

	BranchMatrixFiniteMixture<U>* mix;
	int nrep;
};


#endif