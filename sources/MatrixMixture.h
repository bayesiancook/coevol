#ifndef MATRIXMIXTURE_H
#define MATRIXMIXTURE_H

#include "Mixture.h"
#include "RandomSubMatrix.h"
#include "RandomBranchSitePath.h"
#include "PhyloProcess.h"
#include "Move.h"

class PointerKeeper : public Dnode	{

	public:

	PointerKeeper() : component(0) {
		SetName("pointer keeper");
		flag=true;
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
};


template <class U> class MixtureRandomMatrix {

	public:

	MixtureRandomMatrix() {}
	virtual ~MixtureRandomMatrix() {
		delete matrix;
		delete rvar;
	}

	virtual Rvar<U>* CreateRandomVariable() = 0;
	virtual RandomSubMatrix* CreateRandomSubMatrix() = 0;

	Rvar<U>* GetRandomVariable() {return rvar;}
	RandomSubMatrix* GetRandomSubMatrix() {return matrix;}

	protected:
	Rvar<U>* rvar;
	RandomSubMatrix* matrix;
	map<DAGnode*,int> nodemap;
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
			grid[i] = new PointerKeeper();
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
	}

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
			int temp;
			is >> temp;
			SetAllocation(i, temp);
			//is >> *allocation[i];
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

	// temporarily 'unprotected'
	double SetAllocation(int i, int k)	{
		//Corrupt(i);
		Corrupt(i,true);
		grid[i]->DetachFrom();
		allocation[i]->setval(k);
		grid[i]->AttachTo(GetComponentMatrix(k));
		return Update(i);
	}

	protected:

	void Corrupt(int i, bool bk=false)	{
		grid[i]->Corrupt(bk);
		allocation[i]->Corrupt(bk);
	}

	double Update(int i)	{
		return grid[i]->Update() + allocation[i]->Update();
	}

	virtual MixtureRandomMatrix<U>* CreateComponent(int k) = 0;

	void drawSample()	{
		SampleWeights();
		SampleAllocation();
		SampleValues();
	}

	int size;
	int Ncomponent;
	FiniteDiscrete** allocation;
	MixtureRandomMatrix<U>** component;
	PointerKeeper** grid;
	Rvar<Profile>* weight;

};

//-------------------------------------------
//	Nic Dirichlet process
//-------------------------------------------


template <class U> class MatrixInfiniteMixture  : public MCMC {

	public:

	MatrixInfiniteMixture(int insize, int incomponentnumber)	{
		size = insize;
		Ncomponent = incomponentnumber;
		NmaxComponents = size + 5; // arbitrarily set to five more than size, for now...
		alpha = 10;
		alphaMin = 0.01;
		alphaMax = 1000;
	}



	void Register(DAGnode* innode, int site)	{
		innode->Register(grid[site]);
	}

	int GetComponentNumber()	{
		return Ncomponent;
	}

	void IncreaseComponentNumber()	{
		Ncomponent++;
	}

	void DecreaseComponentNumber()	{
		Ncomponent--;
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
		allocation = new int[size];
		for (int i=0; i<this->GetSize(); i++)	{
			allocation[i] = 0;
		}
		//component = new MixtureRandomMatrix<U>*[GetComponentNumber()];
		component = new MixtureRandomMatrix<U>*[NmaxComponents];
		for (int k=0; k<GetComponentNumber(); k++)	{
			component[k] = CreateComponent(k);
		}
		grid = new PointerKeeper*[this->GetSize()];
		for (int i=0; i<this->GetSize(); i++)	{
			grid[i] = new PointerKeeper();
			grid[i]->AttachTo(GetComponentMatrix(0));
			// Corrupt(i,true);
			// Update(i);
		}
	}

	int GetAllocation(int i)	{
		//return allocation[i]->val();
		return allocation[i];
	}

	//double GetAllocationLogProb(int i)	{
	//	return allocation[i]->logProb();
	//}


	//double GetWeight(int k)	{
	//	return (*weight)[k];
	//}

	//Var<Profile>* GetWeightVector()	{
	//	return weight;
	//}

	//double GetWeightLogProb(int k)	{
	//	return weight->GetLogProb();
	//}

	//void SampleWeights()	{
	//	weight->Sample();
	//}

	/*
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
	*/

	void SampleValues()	{
		for (int k=0; k<GetComponentNumber(); k++)	{
			GetComponent(k)->Sample();
		}
	}

	double MoveAllocation(int nrep = 1, int nauxiliary = 5)	{

		//double* logp = new double[GetComponentNumber()];
		//double* p = new double[GetComponentNumber()];
		double* logp = new double[GetSize() + nauxiliary];
		double* p = new double[GetSize() + nauxiliary];
		double templogp;
		int tempallocstat;
		int h;
		int bkallocstat;
		int nsuccess = 0;

		for (int rep=0; rep<nrep; rep++)	{
			for (int i=0; i<this->GetSize(); i++)	{

				// check not beyond NmaxComponents
				if ((GetComponentNumber() + nauxiliary) > NmaxComponents)	{
					cerr << "beyond NmaxComponent...\n";
					exit(1);
				}

				int bk = GetAllocation(i);
				bkallocstat = GetAllocationStatistic(bk);
				if (bkallocstat > 1)	{
					h = GetComponentNumber() + nauxiliary;
				}
				else	{
					h = GetComponentNumber() + nauxiliary - 1;
				}

				// h = GetComponentNumber();

				for (int j=GetComponentNumber(); j<h; j++)	{
					component[j] = CreateComponent(j);
				}

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
				//for (int k=0; k<GetComponentNumber(); k++)	{
				for (int k=0; k<h; k++)	{
					SetAllocation(i,bk);
					templogp = SetAllocation(i,k);
					tempallocstat = GetAllocationStatistic(k)-1;
					if (tempallocstat) {
						logp[k] = log((double)(tempallocstat)) + templogp;
					}
					else {
						logp[k] = log(alpha/nauxiliary) + templogp;
					}

					if ((!k) || (max < logp[k]))	{
						max = logp[k];
					}
				}
				for (int k=0; k<h; k++)	{
					p[k] = exp(logp[k] - max);
				}

				// gibbs choice
				int choose = Random::FiniteDiscrete(h,p);
				if (choose != bk)	{
					nsuccess++;
				}

				SetAllocation(i,choose);

				// new component
				if (choose >= GetComponentNumber())	{
					if (choose > GetComponentNumber())	{
						SwapComponents(choose, GetComponentNumber());
					}
					IncreaseComponentNumber();
					if (GetComponentNumber() == this->GetSize() + 1)	{
						cerr << "more components than observations...\n";
						exit(1);
					}
				}

				// delete unchosen auxiliary components
				for (int k=GetComponentNumber(); k<h; k++)	{
					delete component[k];
				}

				// dropped component
				if ((choose != bk) && (bkallocstat == 1))	{
					SwapComponents(bk, GetComponentNumber()-1);
					delete component[GetComponentNumber()-1];
					DecreaseComponentNumber();
				}

				/*
				cerr << "swap 0 1\n";
				int n1 = (int) (GetComponentNumber() * Random::Uniform());
				int n2 = (int) ((GetComponentNumber() -1) * Random::Uniform());
				if (n2 >= n1)	{
					n2++;
				}
				SwapComponents(n1,n2);
				*/
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
	}


	double MoveAlpha(double tuning, int nrep)	{
		int naccepted = 0;
		double bkalpha;
		double h, e;
		double deltaLogPrior;
		double deltaLog;
		for (int rep=0; rep<nrep; rep++)	{
			bkalpha = alpha;
			deltaLogPrior = -GetAlphaLogPrior();
			h = (Random::Uniform() - 0.5) * tuning;
			e = exp(h);
			alpha *= e;
			if (alpha < alphaMin) {
				alpha = 2 * alphaMin - alpha;
			}
			if (alpha > alphaMax)   {
				alpha = 2 * alphaMax - alpha;
			}
			deltaLogPrior += GetAlphaLogPrior();
			deltaLog = GetComponentNumber() * log(bkalpha);
			deltaLog -= GetComponentNumber() * log(alpha);
			for (int i=0; i<GetSize(); i++)	{
				deltaLog -= log(bkalpha + i);
				deltaLog += log(alpha + i);
			}
			deltaLog += deltaLogPrior;
			int accepted = (-log(Random::Uniform()) > deltaLog);
			if (accepted)	{
				naccepted++;
			}
			else	{
				alpha = bkalpha;
			}
		}
		return (double)(naccepted)/nrep;
	}

	double GetAlphaLogPrior()	{
		// exponential prior
		return alpha;
	}

	double GetAlpha()	{
		return alpha;
	}

	int GetAllocationStatistic(int k)	{ // occupation for component k
		int tot = 0;
		for (int i=0; i<GetSize(); i++)	{
			if (GetAllocation(i) == k)	{
				tot++;
			}
		}
		return tot;
	}

	void SwapComponents(int c1, int c2)	{
		//exit(1);
		MixtureRandomMatrix<U>* tempComponent = component[c1];
		component[c1] = component[c2];
		component[c2] = tempComponent;

		for (int i=0; i<GetSize(); i++)	{
			if (allocation[i] == c1)	{
				// SetAllocation(i,c2);
				allocation[i] = c2;
			}
			else if (allocation[i] == c2)	{
				// SetAllocation(i,c1);
				allocation[i] = c1;
			}
		}
	}

	/*
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
	*/

	double Move(double tuning)	{
		double tot = 0;
		tot += MoveAllocation();
		tot += MoveValues(tuning);
		//ResampleWeights();
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
		//for (int i=0; i<this->GetSize(); i++)	{
		//	total += allocation[i]->GetLogProb();
		//}
		//total += weight->GetLogProb();
		return total;
	}
	/*
	double GetWeightEntropy()	{
		return weight->GetEntropy();
	}

	double GetEffSize()	{
		return exp(weight->GetEntropy());
	}
	*/

	void ToStream(ostream& os)	{ //check these...
		//os << *weight << '\n';
		os << alpha << '\n';
		os << GetComponentNumber() << '\n';
		for (int k=0; k<GetComponentNumber(); k++)	{
			os << *(GetComponent(k)) << '\t';
		}
		os << '\n';
		for (int i=0; i<GetSize(); i++)	{
			//os << *allocation[i] << '\t';
			os << allocation[i] << '\t';
		}
		os << '\n';
	}

	void FromStream(istream& is)	{ // check these...

		is >> alpha;
		int initialNcomponent = Ncomponent;
		is >> Ncomponent;
		// make sure that the number of compoents created and the number to be read match
		//
		//cout << "initialNcomponent: " << initialNcomponent << " new Ncomponent: " << Ncomponent << "\n";
		//cout.flush();
		if (initialNcomponent != Ncomponent)	{
			if (initialNcomponent < Ncomponent)	{
				for (int i = initialNcomponent; i < Ncomponent; i++)	{
					component[i] = CreateComponent(i);
				}
			}
			else	{
				for (int i = initialNcomponent-1; i > Ncomponent; i--)	{
					delete component[i];
				}
			}
		}

		//cout << "done creating components\n";
		//cout.flush();

		for (int k=0; k<GetComponentNumber(); k++)	{
			is >> *(GetComponent(k));
		}

		//cout << "done reading components\n";
		//cout.flush();

		for (int i=0; i<GetSize(); i++)	{
			//int temp;			//
			//is >> temp;			// for plain reading, these lines take too much time for nothing
			//cout << "before set allocation\n";
			//cout << "temp: " << temp << ", Ncomponent is: " << Ncomponent << "\n";
			//cout.flush();
			//SetAllocation(i, temp); 	//
			//cout << "after set allocation\n";
			//cout.flush();
			is >> allocation[i];
		}

		//cout << "done matrixmixture FromStream\n";
		//cout.flush();
	}

	friend ostream& operator<<(ostream& os, MatrixInfiniteMixture& mix)	{
		mix.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, MatrixInfiniteMixture& mix)	{
		mix.FromStream(is);
		return is;
	}

	//protected:


	double SetAllocation(int i, int k)	{
		Corrupt(i,true);
		grid[i]->DetachFrom();
		allocation[i] = k;
		grid[i]->AttachTo(GetComponentMatrix(k));
		double temp = Update(i);
		return temp;
		//return Update(i);

	}

	protected:

	void Corrupt(int i, bool bk=false)	{
		grid[i]->Corrupt(bk);
		//allocation[i]->Corrupt(bk);
	}

	double Update(int i)	{
		//return grid[i]->Update() + allocation[i]->Update();
		return grid[i]->Update();
	}
	virtual MixtureRandomMatrix<U>* CreateComponent(int k) = 0;

	void drawSample()	{
		//SampleWeights();
		//SampleAllocation();
		SampleValues();
	}


	int size;
	int NmaxComponents;
	int Ncomponent;
	int* allocation;
	double alpha;
	double alphaMin;
	double alphaMax;
	//FiniteDiscrete** allocation;
	MixtureRandomMatrix<U>** component;
	PointerKeeper** grid;
	//double GetAlphaPrior();
	//Rvar<Profile>* weight;

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


template<class U> class MatInfMixAllocMove : public MCUpdate	{

	public:

	MatInfMixAllocMove(MatrixInfiniteMixture<U>* inmix, int innauxiliary, int innrep) : mix(inmix), nauxiliary(innauxiliary), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double ret = mix->MoveAllocation(nrep, nauxiliary);
		return ret;
	}

	protected:

	MatrixInfiniteMixture<U>* mix;
	int nauxiliary;
	int nrep;
};

template<class U> class MatInfMixAlphaMove : public MCUpdate	{

	public:

	MatInfMixAlphaMove(MatrixInfiniteMixture<U>* inmix, double intuning, int innrep) : mix(inmix), tuning(intuning), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double ret = mix->MoveAlpha(tuning * tuning_modulator, nrep);
		return ret;
	}

	protected:

	MatrixInfiniteMixture<U>* mix;
	double tuning;
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
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		RandomBranchSitePath* path = new  RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(site),0);
		matmix->Register(path,site);
		return path;
	}

	RandomSubMatrix* 	GetMatrix(int site) {return matmix->GetMatrix(site);}

	protected:
	MatrixFiniteMixture<U>* matmix;
};


template<class U> class MatrixInfiniteMixturePhyloProcess : public PhyloProcess	{

	protected:

	public:

	MatrixInfiniteMixturePhyloProcess(LengthTree* intree, MatrixInfiniteMixture<U>* inmatmix,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matmix = inmatmix;
	}

	// the CreateRandomBranchSitePath function tells the PhyloProcess class
	// how to create the substitution process (and the associated substitution path)
	// for a given branch (accessible through link), and a given site
	//
	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		RandomBranchSitePath* path = new  RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()),0,GetMatrix(site),0);
		matmix->Register(path,site);
		return path;
	}

	RandomSubMatrix* 	GetMatrix(int site) {return matmix->GetMatrix(site);}

	protected:
	MatrixInfiniteMixture<U>* matmix;
};




#endif
