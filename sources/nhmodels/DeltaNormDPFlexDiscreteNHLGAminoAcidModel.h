
#include "Random.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "BranchProcess.h"
#include "GTRSubMatrix.h"
#include "PrecisionNormalTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "StatProcess.h"
// #include "ACGTProcess.h"
#include "OneMatrixPhyloProcess.h"


class IIDNormal : public virtual Rvar<RealVector>	{

	public:

	IIDNormal(int dim, Var<Real>* inmean, Var<PosReal>*  invariance)	{
		setval(RealVector(dim));
		bkvalue = RealVector(dim);
		mean = inmean;
		meanvector = 0;
		variance = invariance;
		Register(mean);
		Register(variance);
	}

	IIDNormal(Var<RealVector>* inmeanvector, Var<PosReal>*  invariance)	{
		setval(RealVector(inmeanvector->GetDim()));
		bkvalue = RealVector(inmeanvector->GetDim());
		meanvector = inmeanvector;
		mean = 0;
		variance = invariance;
		Register(meanvector);
		Register(variance);
	}

	double logProb()	{
		double total = 0;
		if (mean)	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - mean->val();
				total += tmp * tmp;
			}
		}
		else	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - (*meanvector)[i];
				total += tmp * tmp;
			}
		}
		return -0.5 * GetDim() * log(2 * Pi * variance->val()) -0.5 * total / variance->val();
	}

	double ProposeMove(double tuning, int n)	{
		if ((n<=0) || (n > GetDim()-1))	{
			n = GetDim() -1;
		}
		int* indices = new int[n];
		Random::DrawFromUrn(indices,n,dim);
		for (int i=0; i<n; i++)	{
			vec[indices[i] + 1] += tuning * (Random::Uniform() - 0.5);
		}
		delete[] indices;
		return 0;
	}

	virtual double	Move(double tuning, int m)	{
		if (! isClamped())	{
			// Metropolis Hastings here
			Corrupt(true);
			double logHastings = ProposeMove(tuning, m);
			double deltaLogProb = Update();
			double logRatio = deltaLogProb + logHastings;
			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
			}
			return (double) accepted;
		}
		return 1;
	}

	double GetMean()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += (*this)[i];
		}
		return total / GetDim();
	}
	
	double GetVar()	{
		double mean = 0;
		for (int i=0; i<GetDim(); i++)	{
			mean += (*this)[i];
		}
		double var = 0;
		for (int i=0; i<GetDim(); i++)	{
			var += (*this)[i] * (*this)[i];
		}
		mean /= GetDim();
		var /= GetDim();
		var -= mean * mean;
		return var;
	}

	protected:

	void drawSample()	{
		(*this)[0] = 0;
		if (mean)	{
			for (int i=1; i<GetDim(); i++)	{
				(*this)[i] = Random::sNormal() * sqrt(variance->val()) + mean->val();
			}
		}
		else	{
			for (int i=1; i<GetDim(); i++)	{
				(*this)[i] = Random::sNormal() * sqrt(variance->val()) + (*meanvector)[i];
			}
		}
	}
	private:

	Var<Real>* mean;
	Var<RealVector>* meanvector;
	Var<PosReal>* variance;

};

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

	double GetMean(int index)	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += (*(*this)[i])[index];
		}
		mean /= GetSize();
		return mean;
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

class BinaryAllocationTree : public Rnode {

	public:

	BinaryAllocationTree(BranchVarTree<RealVector>* indeltalogstattree, Var<UnitReal>* inswitchprob, Var<UnitReal>* instat, Var<RealVector>* inmeanvector, Var<PosReal>* invar)	{
		switchprob = inswitchprob;
		stat = instat;
		meanvector = inmeanvector;
		var = invar;
		deltalogstattree = indeltalogstattree;
		Register(switchprob);
		Register(stat);
		Register(var);
		Register(meanvector);
		RecursiveRegister(GetRoot());
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return deltalogstattree->GetTree();}
	Link* GetRoot() {return GetTree()->GetRoot();}

	int GetBranchAllocation(const Branch* branch)	{
		return alloc[branch];
	}


	Var<UnitReal>* GetStat() {return stat;}
	Var<UnitReal>* GetSwitchProb() {return switchprob;}

	void SetAndClamp(ContinuousData* data, int pos) {
		RecursiveSetAndClamp(GetRoot(), data, pos);
	}

	double ProposeMove(double tuning = 1)	{
		return 0;
	}

	double Move(double tuning = 1)	{
		return 1;
	}

	int GetTotalAlloc()	{
		SampleAlloc();
		return RecursiveGetTotalAlloc(GetRoot());
	}

	protected:
	
	void drawSample()	{
	}

	double logProb()	{
		double tot  = 0;
		double* L = condL[0];
		if (! clamp[0])	{
			double max = 0;
			for (int l=0; l<2; l++)	{
				L[l] = Backward(GetRoot(),l);
				if ((!l) || (max < L[l]))	{
					max = L[l];
				}
			}
			tot = (1 - stat->val()) * exp(L[0] - max) + stat->val() * exp(L[1] - max) + max;
		}
		else	{
			if (alloc[0])	{
				tot = log(stat->val()) +  Backward(GetRoot(),0);
			}
			else	{
				tot = log(1 - stat->val()) + Backward(GetRoot(),1);
			}
		}
		if (isinf(tot))	{
			cerr << "null prob : " << stat->val() << '\t' << L[0] << '\t' << L[1] << '\n';
			exit(1);
		}
		if (isnan(tot))	{
			cerr << "null prob : " << stat->val() << '\t' << L[0] << '\t' << L[1] << '\n';
			exit(1);
		}
		return tot;
	}

	double GetTransitionLogProb(int k,int l)	{
		double total = 0;
		if (k==l)	{
			total += 1 - switchprob->val();
		}
		if (l)	{
			total += switchprob->val() * stat->val();
		}
		else	{
			total += switchprob->val() * (1 - stat->val());
		}
		if (total <=0)	{
			cerr << "negative transition prob : " << total << '\n';
			cerr << switchprob->val() << '\t' << stat->val() << '\n';
			exit(1);
		}
		return log(total);
	}

	double GetEmissionLogProb(const Link* from, int k, Var<RealVector>* deltalogstat)	{
		double total = 0;
		if (clamp[from->GetBranch()] && (k != alloc[from->GetBranch()]))	{
			return -200;
		}
		if (k)	{
			// normal centered on deltalogstat, and var var
			for (int i=0; i<meanvector->GetDim(); i++)	{
				double tmp = (*deltalogstat)[i] - (*meanvector)[i];
				total += tmp * tmp;
			}
		}
		else	{
			// normal centered on 0
			for (int i=0; i<meanvector->GetDim(); i++)	{
				double tmp = (*deltalogstat)[i];
				total += tmp * tmp;
			}
		}
		double tmp = -0.5 * meanvector->GetDim() * log(var->val()) -0.5 * total / var->val();
		return tmp;
		// return exp(tmp);
	}

	void SampleAlloc()	{
		if (! clamp[0])	{
			double* L = condL[0];
			double p0 = (1 - stat->val()) * exp(L[0]);
			double p1 = (stat->val()) * exp(L[1]);
			double u = (p0 + p1) * Random::Uniform();
			if (u<p0)	{
				alloc[0] = 0;
			}
			else	{
				alloc[0] = 1;
			}
		}
		Forward(GetRoot(),alloc[0]);
	}

	void Forward(const Link* from, int k)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! clamp[link->GetBranch()])	{
				double* L = condL[link->GetBranch()];
				double p0 = exp(GetTransitionLogProb(k,0) + L[0]);
				double p1 = exp(GetTransitionLogProb(k,1) + L[1]);
				double u = (p0 + p1) * Random::Uniform();
				if (u<p0)	{
					alloc[link->GetBranch()] = 0;
				}
				else	{
					alloc[link->GetBranch()] = 1;
				}
			}
			Forward(link->Out(),alloc[link->GetBranch()]);
		}
	}

	double Backward(const Link* from, int k)	{
		double total = GetEmissionLogProb(from,k,deltalogstattree->GetBranchVal(from->GetBranch()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			double* L = condL[link->GetBranch()];
			double max = 0;
			for (int l=0; l<2; l++)	{
				L[l] = GetTransitionLogProb(k,l) + Backward(link->Out(),l);
				if ((!l) || (max < L[l]))	{
					max = L[l];
				}
			}
			double tot = 0;
			for (int l=0; l<2; l++)	{
				L[l] -= max;
				tot += exp(L[l]);
			}
			total += log(tot) + max;
		}
		return total;
	}

	int RecursiveGetTotalAlloc(const Link* from)	{
		int total = alloc[from->GetBranch()];
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotalAlloc(link->Out());
		}
		return total;
	}

	/*
	double Backward(const Link* from, int k)	{
		// offset[from->GetBranch()] = 1;
		double total = GetEmissionProb(from,k,deltalogstattree->GetBranchVal(from->GetBranch()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			double* L = condL[link->GetBranch()];
			for (int l=0; l<2; l++)	{
				L[l] = Backward(link->Out(),l);
			}
			double min = 0;
			if (L[0] > 0)	{
				min = L[0];
			}
			if (L[1] > 0)	{
				if (min > L[1])	{
					min = L[1];
				}
			}
			offset[link->GetBranch()] *= min;
			double tot = 0;
			for (int l=0; l<2; l++)	{
				L[l] /= min;
				L[l] *= GetTransitionProb(k,l);
				tot += L[l];
			}
			total *= tot;
		}
		return total;
	}
	*/

	void RecursiveSetAndClamp(const Link* from, ContinuousData* data, int pos){
		if(from->isLeaf()){
			int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			cerr << "set and clamp : " << from->GetNode()->GetName() << '\t' << tax << '\t' << pos << '\t';
			if (tax != -1)	{
				double tmp = data->GetState(tax,pos);
				cerr << tmp << '\n';
				if (tmp != -1)	{
					if ((tmp != 0) && (tmp != 1))	{
						cerr << "error in recursive clamp at: negative cont data\n";
						exit(1);
					}
					clamp[from->GetBranch()] = true;
					alloc[from->GetBranch()] = (int) tmp;
				}
			}
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAndClamp(link->Out(), data, pos);
		}
	}

	void RecursiveRegister(const Link* from)	{
		Register(deltalogstattree->GetBranchVal(from->GetBranch()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	void RecursiveCreate(const Link* from)	{
		if (from->isRoot())	{
			alloc[0] = -1;
			clamp[0] = false;
			condL[0] = new double[2];
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			alloc[link->GetBranch()] = -1;
			clamp[link->GetBranch()] = false;
			condL[link->GetBranch()] = new double[2];
			RecursiveCreate(link->Out());
		}
	}

	Var<UnitReal>* switchprob;
	Var<UnitReal>* stat;
	BranchVarTree<RealVector>* deltalogstattree;
	Var<RealVector>* meanvector;
	Var<PosReal>* var;
	map<const Branch*,int> alloc;
	map<const Branch*,bool> clamp;
	map<const Branch*,double*> condL;
	map<const Branch*,double> postdec;
};

class PseudoIIDNormal : public Rvar<RealVector>	{

	public: 

	PseudoIIDNormal(int indim)	{
		setval(RealVector(indim));
		bkvalue = *this;
	}

	double logProb()	{
		return 0;
	}

	double ProposeMove(double tuning, int n)	{
		if ((n<=0) || (n > GetDim()-1))	{
			n = GetDim() -1;
		}
		int* indices = new int[n];
		Random::DrawFromUrn(indices,n,dim);
		for (int i=0; i<n; i++)	{
			vec[indices[i] + 1] += tuning * (Random::Uniform() - 0.5);
		}
		delete[] indices;
		return 0;
	}

	virtual double	Move(double tuning, int m)	{
		if (! isClamped())	{
			// Metropolis Hastings here
			Corrupt(true);
			double logHastings = ProposeMove(tuning, m);
			double deltaLogProb = Update();
			double logRatio = deltaLogProb + logHastings;
			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
			}
			return (double) accepted;
		}
		return 1;
	}

	double GetMean()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += (*this)[i];
		}
		return total / GetDim();
	}
	
	double GetVar()	{
		double mean = 0;
		for (int i=0; i<GetDim(); i++)	{
			mean += (*this)[i];
		}
		double var = 0;
		for (int i=0; i<GetDim(); i++)	{
			var += (*this)[i] * (*this)[i];
		}
		mean /= GetDim();
		var /= GetDim();
		var -= mean * mean;
		return var;
	}


	protected:

	void drawSample()	{
		(*this)[0] = 0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = 0;
		}
	}
};

class DeltaLogStatTree : public BranchProcess<RealVector> {

	public:

	DeltaLogStatTree(Tree* intree, int indim) : BranchProcess<RealVector>(intree,true) {
		dim = indim;
		RecursiveCreate(GetRoot());
	}

	double Move(double tuning, int k)	{
		int n = 0;
		double tot = RecursiveMove(this->GetRoot(),tuning,k,n);
		return tot / n;
	}

	double GetGrandMean()	{
		int n = 0;
		double tmp = RecursiveGetGrandMean(GetRoot(),n);
		return tmp / n;
	}

	double GetGrandVar()	{
		int n = 0;
		double tmp = RecursiveGetGrandVar(GetRoot(),n);
		return tmp / n;
	}

	private:

	double RecursiveMove(const Link* from, double tuning, int k, int& count)	{
		double total = GetIIDNormal(from)->Move(tuning, k);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveMove(link->Out(),tuning,k,count);
		}
		return total;
	}

	double RecursiveGetGrandVar(const Link* from, int& tot)	{
		double total = GetIIDNormal(from)->GetVar();
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetGrandVar(link->Out(), tot);
		}
		return total;
	}

	double RecursiveGetGrandMean(const Link* from, int& tot)	{
		double total = GetIIDNormal(from)->GetMean();
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetGrandMean(link->Out(), tot);
		}
		return total;
	}

	Rvar<RealVector>* CreateBranchVal(const Link* link)	{
		return new PseudoIIDNormal(dim);
	}

	PseudoIIDNormal* GetIIDNormal(const Link* link) {
		PseudoIIDNormal* tmp = dynamic_cast<PseudoIIDNormal*>(GetBranchVal(link->GetBranch()));
		if (! tmp)	{
			cerr << "error in delta log stat tree : null pointer\n";
			exit(1);
		}
		return tmp;
	}

	private:
	int dim;
};

class DeltaLogStatTreeMove : public MCUpdate	{

	public:
	
	DeltaLogStatTreeMove(DeltaLogStatTree* intree, double intuning, int inm) {
		tree = intree;
		tuning = intuning;
		m = inm;
	}
	
	double Move(double tuning_modulator = 1)	{
		return tree->Move(tuning * tuning_modulator,m);
	}

	private:

	DeltaLogStatTree* tree;
	double tuning;
	int m;
};

class DiscreteStat : public Dvar<Profile>	{

	public:
	
	DiscreteStat(Var<Profile>* inrefstat, Var<RealVector>** inlist, int inP)	{
		setval(Profile(inrefstat->GetDim()));
		bkvalue = *this;
		P = inP;
		refstat = inrefstat;
		Register(refstat);
		logstatlist = new Var<RealVector>*[P];
		for (int p=0; p<P; p++)	{
			logstatlist[p] = inlist[p];
			Register(logstatlist[p]);
		}
	}

	void specialUpdate()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			double tmp = 0;
			for (int p=0; p<P; p++)	{
				tmp += (*(logstatlist[p]))[i];
			}
			(*this)[i] = (*refstat)[i] * exp(tmp);
			total += (*this)[i];
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}

	private:

	int P;
	Var<RealVector>** logstatlist;
	Var<Profile>* refstat;

};

class DiscreteStatTree : public BranchValPtrTree<Dvar<Profile> > {

	public:

	DiscreteStatTree(Var<Profile>* inrefstat, DeltaLogStatTree** inlist, int inP)	{
		SetWithRoot(true);
		P = inP;
		logstattreelist = inlist;
		refstat = inrefstat;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return logstattreelist[0]->GetTree();}

	private:
	
	Dvar<Profile>* CreateBranchVal(const Link* link)	{
		Var<RealVector>** tmp = new Var<RealVector>*[P];
		for (int p=0; p<P; p++)	{
			tmp[p] = logstattreelist[p]->GetBranchVal(link->GetBranch());
		}
		DiscreteStat* ret = new DiscreteStat(refstat,tmp,P);
		delete[] tmp;
		return ret;
	}

	int P;
	Var<Profile>* refstat;
	DeltaLogStatTree** logstattreelist;
	
};


class LGMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	LGMatrixTree(BranchVarTree<Profile>* instattree) {
		SetWithRoot(true);
		stattree = instattree;
		RecursiveCreate(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new LGRandomSubMatrix(stattree->GetBranchVal(link->GetBranch()),true);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;

};

class BranchMatrixRASPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
		rate = inrate;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), GetRate(site), matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:

	Var<PosReal>* GetRate(int site)	{
		if (! rate)	{
			return 0;
		}
		return rate->GetVal(site);
	}

	BranchValPtrTree<RandomSubMatrix>* matrixtree;
	VarArray<PosReal>* rate;
};


class FlexDiscreteNHAminoAcidModel: public ProbModel {

	bool withsep;

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	ProteinSequenceAlignment* proteindata;
	ContinuousData* contdata;
	GCContinuousData* gcdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ntaxa;

	public:
	int Ncont;

	private:
	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	// tree and branch lengths
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	GammaTree* gamtree;
	
	Gamma* alpha;
	GammaIIDArray* rate;

	Gamma* varlogstat;
	NormalIIDArray* deltalogstat;
	BetaIIDArray* binswitch;
	BetaIIDArray* binstat;
	BinaryAllocationTree** alloctree;

	Dirichlet* refstat;
	DeltaLogStatTree** deltalogstattree;
	DiscreteStatTree* discretestattree;
	LGMatrixTree* aamatrixtree;
	
	BranchMatrixRASPhyloProcess* phyloprocess;

	public :

	// constructor
	// this is where the entire graph structure of the model is created

	FlexDiscreteNHAminoAcidModel(string datafile, string treefile, string contdatafile, bool inwithsep, bool sample=true, GeneticCodeType type=Universal)	{

		withsep = inwithsep;

		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();
		Ntaxa = taxonset->GetNtaxa();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "protein data\n";
		proteindata = new ProteinSequenceAlignment(codondata);

		cerr << "gc data\n";
		gcdata = new GCContinuousData(codondata,2);

		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------
		
		Zero= new Const<Real>(0);
		One = new Const<PosReal>(1);

		PriorLambda = new Const<PosReal>(10);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,One,lambda);
		lambda->setval(10);
		gamtree->Sample();

		alpha = new Gamma(One,One);
		rate = new GammaIIDArray(Nsite,alpha,alpha);
		alpha->Sample();
		rate->Sample();

		binswitch= new BetaIIDArray(Ncont,One,One);
		binstat = new BetaIIDArray(Ncont,One,One);
		alloctree = new BinaryAllocationTree*[Ncont];

		varlogstat = new Gamma(One,One);
		deltalogstat = new NormalIIDArray(Naa,Ncont,Zero,One);
		deltalogstattree = new DeltaLogStatTree*[Ncont];
		for (int p=0; p<Ncont; p++)	{
			deltalogstattree[p] = new DeltaLogStatTree(tree,Naa);
			alloctree[p] = new BinaryAllocationTree(deltalogstattree[p],(*binswitch)[p],(*binstat)[p],(*deltalogstat)[p],varlogstat);
		}

		refstat = new Dirichlet(Naa);
		refstat->setuniform();
		cerr << "stattree\n";
		discretestattree = new DiscreteStatTree(refstat,deltalogstattree,Ncont);
		aamatrixtree = new LGMatrixTree(discretestattree);

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				alloctree[i]->SetAndClamp(contdata,i);
			}
		}

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixRASPhyloProcess(gamtree, rate, aamatrixtree, proteindata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "root register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(refstat);
		RootRegister(PriorLambda);
		cerr << "register\n";
		Register();
		cerr << "register ok\n";

		MakeScheduler();

		if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
			cerr << '\n';
			GetTree()->Print(cerr);
			cerr << '\n';
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~FlexDiscreteNHAminoAcidModel() {}

	Tree* GetTree() {return tree;}
	LengthTree* GetGamTree() {return gamtree;}

	int GetNtaxa() {return Ntaxa;}
	int GetNstate() {return Nstate;}

	BinaryAllocationTree* GetAllocationTree(int k)	{return alloctree[k];}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();

		total += alpha->GetLogProb();
		total += rate->GetLogProb();

		total += binswitch->GetLogProb();
		total += binstat->GetLogProb();
		
		total += varlogstat->GetLogProb();
		total += deltalogstat->GetLogProb();
		for (int p=0; p<Ncont; p++)	{
			total += deltalogstattree[p]->GetLogProb();
			total += alloctree[p]->GetLogProb();
		}

		total += refstat->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),30,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.1),30,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.01),30,"gamtree");

		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(rate,3),10,"rates across sites");
		scheduler.Register(new SimpleMove(rate,0.3),10,"rates across sites");

		scheduler.Register(new SimpleMove(binswitch,1),10,"binary switch");
		scheduler.Register(new SimpleMove(binswitch,0.1),10,"binary switch");
		scheduler.Register(new SimpleMove(binswitch,0.01),10,"binary switch");

		scheduler.Register(new SimpleMove(binstat,1),10,"binary stat");
		scheduler.Register(new SimpleMove(binstat,0.1),10,"binary stat");
		scheduler.Register(new SimpleMove(binstat,0.01),10,"binary stat");

		scheduler.Register(new SimpleMove(varlogstat,1),10,"binary stat");
		scheduler.Register(new SimpleMove(varlogstat,0.1),10,"binary stat");
		scheduler.Register(new SimpleMove(varlogstat,0.01),10,"binary stat");

		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.1,2),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.03,2),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.01,5),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.001,10),10,"delta log stat");

		for (int p=0; p<Ncont; p++)	{
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],1,2),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.3,2),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.1,5),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.01,10),10,"delta log stat tree");
		}

		scheduler.Register(new ProfileMove(refstat,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(refstat,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(refstat,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(refstat,0.001),10,"stat");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
		// lambda->Sample();
		cerr << "gamtree\n";
		// gamtree->Sample();

		alpha->Sample();
		alpha->setval(10);
		rate->Sample();

		cerr << "varlogstat\n";
		varlogstat->Sample();
		cerr << *varlogstat << '\n';
		//varlogstat->setval(1);
		binswitch->Sample();
		binstat->Sample();
		/*
		for (int p=0; p<Ncont; p++)	{
			(*binswitch)[p]->setval(0.5);
			(*binstat)[p]->setval(0.5);
		}
		*/
		deltalogstat->Sample();
		deltalogstat->SetAtZero();
		for (int p=0; p<Ncont; p++)	{
			cerr << p << '\n';
			cerr << "deltalogstattree\n";
			deltalogstattree[p]->Sample();
		}
		cerr << "refstat\n";
		refstat->Sample();
		cerr << "phyloprocess\n";
		phyloprocess->Sample();
		cerr << "ok\n";
	}

	// various summary statistics
	// used to check mcmc convergence

	double GetLength()	{
		return gamtree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\talpha\tstatent\tvarlogstat";
		for (int p=0; p<Ncont; p++)	{
			os << "\tswitch" << p;
			os << "\tstat" << p;
			os << "\ttotalloc" << p;
			os << "\tmean" << p;
			os << "\tvar" << p;
			os << "\ttreemean" << p;
			os << "\ttreevar" << p;
		}
		os << '\n';
		os.flush();
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << *alpha;
		os << '\t' << refstat->GetEntropy();
		os << '\t' << varlogstat->val();
		for (int p=0; p<Ncont; p++)	{
			os << '\t' << *(*binswitch)[p];
			os << '\t' << *(*binstat)[p];
			os << '\t' << alloctree[p]->GetTotalAlloc();
			os << '\t' << (*deltalogstat)[p]->GetMean();
			os << '\t' << (*deltalogstat)[p]->GetVar();
			os << '\t' << deltalogstattree[p]->GetGrandMean();
			os << '\t' << deltalogstattree[p]->GetGrandVar();
		}
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *alpha << '\n';
		os << *rate << '\n';
		os << *varlogstat << '\n';
		os << *binswitch << '\n';
		os << *binstat << '\n';
		os << *deltalogstat << '\n';
		for (int p=0; p<Ncont; p++)	{
			// os << *alloctree[p] << '\n';
			os << *deltalogstattree[p] << '\n';
		}
		os << *refstat << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *alpha;
		is >> *rate;
		is >> *varlogstat;
		is >> *binswitch;
		is >> *binstat;
		is >> *deltalogstat;
		for (int p=0; p<Ncont; p++)	{
			// is >> *alloctree[p];
			is >> *deltalogstattree[p];
		}
		is >> *refstat;
		
	}

};


