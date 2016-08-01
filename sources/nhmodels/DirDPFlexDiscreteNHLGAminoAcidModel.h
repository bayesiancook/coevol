
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
#include "NormalIIDArray.h"

class BinaryAllocationTree : public Rnode {

	public:

	BinaryAllocationTree(BranchVarTree<Profile>* instattree, Var<UnitReal>* inswitchprob, Var<UnitReal>* instat, Dirichlet** inmeanvector, Var<PosReal>* invar)	{
		switchprob = inswitchprob;
		stat = instat;
		meanvector = inmeanvector;
		var = invar;
		stattree = instattree;

		totalweight = new double[2];
		loggammatotalweight = new double[2];
		weight = new double*[2];
		loggammaweight = new double*[2];
		for (int k=0; k<2; k++)	{
			weight[k] = new double[meanvector[0]->GetDim()];
			loggammaweight[k] = new double[meanvector[0]->GetDim()];
		}

		Register(switchprob);
		Register(stat);
		Register(var);
		Register(meanvector[0]);
		Register(meanvector[1]);
		RecursiveRegister(GetRoot());
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return stattree->GetTree();}
	Link* GetRoot() {return GetTree()->GetRoot();}

	int GetBranchAllocation(const Branch* branch)	{
		return alloc[branch];
	}


	double LogProb() { return logProb();}

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

	double GetMeanAlloc()	{
		SampleAlloc();
		int tot = 0;
		double tmp =  RecursiveGetMeanAlloc(GetRoot(),tot);
		return tmp / tot;
	}

	int GetAlloc(const Branch* branch)	{
		return alloc[branch];
	}

	void SampleAlloc()	{
		ComputeLogGammaWeights();
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

	protected:
	
	void drawSample()	{
	}

	double logProb()	{
		ComputeLogGammaWeights();
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

	void ComputeLogGammaWeights()	{

		for (int k=0; k<2; k++)	{
			totalweight[k] = 0;
			for (int i=0; i<meanvector[0]->GetDim(); i++)	{
				weight[k][i] = var->val() * (*(meanvector[k]))[i];
				loggammaweight[k][i] = Random::logGamma(weight[k][i]);
				totalweight[k] += weight[k][i];
			}
			loggammatotalweight[k] = Random::logGamma(totalweight[k]);
		}
	}

	double GetEmissionLogProb(const Link* from, int k, Var<Profile>* stat)	{
		if (clamp[from->GetBranch()] && (k != alloc[from->GetBranch()]))	{
			return -200;
		}
		double total = 0;
		for (int i=0; i<meanvector[0]->GetDim(); i++)	{
			total +=  (weight[k][i] - 1) * log((*stat)[i]);
			total -= loggammaweight[k][i];
		}
		total += loggammatotalweight[k];
		return total;
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
		double total = GetEmissionLogProb(from,k,stattree->GetBranchVal(from->GetBranch()));
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

	int RecursiveGetMeanAlloc(const Link* from, int& tot)	{
		int total = alloc[from->GetBranch()];
		tot++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMeanAlloc(link->Out(),tot);
		}
		return total;
	}

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
		Register(stattree->GetBranchVal(from->GetBranch()));
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
	BranchVarTree<Profile>* stattree;
	Dirichlet** meanvector;
	Var<PosReal>* var;
	map<const Branch*,int> alloc;
	map<const Branch*,bool> clamp;
	map<const Branch*,double*> condL;
	map<const Branch*,double> postdec;
	double** weight;
	double** loggammaweight;
	double* totalweight;
	double* loggammatotalweight;
};

class PseudoDirichlet : public Rvar<Profile>	{

	public: 

	PseudoDirichlet(int indim)	{
		setval(Profile(indim));
		bkvalue = *this;
	}

	double logProb()	{
		return 0;
	}

	virtual double	Move(double tuning, int m)	{
		if (! isClamped())	{
			// Metropolis Hastings here
			Corrupt(true);
			double logHastings = Profile::ProposeMove(tuning, m);
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

	protected:

	void drawSample()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = Random::sGamma(1);
			total += (*this)[i];
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}
};

class BranchPseudoStatTree : public BranchProcess<Profile> {

	public:

	BranchPseudoStatTree(Tree* intree, int indim) : BranchProcess<Profile>(intree,true) {
		dim = indim;
		RecursiveCreate(GetRoot());
	}

	double Move(double tuning, int k)	{
		int n = 0;
		double tot = RecursiveMove(this->GetRoot(),tuning,k,n);
		return tot / n;
	}

	double GetMeanEntropy()	{
		int n = 0;
		double tmp = RecursiveGetTotalEntropy(GetRoot(),n);
		return tmp / n;
	}

	private:

	double RecursiveMove(const Link* from, double tuning, int k, int& count)	{
		double total = GetDirichlet(from)->Move(tuning, k);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveMove(link->Out(),tuning,k,count);
		}
		return total;
	}

	double RecursiveGetTotalEntropy(const Link* from, int& tot)	{
		double total = GetDirichlet(from)->GetEntropy();
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotalEntropy(link->Out(), tot);
		}
		return total;
	}

	Rvar<Profile>* CreateBranchVal(const Link* link)	{
		return new PseudoDirichlet(dim);
	}

	PseudoDirichlet* GetDirichlet(const Link* link) {
		PseudoDirichlet* tmp = dynamic_cast<PseudoDirichlet*>(GetBranchVal(link->GetBranch()));
		if (! tmp)	{
			cerr << "error in stat tree : null pointer\n";
			exit(1);
		}
		return tmp;
	}

	private:
	int dim;
};

class StatTreeMove : public MCUpdate	{

	public:
	
	StatTreeMove(BranchPseudoStatTree* intree, double intuning, int inm) {
		tree = intree;
		tuning = intuning;
		m = inm;
	}
	
	double Move(double tuning_modulator = 1)	{
		return tree->Move(tuning * tuning_modulator,m);
	}

	private:

	BranchPseudoStatTree* tree;
	double tuning;
	int m;
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

	Const<PosReal>* One;
	Const<PosReal>* NAA;

	// tree and branch lengths
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	GammaTree* gamtree;
	
	Gamma* alpha;
	GammaIIDArray* rate;

	Exponential* concentration;
	Dirichlet** center;
	Beta* binswitch;
	Beta* binstat;

	BranchPseudoStatTree* stattree;
	BinaryAllocationTree* alloctree;

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
		
		One = new Const<PosReal>(1.0);
		NAA = new Const<PosReal>(200.0);

		PriorLambda = new Const<PosReal>(10);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,One,lambda);
		lambda->setval(10);
		gamtree->Sample();

		alpha = new Gamma(One,One);
		rate = new GammaIIDArray(Nsite,alpha,alpha);
		alpha->Sample();
		rate->Sample();

		binswitch= new Beta(One,One);
		binstat = new Beta(One,One);

		concentration = new Exponential(NAA,Exponential::MEAN);
		concentration->Sample();
		center = new Dirichlet*[2];
		for (int k=0; k<2; k++)	{
			center[k] = new Dirichlet(Naa);
			center[k]->Sample();
		}
		
		stattree = new BranchPseudoStatTree(tree,Naa);
		alloctree = new BinaryAllocationTree(stattree,binswitch,binstat,center,concentration);
		aamatrixtree = new LGMatrixTree(stattree);

		if (contdata)	{
			alloctree->SetAndClamp(contdata,0);
		}

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixRASPhyloProcess(gamtree, rate, aamatrixtree, proteindata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "root register\n";
		RootRegister(NAA);
		RootRegister(One);
		RootRegister(center[0]);
		RootRegister(center[1]);
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

	BinaryAllocationTree* GetAllocationTree()	{return alloctree;}
	double GetCenter(int k, int i) {return (*center[k])[i];}


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
		
		total += concentration->GetLogProb();
		total += center[0]->GetLogProb();
		total += center[1]->GetLogProb();
		total += stattree->GetLogProb();
		total += alloctree->GetLogProb();

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

		scheduler.Register(new ProfileMove(center[0],0.1,2),10,"center0");
		scheduler.Register(new ProfileMove(center[0],0.03,2),10,"center0");
		scheduler.Register(new ProfileMove(center[0],0.01,5),10,"center0");
		scheduler.Register(new ProfileMove(center[0],0.001,10),10,"center0");

		scheduler.Register(new ProfileMove(center[1],0.1,2),10,"center1");
		scheduler.Register(new ProfileMove(center[1],0.03,2),10,"center1");
		scheduler.Register(new ProfileMove(center[1],0.01,5),10,"center1");
		scheduler.Register(new ProfileMove(center[1],0.001,10),10,"center1");

		scheduler.Register(new StatTreeMove(stattree,1,2),3,"stat tree");
		scheduler.Register(new StatTreeMove(stattree,0.3,2),3,"stat tree");
		scheduler.Register(new StatTreeMove(stattree,0.1,5),3,"stat tree");
		scheduler.Register(new StatTreeMove(stattree,0.01,10),3,"stat tree");

		scheduler.Register(new SimpleMove(binswitch,1),10,"binary switch");
		scheduler.Register(new SimpleMove(binswitch,0.1),10,"binary switch");
		scheduler.Register(new SimpleMove(binswitch,0.01),10,"binary switch");

		scheduler.Register(new SimpleMove(binstat,1),10,"binary stat");
		scheduler.Register(new SimpleMove(binstat,0.1),10,"binary stat");
		scheduler.Register(new SimpleMove(binstat,0.01),10,"binary stat");

		scheduler.Register(new SimpleMove(concentration,1),10,"concentration");
		scheduler.Register(new SimpleMove(concentration,0.1),10,"concentration");
		scheduler.Register(new SimpleMove(concentration,0.01),10,"concentration");

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

		concentration->Sample();
		binswitch->Sample();
		binstat->Sample();
		center[0]->Sample();
		center[1]->Sample();
		stattree->Sample();

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
		os << "#logprior\tlnL\tlength\talpha\tconcentration";
		os << "\tswitch";
		os << "\tstat";
		os << "\ttotalloc";
		os << "\tstatent0";
		os << "\tstatent1";
		os << "\ttreestatent";
		os << '\n';
		os.flush();
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << *alpha;
		os << '\t' << concentration->val();
		os << '\t' << binswitch->val();
		os << '\t' << binstat->val();
		os << '\t' << alloctree->GetTotalAlloc();
		os << '\t' << center[0]->GetEntropy();
		os << '\t' << center[1]->GetEntropy();
		os << '\t' << stattree->GetMeanEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *alpha << '\n';
		os << *rate << '\n';
		os << *concentration << '\n';
		os << *binswitch << '\n';
		os << *binstat << '\n';
		os << *center[0] << '\n';
		os << *center[1] << '\n';
		os << *stattree << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *alpha;
		is >> *rate;
		is >> *concentration;
		is >> *binswitch;
		is >> *binstat;
		is >> *center[0];
		is >> *center[1];
		is >> *stattree;
		
	}

};


