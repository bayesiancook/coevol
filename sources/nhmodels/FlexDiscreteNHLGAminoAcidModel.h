
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

class BinarySwitch : public Rvar<Int>	{

	public:

	BinarySwitch(BinarySwitch* inup, Var<UnitReal>* inswitchprob, Var<UnitReal>* instat)	{
		up = inup;
		switchprob = inswitchprob;
		stat = instat;
		if (up)	{
			Register(up);
		}
		Register(switchprob);
		Register(stat);
		p = new double[2];
	}

	void drawSample()	{
		updateTransitionProb();
		if (Random::Uniform() > p[0])	{
			setval(1);
		}
		else	{
			setval(0);
		}
	}

	double Move(double tuning = 1)	{
		if (! isClamped())	{
			Corrupt(true);
			setval(1 - val());
			double logratio = Update();
			bool accepted = (log(Random::Uniform()) < logratio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
				return 0;
			}
			return 1;
		}
		return 0;
	}

	private:

	void updateTransitionProb()	{
		if (up)	{
			p[0] = switchprob->val() * (1 - stat->val());	
			p[1] = switchprob->val() * stat->val();
			p[up->val()] += 1 - switchprob->val();
		}
		else	{
			p[0] = 1 - stat->val();	
			p[1] = stat->val();
		}
	}

	double logProb()	{
		updateTransitionProb();
		return log(p[val()]);
	}

	BinarySwitch* up;
	Var<UnitReal>* switchprob;
	Var<UnitReal>* stat;
	double* p;
};

class BinaryAllocationTree : public BranchProcess<Int> {

	public:

	BinaryAllocationTree(Tree* intree, Var<UnitReal>* inswitchprob, Var<UnitReal>* instat)	: BranchProcess<Int>(intree,true) {
		switchprob = inswitchprob;
		stat = instat;
		RecursiveCreate(GetRoot());
	}

	int GetBranchAllocation(const Branch* branch)	{
		return GetBranchVal(branch)->val();
	}

	Var<UnitReal>* GetStat() {return stat;}
	Var<UnitReal>* GetSwitchProb() {return switchprob;}

	void SetAndClamp(ContinuousData* data, int pos) {
		RecursiveSetAndClamp(GetRoot(), data, pos);
	}

	double GetTotal()	{
		return RecursiveGetTotal(GetRoot());
	}

		

	BinarySwitch* GetBinarySwitch(const Branch* branch)	{
		BinarySwitch* tmp = dynamic_cast<BinarySwitch*>(GetBranchVal(branch));
		if (! tmp)	{
			cerr << "error in BinaryAllocationTree: null pointer\n";
			cerr << GetBranchVal(branch) << '\n';
			exit(1);
		}
		return tmp;
	}
	
	protected:
	
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
					GetBinarySwitch(from->GetBranch())->setval((int) tmp);
					GetBinarySwitch(from->GetBranch())->Clamp();
				}
			}
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAndClamp(link->Out(), data, pos);
		}
	}

	double RecursiveGetTotal(const Link* from)	{
		double total = GetBranchVal(from->GetBranch())->val();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotal(link->Out());
		}
		return total;
	}

	Rvar<Int>* CreateBranchVal(const Link* link)	{
		cerr << "error : should not be in create branch val\n";
		exit(1);
		return new BinarySwitch(0,switchprob,stat);
		if (link->isRoot())	{
			return new BinarySwitch(0,switchprob,stat);
		}
		return new BinarySwitch(GetBinarySwitch(link->Out()->GetBranch()),switchprob,stat);
	}

	void RecursiveCreate(const Link* from)	{
		if (from->isRoot())	{
			branchval[0] = CreateBranchAllocation(0,0);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			branchval[link->GetBranch()] = CreateBranchAllocation(from->GetBranch(), link->GetBranch());
			RecursiveCreate(link->Out());
		}
	}

	Rvar<Int>* CreateBranchAllocation(const Branch* branchup, const Branch* branch)	{
		if (! branch)	{
			return new BinarySwitch(0,switchprob,stat);
		}
		return new BinarySwitch(GetBinarySwitch(branchup),switchprob,stat);
	}

	Var<UnitReal>* switchprob;
	Var<UnitReal>* stat;
	Tree* tree;
	
};

class LogStatBinarySwitch : public Dvar<RealVector>	{

	public:

	LogStatBinarySwitch(Var<RealVector>* indelta, Var<Int>* inbinswitch)	{
		setval(RealVector(indelta->GetDim()));
		bkvalue = *this;
		binswitch = inbinswitch;
		delta = indelta;
		Register(delta);
		Register(binswitch);
	}

	void specialUpdate()	{
		// int epsilon = binswitch->val() ? 1 : 0;
		int epsilon = binswitch->val() ? 1 : -1;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = epsilon * (*delta)[i];
		}
	}

	private:
	Var<RealVector>* delta;
	Var<Int>* binswitch;

};

class MeanDeltaLogStatTree : public  BranchValPtrTree< Dvar<RealVector> >	{

	public:

	MeanDeltaLogStatTree(BinaryAllocationTree* inalloctree, Var<RealVector>* indelta)	{
		SetWithRoot(true);
		alloctree = inalloctree;
		delta = indelta;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return alloctree->GetTree();}

	private:

	Dvar<RealVector>* CreateBranchVal(const Link* link)	{
		return new LogStatBinarySwitch(delta,alloctree->GetBranchVal(link->GetBranch()));
	}

	BinaryAllocationTree* alloctree;
	Var<RealVector>* delta;
};

class DeltaLogStatTree : public BranchProcess<RealVector> {

	public:

	DeltaLogStatTree(MeanDeltaLogStatTree* inmeantree, Var<PosReal>* invar) : BranchProcess<RealVector>(inmeantree->GetTree(),true) {
		meantree = inmeantree;
		var = invar;
		RecursiveCreate(GetRoot());
	}

	double PiecewiseTranslationMove(double tuning, int index, int k)	{
		int n = 0;
		double tot = RecursivePiecewiseTranslationMove(this->GetRoot(),tuning,index,k,n);
		return tot / n;
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

	double RecursivePiecewiseTranslationMove(const Link* from, double tuning, int index, int k, int& count)	{
		double total = GetIIDNormal(from)->PiecewiseTranslationMove(tuning, index, k);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursivePiecewiseTranslationMove(link->Out(),tuning,index,k,count);
		}
		return total;
	}

	double RecursiveMove(const Link* from, double tuning, int k, int& count)	{
		double total = GetIIDNormal(from)->Move(tuning, k);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveMove(link->Out(),tuning,k,count);
		}
		return total;
	}

	Rvar<RealVector>* CreateBranchVal(const Link* link)	{
		return new IIDNormal(meantree->GetBranchVal(link->GetBranch()),var);
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

	IIDNormal* GetIIDNormal(const Link* link) {
		IIDNormal* tmp = dynamic_cast<IIDNormal*>(GetBranchVal(link->GetBranch()));
		if (! tmp)	{
			cerr << "error in delta log stat tree : null pointer\n";
			exit(1);
		}
		return tmp;
	}
	
	Var<PosReal>* var;
	MeanDeltaLogStatTree* meantree;
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

class DeltaLogStatTreePiecewiseTranslationMove : public MCUpdate	{

	public:
	
	DeltaLogStatTreePiecewiseTranslationMove(DeltaLogStatTree* intree, double intuning, int inindex, int inm) {
		tree = intree;
		tuning = intuning;
		index = inindex;
		m = inm;
	}
	
	double Move(double tuning_modulator = 1)	{
		return tree->PiecewiseTranslationMove(tuning * tuning_modulator,index,m);
	}

	private:

	DeltaLogStatTree* tree;
	double tuning;
	int m;
	int index;
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
	MeanDeltaLogStatTree** meandeltalogstattree;
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
		for (int p=0; p<Ncont; p++)	{
			alloctree[p] = new BinaryAllocationTree(tree,(*binswitch)[p],(*binstat)[p]);
			// alloctree[p]->GetBinarySwitch(0)->setval(0);
			// alloctree[p]->GetBinarySwitch(0)->Clamp();
		}

		varlogstat = new Gamma(One,One);
		deltalogstat = new NormalIIDArray(Naa,Ncont,Zero,One);
		meandeltalogstattree = new MeanDeltaLogStatTree*[Ncont];
		deltalogstattree = new DeltaLogStatTree*[Ncont];
		for (int p=0; p<Ncont; p++)	{
			meandeltalogstattree[p] = new MeanDeltaLogStatTree(alloctree[p],(*deltalogstat)[p]);
			deltalogstattree[p] = new DeltaLogStatTree(meandeltalogstattree[p],varlogstat);
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
			total += alloctree[p]->GetLogProb();
			total += deltalogstattree[p]->GetLogProb();
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

		scheduler.Register(new NormalIIDArrayPiecewiseTranslationMove(deltalogstat,1,0,20),10,"delta log stat translation");
		scheduler.Register(new NormalIIDArrayPiecewiseTranslationMove(deltalogstat,0.1,0,20),10,"delta log stat translation");

		for (int p=0; p<Ncont; p++)	{
			scheduler.Register(new SimpleMove(alloctree[p],1),10,"allocations");

			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],1,2),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.3,2),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.1,5),10,"delta log stat tree");
			scheduler.Register(new DeltaLogStatTreeMove(deltalogstattree[p],0.01,10),10,"delta log stat tree");

			scheduler.Register(new DeltaLogStatTreePiecewiseTranslationMove(deltalogstattree[p],10,0,20),10,"delta log stat tree translation");
			scheduler.Register(new DeltaLogStatTreePiecewiseTranslationMove(deltalogstattree[p],1,0,20),10,"delta log stat tree translation");
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
		binswitch->Sample();
		binstat->Sample();
		deltalogstat->Sample();
		deltalogstat->SetAtZero();
		for (int p=0; p<Ncont; p++)	{
			cerr << p << '\n';
			cerr << "alloc\n";
			alloctree[p]->Sample();
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
			os << '\t' << alloctree[p]->GetTotal();
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
			os << *alloctree[p] << '\n';
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
			is >> *alloctree[p];
			is >> *deltalogstattree[p];
		}
		is >> *refstat;
		
	}

};


