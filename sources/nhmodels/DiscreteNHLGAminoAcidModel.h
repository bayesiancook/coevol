
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
			branchval[0] = CreateBranchAllocation(0);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			branchval[link->GetBranch()] = CreateBranchAllocation(from->GetBranch());
			RecursiveCreate(link->Out());
		}
	}

	Rvar<Int>* CreateBranchAllocation(const Branch* branchup)	{
		if (! branchup)	{
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

	LogStatBinarySwitch(Var<RealVector>* indelta, Var<UnitReal>* instat, Var<Int>* inbinswitch)	{
		setval(RealVector(indelta->GetDim()));
		bkvalue = *this;
		binswitch = inbinswitch;
		stat = instat;
		delta = indelta;
		Register(delta);
		Register(stat);
		Register(binswitch);
	}

	void specialUpdate()	{
		int epsilon =binswitch->val() ? 1 : -1;
		double scale = binswitch->val() ? 1.0 / stat->val() : 1.0 / (1 - stat->val());
		scale *= epsilon;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = scale * (*delta)[i];
		}
	}

	private:
	Var<RealVector>* delta;
	Var<UnitReal>* stat;
	Var<Int>* binswitch;

};

class DeltaLogStatTree : public  BranchValPtrTree< Dvar<RealVector> >	{

	public:

	DeltaLogStatTree(BinaryAllocationTree* inalloctree, Var<RealVector>* indelta)	{
		SetWithRoot(true);
		alloctree = inalloctree;
		delta = indelta;
		stat = inalloctree->GetStat();
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return alloctree->GetTree();}

	private:

	Dvar<RealVector>* CreateBranchVal(const Link* link)	{
		return new LogStatBinarySwitch(delta,stat,alloctree->GetBranchVal(link->GetBranch()));
	}

	BinaryAllocationTree* alloctree;
	Var<RealVector>* delta;
	Var<UnitReal>* stat;
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


class DiscreteNHAminoAcidModel: public ProbModel {

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

	DiscreteNHAminoAcidModel(string datafile, string treefile, string contdatafile, bool inwithsep, bool sample=true, GeneticCodeType type=Universal)	{

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

		deltalogstat = new NormalIIDArray(Naa,Ncont,Zero,One);
		deltalogstattree = new DeltaLogStatTree*[Ncont];
		for (int p=0; p<Ncont; p++)	{
			deltalogstattree[p] = new DeltaLogStatTree(alloctree[p],(*deltalogstat)[p]);
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

	~DiscreteNHAminoAcidModel() {}

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
		total += deltalogstat->GetLogProb();
		for (int p=0; p<Ncont; p++)	{
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

		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.1,2),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.03,2),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.01,5),10,"delta log stat");
		scheduler.Register(new NormalIIDArrayMove(deltalogstat,0.001,10),10,"delta log stat");

		scheduler.Register(new NormalIIDArrayPiecewiseTranslationMove(deltalogstat,1,0,20),10,"delta log stat translation");
		scheduler.Register(new NormalIIDArrayPiecewiseTranslationMove(deltalogstat,0.1,0,20),10,"delta log stat translation");

		for (int p=0; p<Ncont; p++)	{
			scheduler.Register(new SimpleMove(alloctree[p],1),10,"allocations");
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

		binswitch->Sample();
		binstat->Sample();
		deltalogstat->Sample();
		deltalogstat->SetAtZero();
		for (int p=0; p<Ncont; p++)	{
			alloctree[p]->Sample();
		}
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
		os << "#logprior\tlnL\tlength\talpha\tstatent";
		for (int p=0; p<Ncont; p++)	{
			os << "\tswitch" << p;
			os << "\tstat" << p;
			os << "\tmean" << p;
			os << "\tvar" << p;
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
		for (int p=0; p<Ncont; p++)	{
			os << '\t' << *(*binswitch)[p];
			os << '\t' << *(*binstat)[p];
			os << '\t' << (*deltalogstat)[p]->GetMean();
			os << '\t' << (*deltalogstat)[p]->GetVar();
		}
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *alpha << '\n';
		os << *rate << '\n';
		os << *binswitch << '\n';
		os << *binstat << '\n';
		os << *deltalogstat << '\n';
		for (int p=0; p<Ncont; p++)	{
			os << *alloctree[p] << '\n';
		}
		os << *refstat << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *alpha;
		is >> *rate;
		is >> *binswitch;
		is >> *binstat;
		is >> *deltalogstat;
		for (int p=0; p<Ncont; p++)	{
			is >> *alloctree[p];
		}
		is >> *refstat;
		
	}

};


