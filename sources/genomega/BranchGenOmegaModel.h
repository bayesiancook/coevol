
// a chronogram
// 3 independent log normal processes
// - log tau (generation time)
// - log Ks per generation
// - log omega = log Kn/Ks
//
// first 2 components are combined to obtain log Ks per Myr

#ifndef BRANCHOMEGA_H
#define BRANCHOMEGA_H


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "Chronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"


#include "MeanExpTree.h"


class LogSynRate : public Dvar<Real>	{

	public:

	LogSynRate(Var<Real>* inlogmu, Var<Real>* inloggentime)	{
		logmu = inlogmu;
		loggentime = inloggentime;
		Register(logmu);
		Register(loggentime);
	}

	protected:

	void specialUpdate()	{
		setval(logmu->val() - 0.5 * loggentime->val());
	}

	private:

	Var<Real>* logmu;
	Var<Real>* loggentime;
};

class LogSynRateTree : public NodeValPtrTree<Dvar<Real> >	{

	public:

	LogSynRateTree(LogNormalTreeProcess* inmutree, LogNormalTreeProcess* ingentimetree) {
		mutree = inmutree;
		gentimetree = ingentimetree;
		RecursiveCreate(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	protected:

	void specialUpdate(Link* from)	{
		GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
		
	Dvar<Real>* CreateNodeVal(const Link* link)	{
		return new LogSynRate(mutree->GetNodeVal(link->GetNode()), gentimetree->GetNodeVal(link->GetNode()));
	}

	Tree* GetTree() {return mutree->GetTree();}

	private:

	LogNormalTreeProcess* mutree;
	LogNormalTreeProcess* gentimetree;

};

class MatrixTree : public BranchValPtrTree<RandomMGOmegaCodonSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = new Const<PosReal>(1.0);
		RecursiveCreate(GetRoot());
	}

	Var<PosReal>* GetRootOmega() {return rootomega;}

	protected:

	RandomMGOmegaCodonSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
	Const<PosReal>* rootomega;

};

class BranchMatrixPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchMatrixPhyloProcess(LengthTree* intree, BranchValPtrTree<RandomMGOmegaCodonSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:
	BranchValPtrTree<RandomMGOmegaCodonSubMatrix>* matrixtree;
};


// Compensatory compensatory move functions
class LogNormalAntiCompensatoryMove : public MCUpdate, public NodeValPtrTree<Mnode> {

	public:

	LogNormalAntiCompensatoryMove(LogNormalTreeProcess* ina1, LogNormalTreeProcess* ina2, double intuning)	{
		tree = ina1->GetTree();
		tuning = intuning;
		a1 = ina1;
		a2 = ina2;
		RecursiveCreate(GetRoot());
		RecursiveRegister(GetRoot());
	}

	Tree* GetTree() {return tree;}

	double Move(double tuning_modulator = 1)	{
		int tot = 0;
		double total = RecursiveMove(GetRoot(),tuning*tuning_modulator,tot);
		return total/tot;
	}

	Mnode* CreateNodeVal(const Link* from)	{
		return new Mnode;
	}

	private:

	double RecursiveMove(Link* from, double tuning, int& tot)	{

		tot++;
		double totaccept = 0;

		Mnode* m = GetNodeVal(from->GetNode());
		Rvar<Real>* b1 = a1->GetNodeVal(from->GetNode());
		Rvar<Real>* b2 = a2->GetNodeVal(from->GetNode());

		m->Corrupt(true);

		double h = tuning * (Random::Uniform() - 0.5);
		
		if (b1->isClamped())	{
			b1->ScalarAddition(0);
		}
		else	{
			b1->ScalarAddition(h);
		}
		if (b2->isClamped())	{
			b2->ScalarAddition(0);
		}
		else	{
			b2->ScalarAddition(h);
		}
		
		double logratio = m->Update();
		bool accepted = (log(Random::Uniform()) < logratio);

		if (! accepted)	{
			m->Corrupt(false);
			m->Restore();
			totaccept++;
		}
		else	{
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			totaccept += RecursiveMove(link->Out(),tuning,tot);
		}
	
		return (double) totaccept;
	}

	void RecursiveRegister(Link* from)	{
		a1->GetNodeVal(from->GetNode())->Register(GetNodeVal(from->GetNode()));	
		a2->GetNodeVal(from->GetNode())->Register(GetNodeVal(from->GetNode()));	
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	LogNormalTreeProcess* a1;
	LogNormalTreeProcess* a2;
	Tree* tree;
	double tuning;
};

class BranchGenOmegaModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	SequenceAlignment* codondata;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;

	// chronogram 
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;
	
	// mutation rate (per generation)
	Dvar<PosReal>* PriorSigma;
	Gamma* sigma;
	LogNormalTreeProcess* synratepergentree;

	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* gentimetree;

	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* omegatree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;
	
	LogSynRateTree* logsynratetree;
	MeanExpTree* synratetree;
	MatrixTree* matrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	BranchGenOmegaModel(string datafile, string treefile, string contdatafile, bool sample=true)	{
		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		contdata = new FileContinuousData(contdatafile);

		// ----------
		// construction of the graph
		// ----------
		
		One = new Const<PosReal>(1);

		// a tree with all branch lengths iid from an exponential distribution of mean meanlength
		// meanlength is itself endowed with an exponential prior of mean 0.1
		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu); 
		mu->ClampAt(1);
		chronogram = new Chronogram(tree,mu);

		// a log normal process on that tree for the variations of the mutation rate
		PriorSigma = new Const<PosReal>(1);		
		sigma = new Gamma(One,PriorSigma);
		gentimetree = new LogNormalTreeProcess(chronogram,sigma,MEAN);
		gentimetree->ClampAt(contdata,0);
		
		// a log normal process on that tree for the variations of the mutation rate
		PriorTheta = new Const<PosReal>(1);		
		theta = new Gamma(One,PriorTheta);
		synratepergentree = new LogNormalTreeProcess(chronogram,theta,MEAN);

		logsynratetree = new LogSynRateTree(synratepergentree,gentimetree);
		synratetree = new MeanExpTree(logsynratetree,chronogram,INTEGRAL);
		
		// another log normal process for the variations of pop size
		PriorTau = new Const<PosReal>(1);		
		tau = new Gamma(One,PriorTau);
		omegatree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		
		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree);

		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(PriorTheta);
		RootRegister(PriorTau);
		RootRegister(synratepergentree->GetRootRate());
		RootRegister(gentimetree->GetRootRate());
		RootRegister(omegatree->GetRootRate());
		RootRegister(matrixtree->GetRootOmega());
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		MakeScheduler();

		if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~BranchGenOmegaModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetOmegaTree() {return omegatree;}
	LogNormalTreeProcess* GetSynRatePerGenTree() {return synratepergentree;}
	LogSynRateTree* GetLogSynRateTree() {return logsynratetree;}
	MeanExpTree* GetSynRateTree() {return synratetree;}
	LogNormalTreeProcess* GetGenTimeTree() {return gentimetree;}
	LengthTree* GetChronogram() {return chronogram;}

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

		total += mu->GetLogProb();
		total += chronogram->GetLogProb();

		total += sigma->GetLogProb();
		total += gentimetree->GetLogProb();

		total += theta->GetLogProb();
		total += synratepergentree->GetLogProb();

		total += tau->GetLogProb();
		total += omegatree->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{
	
		scheduler.Register(new SimpleMove(mu,1),10,"mu");
		scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(gentimetree,10),40,"loggentime");
		scheduler.Register(new SimpleMove(gentimetree,1),40,"loggentime");
		// scheduler.Register(new SimpleMove(gentimetree,0.1),10,"loggentime");

		scheduler.Register(new SimpleMove(theta,1),100,"sigma");
		scheduler.Register(new SimpleMove(theta,0.1),100,"sigma");
		// scheduler.Register(new SimpleMove(synratepergentree,10),10,"lograte");
		scheduler.Register(new SimpleMove(synratepergentree,1),40,"lograte");
		scheduler.Register(new SimpleMove(synratepergentree,0.1),40,"lograte");

		scheduler.Register(new LogNormalAntiCompensatoryMove(synratepergentree,gentimetree,1),10,"rategentime");
		scheduler.Register(new LogNormalAntiCompensatoryMove(synratepergentree,gentimetree,0.1),10,"rategentime");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		// scheduler.Register(new SimpleMove(omegatree,10),10,"logomega");
		scheduler.Register(new SimpleMove(omegatree,1),10,"logomega");
		scheduler.Register(new SimpleMove(omegatree,0.1),10,"logomega");

		scheduler.Register(new SimpleMove(relrate,0.03),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		// scheduler.Register(new ProfileMove(stationary,0.1,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

		scheduler.Register(new AdditiveAntiCompensatoryMove(gentimetree,synratepergentree,0.3),10,"rategentime");
		scheduler.Register(new AdditiveAntiCompensatoryMove(gentimetree,synratepergentree,0.03),30,"rategentime");
		scheduler.Register(new AdditiveAntiCompensatoryMove(gentimetree,synratepergentree,0.005),30,"rategentime");
		// scheduler.Register(new AdditiveAntiCompensatoryMove(gentimetree,synratepergentree,0.01),4,"rategentime");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	// Metropolis Hastings, successively called on each component of the model
	// old fashioned move
	//
	
	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
		mu->Sample();
		// chronogram->Sample();
		sigma->Sample();
		sigma->setval(10);
		gentimetree->Sample();

		theta->Sample();
		theta->setval(10);
		synratepergentree->Sample();

		tau->Sample();
		tau->setval(10);
		omegatree->Sample();

		stationary->Sample();
		relrate->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanOmega()	{
		return omegatree->GetMeanRate();
	}

	double GetVarOmega()	{
		return omegatree->GetVarRate();
	}

	double GetMeanSynRatePerGen()	{
		return synratepergentree->GetMeanRate();
	}
	
	double GetVarSynRatePerGen()	{
		return synratepergentree->GetVarRate();
	}

	double GetMeanGenTime()	{
		return gentimetree->GetMeanRate();
	}
	
	double GetVarGenTime()	{
		return gentimetree->GetVarRate();
	}

	double GetMeanLogGenTime()	{
		return gentimetree->GetMeanLogRate();
	}
	
	double GetVarLogGenTime()	{
		return gentimetree->GetVarLogRate();
	}

	double GetLength()	{
		return synratetree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength(Ls)\tdS/gen\tstderr\tomega\tstderr\tgentime\tstderr\tloggentime\tsigma\ttheta\ttau\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << GetMeanSynRatePerGen();
		os << '\t' << sqrt(GetVarSynRatePerGen());
		os << '\t' << GetMeanOmega();
		os << '\t' << sqrt(GetVarOmega());
		os << '\t' << GetMeanLogGenTime();
		os << '\t' << sqrt(GetVarLogGenTime());
		os << '\t' << GetMeanGenTime();
		os << '\t' << sigma->val();
		os << '\t' << theta->val();
		os << '\t' << tau->val();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *synratepergentree << '\n';
		os << *theta << '\n';
		os << *gentimetree << '\n';
		os << *tau << '\n';
		os << *omegatree << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *synratepergentree;
		is >> *theta;
		is >> *gentimetree;
		is >> *tau;
		is >> *omegatree;
		is >> *relrate;
		is >> *stationary;
	}

};

#endif

