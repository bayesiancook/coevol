
#ifndef BRANCHOMEGAIID_H
#define BRANCHOMEGAIID_H


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"


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

class BranchOmegaIIDModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	SequenceAlignment* codondata;
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
	GammaTree* gamtree;
	
	Dvar<PosReal>* PriorSigma;
	Gamma* sigma;
	GammaTree* omegatree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;
	
	MatrixTree* matrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	BranchOmegaIIDModel(string datafile, string treefile, bool sample=true)	{
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

		// ----------
		// construction of the graph
		// ----------
		
		One = new Const<PosReal>(1);

		// a tree with all branch lengths iid from an exponential distribution of mean meanlength
		// meanlength is itself endowed with an exponential prior of mean 0.1
		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu); 
		gamtree = new GammaTree(tree,One,mu);

		// a log normal process on that tree for the variations of the mutation rate
		PriorSigma = new Const<PosReal>(1);		
		sigma = new Gamma(One,PriorSigma);
		omegatree = new GammaTree(tree,One,sigma);
		
		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree);

		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(gamtree, matrixtree, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
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

	~BranchOmegaIIDModel() {}

	Tree* GetTree() {return tree;}
	GammaTree* GetOmegaTree() {return omegatree;}
	LengthTree* GetLengthTree() {return gamtree;}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += mu->GetLogProb();
		total += gamtree->GetLogProb();

		total += sigma->GetLogProb();
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
		scheduler.Register(new SimpleMove(gamtree,1),10,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"gamtree");

		scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(omegatree,1),10,"omega");
		scheduler.Register(new SimpleMove(omegatree,0.1),10,"omega");
		scheduler.Register(new SimpleMove(omegatree,0.01),10,"omega");

		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new ProfileMove(stationary,1,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.3,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.1,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.1),10,"stat");
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		/*
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,0.1),1,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,0.1),1,"lengthrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,1),1,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,1),1,"lengthrate");
		*/
	}
	
	// Metropolis Hastings, successively called on each component of the model
	// old fashioned move
	//
	
	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
		mu->Sample();
		mu->setval(10);
		gamtree->Sample();
		sigma->Sample();
		sigma->setval(10);
		omegatree->Sample();
		stationary->Sample();
		relrate->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}
	
	double GetVarOmega()	{
		return omegatree->GetVar();
	}

	double GetLength()	{
		return gamtree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tomega\tstderr\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << GetMeanOmega();
		os << '\t' << sqrt(GetVarOmega());
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *gamtree<< '\n';
		os << *sigma << '\n';
		os << *omegatree << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *gamtree;
		is >> *sigma;
		is >> *omegatree;
		is >> *relrate;
		is >> *stationary;
	}

};

#endif

