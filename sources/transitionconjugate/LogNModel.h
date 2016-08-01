
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "LogNormalTreeProcess.h"
#include "Chronogram.h"
#include "GeneralTransitionConjugatePath.h"
#include "F81TransitionMatrix.h"

class F81MatrixTree : public BranchValPtrTree<RandomTransitionMatrix>	{

	public:

	F81MatrixTree(LengthTree* inlengthtree, Var<Profile>* instationary) {
		SetWithRoot(true);
		lengthtree = inlengthtree;
		stationary = instationary;
		cerr << "create f81 matrix tree\n";
		RecursiveCreate(GetRoot());
		cerr << "create f81 matrix tree ok\n";
	}

	~F81MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomTransitionMatrix* CreateBranchVal(const Link* link)	{
		return new RandomF81TransitionMatrix(stationary,lengthtree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	private:

	LengthTree* lengthtree;
	Var<Profile>* stationary;

};


class BranchTransitionMatrixPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchTransitionMatrixPhyloProcess(LengthTree* inlengthtree, BranchValPtrTree<RandomTransitionMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(inlengthtree,indata,false)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		RandomBranchSitePath* tmp = new RandomBranchSitePath(this, matrixtree->GetBranchVal(link->GetBranch()), 0);
		return tmp;
	}

	protected:
	BranchValPtrTree<RandomTransitionMatrix>* matrixtree;
};

class GTRLogNormalModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------
	Const<Real>* Zero;
	Const<PosReal>* One;

	// chronogram 
	Chronogram* chronogram;

	// autocorrelated process
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
	Dirichlet* stationary;
	F81MatrixTree* matrixtree;
	
	// phylo process
	TransitionPathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;
	// OneMatrixPhyloProcess* phyloprocess;
	bool pathconjugate;
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRLogNormalModel(string datafile, string treefile, bool inpathconjugate, bool sample)	{

		pathconjugate = inpathconjugate;

		// fetch data from file
		data = new FileSequenceAlignment(datafile);
		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)

		taxonset = data->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------
		
		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		chronogram = new Chronogram(tree,One);

		// a log normal process on that tree
		sigma = new Gamma(One,One);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		
		// substitution matrix
		stationary = new Dirichlet(Nstate);
		stationary->setuniform();
		// relrate->Clamp();
		// stationary->Clamp();
		matrixtree = new F81MatrixTree(lognormaltree,stationary);

		// a phylogenetic process
		if (pathconjugate)	{
			pathconjtree = new BranchMatrixTransitionPathConjugateTree(chronogram,matrixtree,data);
			phyloprocess = new TransitionPathConjugatePhyloProcess(pathconjtree);
		}
		else	{
			pathconjtree = 0;
			phyloprocess = new BranchTransitionMatrixPhyloProcess(lognormaltree,matrixtree,data);
		}
		cerr << "unfold\n";
		phyloprocess->Unfold();
		if (sample)	{
			cerr << "sample\n";
			phyloprocess->Sample();
		}

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(stationary);

		Register();

		MakeScheduler();

		Update();
		TraceHeader(cerr);
		Trace(cerr);
		
		cerr << "model created\n";

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GTRLogNormalModel() {}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	LogNormalTreeProcess* GetLogNormalProcess()	{
		return lognormaltree;
	}

	Tree* GetTree()	{
		return tree;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	Chronogram* GetChronogram()	{
		return chronogram;
	}

	double GetLogPrior()	{
		double total = 0;
		total += chronogram->GetLogProb();
		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		if (pathconjugate)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		int nrep = 100;
		for (int rep=0; rep<nrep; rep++)	{
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.01),10,"sigma");

		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}

	}
	
	// Metropolis Hastings, successively called on each component of the model
	// old fashioned move
	//
	
	// Draw a sample from the prior

	void drawSample()	{
		chronogram->Sample();
		sigma->Sample();
		lognormaltree->Sample();
		stationary->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}
	
	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tsigma\tmeanrho\tvarrho";
		os << '\n';

	}
	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << sigma->val();
		os << '\t' << GetMeanRho() << '\t' << GetVarRho();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *stationary;
	}

};
