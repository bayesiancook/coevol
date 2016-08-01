
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugateOneMatrixPhyloProcess.h"
#include "ConjugateGammaTree.h"
#include "Move.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"

class MGModel : public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	MCScheduler scheduler;

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

	// tree and branch lengths
	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	GammaTree* gamtree;

	// nucleotide relative exchange rates
	IIDExp* relrate;
	Dirichlet* stationary;

	// nucleotide substitution matrix is relrate * stationary
	GTRRandomSubMatrix* matrix;

	// one single global omega parameter
	// with an exponential prior of mean 1
	Dvar<PosReal>* PriorOmega;
	Exponential* omega;

	RandomMGOmegaCodonSubMatrix* codonmatrix;

	OneMatrixPhyloProcess* phyloprocess;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	MGModel(string datafile, string treefile) : scheduler(this)	{
		// fetch data from file

		nucdata = new FileSequenceAlignment(datafile);

		codondata = new CodonSequenceAlignment(nucdata, true);

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		cerr << Nsite << '\t' << Nstate << '\n';

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------

		// a tree with all branch lengths iid from an exponential distribution of rate lambda
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,mu,lambda);

		// nucleotide substitution matrix
		relrate = new IIDExp(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		matrix = new GTRRandomSubMatrix(relrate,stationary);

		// omega
		PriorOmega = new Const<PosReal>(1);
		omega = new Exponential(PriorOmega, Exponential::MEAN);

		// codon matrix
		codonmatrix = new RandomMGOmegaCodonSubMatrix((CodonStateSpace*) codondata->GetStateSpace(),matrix,omega, true);

		cerr << "create phylo process\n";
		phyloprocess = new OneMatrixPhyloProcess(gamtree,codonmatrix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "register\n";

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		RootRegister(PriorLambda);
		RootRegister(PriorOmega);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		cerr << "model created\n";
		cerr << "initialise\n";
		Sample();

		Trace(cerr);
		MakeScheduler();
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~MGModel() {}

	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(omega,1),10,"alpha");
		scheduler.Register(new SimpleMove(omega,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(omega,0.01),10,"alpha");
		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"branch lengths");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new ProfileMove(stationary,1,2),10,"stat4");
		scheduler.Register(new SimpleMove(stationary,0.1),10,"stat");
		scheduler.Register(new SimpleMove(stationary,0.01),10,"stat");
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,1),4,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,0.1),4,"lengthrelrate");
	}


	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();
		total += stationary->GetLogProb();
		total += relrate->GetLogProb();
		total += omega->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}


	// Metropolis Hastings, successively called on each component of the model

	double OldFashionedMove(double tuning)	{
		for (int i=0; i<10; i++)	{
			for (int j=0; j<10; j++)	{
				lambda->Move(tuning);
				lambda->Move(0.1*tuning);
			}
			for (int j=0; j<10; j++)	{
				omega->Move(tuning);
				omega->Move(0.1*tuning);
			}
			for (int j=0; j<10; j++)	{
				relrate->Move(tuning/10);
			}

			for (int k=0; k<2; k++)	{
				stationary->Move(tuning,1);
				stationary->Move(tuning/10,1);
				stationary->Move(tuning/50,2);
				stationary->Move(0.1*tuning);
				stationary->Move(0.01*tuning);
				stationary->Move(0.001*tuning);
			}

			gamtree->Move(tuning);
			gamtree->Move(0.1 * tuning);
		}
		phyloprocess->Move(tuning);
		return 1;
	}

	// Draw a sample from the prior

	void drawSample()	{
		lambda->Sample();
		gamtree->Sample();
		stationary->Sample();
		relrate->Sample();
		omega->Sample();
		phyloprocess->Sample();
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetLength()	{
		return gamtree->GetTotalLength();
		// return gamtree->GetTotalLength() * matrix->GetRate();
	}


	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#t1\tt2\tt\tlogprior\tlnL\tomega\tlength\tlambda\ttstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << omega->val();
		os << '\t' << GetLength();
		os << '\t' << lambda->val();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *lambda<< '\n';
		os << *omega << '\n';
		os << *gamtree << '\n';
		os << *stationary << '\n';
		os << *relrate << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *omega;
		is >> *gamtree;
		is >> *stationary;
		is >> *relrate;
	}

};


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	MGModel* model = new MGModel(datafile,treefile);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Trace(os);
	}
}
