
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "OneMatrixPhyloProcess.h"
#include "IID.h"
#include "Move.h"
#include "GTRSubMatrix.h"


class DirichletIIDArray : public IIDArray<Profile>	{

	public:

	DirichletIIDArray(int insize, Var<Profile>* incenter, Var<PosReal>* inconcentration) : IIDArray<Profile>(insize)	{
		center = incenter;
		concentration = inconcentration;
		Create();
	}

	Dirichlet* GetDirichletVal(int site)	{
		return dynamic_cast<Dirichlet*>(GetVal(site));
	}

	double GetMeanEntropy()	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += GetDirichletVal(i)->GetEntropy();
		}
		mean /= GetSize();
		return mean;
	}

	protected:

	Rvar<Profile>* CreateVal(int site)	{
		return new Dirichlet(center, concentration);
	}

	Var<Profile>* center;
	Var<PosReal>* concentration;
};


class MoveStatArray : public MCUpdate	{

	public:

	MoveStatArray(DirichletIIDArray* instatarray, double intuning, int inm) : statarray(instatarray), tuning(intuning), m(inm)	{}

	double Move(double tuning_modulator)	{
		double total = 0;
		for (int i=0; i<statarray->GetSize(); i++)	{
			total += statarray->GetDirichletVal(i)->Move(tuning_modulator * tuning,m);
		}
		return total / statarray->GetSize();
	}

	private:

	DirichletIIDArray* statarray;
	double tuning;
	int m;
};

class SiteStatGTRModel : public ProbModel	{

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

	// tree and branch lengths
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	Const<PosReal>* mu;
	GammaTree* gamtree;

	// nucleotide relative exchange rates
	IIDExp* relrate;

	Dirichlet* center;
	Const<PosReal>* PriorConcentration;
	Exponential* concentration;

	DirichletIIDArray* statarray;
	RandomSubMatrix** matrixarray;

	SiteMatrixPhyloProcess* phyloprocess;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	SiteStatGTRModel(string datafile, string treefile) {
		// fetch data from file

		data = new FileSequenceAlignment(datafile);

		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)

		cerr << Nsite << '\t' << Nstate << '\n';

		taxonset = data->GetTaxonSet();

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
		relrate = new IIDExp(Nstate*(Nstate-1)/2);

		PriorConcentration = new Const<PosReal>(1);
		concentration = new Exponential(PriorConcentration,Exponential::MEAN);
		center = new Dirichlet(Nstate);

		statarray = new DirichletIIDArray(Nsite,center,concentration);
		matrixarray = new RandomSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			matrixarray[i] = new GTRRandomSubMatrix(relrate,statarray->GetVal(i));
		}

		cerr << "create phylo process\n";
		phyloprocess = new SiteMatrixPhyloProcess(gamtree,matrixarray,data);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "register\n";

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		RootRegister(PriorLambda);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(PriorConcentration);
		RootRegister(center);
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

	~SiteStatGTRModel() {}

	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),4,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),4,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),4,"branch lengths");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new SimpleMove(center,0.1),10,"center");
		scheduler.Register(new SimpleMove(concentration,0.1),10,"concentration");
		scheduler.Register(new MoveStatArray(statarray,1,1),4,"stat");
		scheduler.Register(new MoveStatArray(statarray,0.3,2),4,"stat");
		scheduler.Register(new MoveStatArray(statarray,0.1,5),4,"stat");
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
		total += relrate->GetLogProb();
		total += center->GetLogProb();
		total += concentration->GetLogProb();
		total += statarray->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	// Draw a sample from the prior

	void drawSample()	{

		lambda->Sample();
		gamtree->Sample();

		center->Sample();
		concentration->Sample();

		// here I cheat: I make the hyperprior flat
		center->setuniform();
		concentration->setval(1);

		relrate->Sample();
		statarray->Sample();
		phyloprocess->Sample();
	}

	// various summary statistics
	// used to check mcmc convergence

	double GetLength()	{
		return gamtree->GetTotalLength();
	}

	double GetMeanEntropy()	{
		return statarray->GetMeanEntropy();
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "logprior\tlnL\tlength\ttstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();
		os << '\t' << GetMeanEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	// to be implemented
	void ToStream(ostream& os)	{
	}

	void FromStream(istream& is)	{
	}

};


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	SiteStatGTRModel* model = new SiteStatGTRModel(datafile,treefile);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Trace(os);
		ofstream mos((name + ".monitor").c_str());
		model->Monitor(mos);
	}
}
