
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "Mixture.h"
// #include "Mixturepointers.h"
#include "OneMatrixPhyloProcess.h"
#include "Move.h"

class GamMixModel: public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

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
	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	GammaTree* gamtree;
	
	// relative exchange rates of the matrix
	IIDExp* relrate;

	// equilibrium frequencies of the matrix
	Dirichlet* stationary;

	// substitution matrix is relrate * stationary
	GTRRandomSubMatrix* matrix;
	
	// rates across sites
	Dvar<PosReal>* PriorVarRate;
	Exponential* alpha;
	OneMatrixRASPhyloProcess* phyloprocess;
	
	public:

	int K;
	GammaMixture* rate;
	// constructor
	// this is where the entire graph structure of the model is created

	GamMixModel(string datafile, string treefile) {
		cerr << "data and tree file\n";
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

		K = 10;
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
		/*
		gamtree->SetBranchLengths();
		gamtree->Clamp();
		*/

		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1
		PriorVarRate = new Const<PosReal>(PosReal(1));		
		alpha = new Exponential(PriorVarRate,Exponential::MEAN);

		cerr << "rate mixture\n";
		rate = new GammaMixture(Nsite,K, alpha,alpha);
		cerr << "ok\n";
		//
		// a general time reversible (GTR) substitution matrix
		// for any pair of states l,m (distinct)
		//
		// Q_lm = rho_lm pi_m
		//
		// where rho_lm are the relative rates
		// and pi_m is the vector of equilibrium frequencies of the substitution process
		//
		// the proces is reversible, thus rho_lm = rho_ml, and there are Nstate (Nstate-1) / 2 distinct relative rates
		// all iid from an exponential of mean 1
		// the stationary probabilities (pi) are from a uniform distribution (Dirichlet of weights 1,1,1,1,...)
		relrate = new IIDExp(Nstate*(Nstate-1)/2);
		stationary = new Dirichlet(Nstate);
		matrix = new GTRRandomSubMatrix(relrate,stationary);

		/*
		alpha->ClampAt(1);
		for (int i=0; i<Nsite; i++)	{
			(*rate)[i]->ClampAt(1);
		}

		relrate->setall(1);
		relrate->Clamp();

		double* empfreq = new double[Nstate];
		data->GetEmpiricalFreq(empfreq);
		stationary->setarray(empfreq);
		stationary->Clamp();
		*/

		// a phylogenetic process
		// this is where everything comes together
		// the process is of type GTRRAS
		// which means that it gives the same GTR matrix for every site and every branch
		// but gives each site its own rate from the IID vector "rate"
		//
		cerr << "create phylo process\n";
		phyloprocess = new OneMatrixRASPhyloProcess(gamtree,rate,matrix,data);
		// a phyloprocess needs to "unfold"
		// (recursive construction of all branch/site substitution processes)
		cerr << "unfold\n";
		cerr << "one matrix model\n";
		phyloprocess->Unfold();
		cerr << "register\n";

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		RootRegister(PriorVarRate);
		RootRegister(PriorLambda);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(rate->GetWeightVector());
		Register();
	
		cerr << "initialise\n";
		Sample();
		cerr << "ok\n";
		Trace(cerr);

		cerr << "scheduler\n";
		MakeScheduler();

		cerr << "model created\n";
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment


	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += alpha->GetLogProb();
		total += rate->GetLogProb();
		total += gamtree->GetLogProb();
		total += stationary->GetLogProb();
		total += relrate->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(rate,1),10,"rate");
		scheduler.Register(new SimpleMove(rate,0.1),10,"rate");
		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"branch lengths");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new SimpleMove(stationary,1),10,"stat");
		scheduler.Register(new SimpleMove(stationary,0.1),10,"stat");
		scheduler.Register(new SimpleMove(stationary,0.01),10,"stat");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,1),3,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,1),3,"lengthrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,0.1),3,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,0.1),3,"lengthrate");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}

	// Draw a sample from the prior

	void drawSample()	{
		alpha->Sample();
		cerr << "rate\n";
		rate->Sample();
		for (int i=0; i<Nsite; i++)	{
			cerr << (*rate)[i]->val() << '\t';
		}
		cerr << '\n';
		cerr << "ok\n";
		lambda->Sample();
		cerr << "gamtree\n";
		gamtree->Sample();
		stationary->Sample();
		relrate->Sample();
		cerr << "phyloprocess sample\n";
		phyloprocess->Sample();
		cerr << "ok\n";
	}

	// various summary statistics
	// used to check mcmc convergence

	double GetMeanRate()	{
		double total = 0;
		for (int i=0; i<Nsite; i++)	{
			total += (*rate)[i]->val();
		}
		total /= Nsite;
		return total;
	}

	double GetVarRate()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nsite; i++)	{
			double tmp = (*rate)[i]->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= Nsite;
		var /= Nsite;
		var -= mean * mean;
		return var;
	}

	double GetLength()	{
		return gamtree->GetTotalLength() * matrix->GetRate() * GetMeanRate();
	}


	void Monitor(ostream& os)	{
		scheduler.ToStream(os);
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tlambda\talpha\tmeanrate\tvarrate\tstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << lambda->val();
		os << '\t' << alpha->val();
		os << '\t' << GetMeanRate();
		os << '\t' << GetVarRate();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
	}

	void FromStream(istream& is)	{
	}
};


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	// MGModel* model = new MGModel(datafile,treefile);
	GamMixModel* model = new GamMixModel(datafile,treefile);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	ofstream ros((name + ".rate").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Move(0.1);
		model->Move(0.01);
		model->Trace(os);
		model->rate->ToStream(ros);
		ros << '\n';
		ros.flush();
		ofstream mon_os((name + ".monitor").c_str());
		model->Monitor(mon_os);
	}
}
