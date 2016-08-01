
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "IID.h"
#include "ConjugateLogNormalTreeProcess.h"
#include "Chronogram.h"

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

	Const<PosReal>* One;

	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	// autocorrelated process
	Dvar<PosReal>* PriorSigma;
	ConjugateGamma* sigma;
	ConjugateLogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
	IIDExp* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrix* matrix;

	// rates across sites
	Dvar<PosReal>* PriorVarRate;
	Gamma* alpha;
	Gamma* beta;
	GammaIIDArray* rate;

	// phylo process
	OneMatrixRASPhyloProcess* phyloprocess;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRLogNormalModel(string datafile, string treefile)	{
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

		One = new Const<PosReal>(1);

		// a tree with all branch lengths iid from an exponential distribution of mean meanlength
		// meanlength is itself endowed with an exponential prior of mean 0.1
		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu);
		mu->ClampAt(1);
		chronogram = new Chronogram(tree,mu);

		// a log normal process on that tree
		PriorSigma = new Const<PosReal>(1);
		sigma = new ConjugateGamma(One,PriorSigma);
		lognormaltree = new ConjugateLogNormalTreeProcess(chronogram,sigma);

		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1
		PriorVarRate = new Const<PosReal>(1);
		alpha = new Gamma(One,PriorVarRate);
		beta = new Gamma(One,PriorVarRate);
		rate = new GammaIIDArray(Nsite,alpha,beta);

		// substitution matrix
		relrate = new IIDExp(Nstate*(Nstate-1)/2);
		stationary = new Dirichlet(Nstate);
		matrix = new GTRRandomSubMatrix(relrate,stationary);

		// a phylogenetic process
		phyloprocess = new OneMatrixRASPhyloProcess(lognormaltree,rate,matrix,data);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(PriorVarRate);
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		cerr << "sample model\n";
		Sample();

		cerr << "update\n";
		Update();

		TraceHeader(cerr);
		Trace(cerr);

		cerr << "scheduler\n";
		MakeScheduler();
		cerr << "model created\n";

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GTRLogNormalModel() {}

	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;
		total += mu->GetLogProb();
		total += chronogram->GetLogProb();
		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();
		total += alpha->GetLogProb();
		total += beta->GetLogProb();
		total += rate->GetLogProb();
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

		/*
		scheduler.Register(new SemiConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,1,10,1,0),1,"semiconjlogn");
		scheduler.Register(new SemiConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,0.1,10,1,0),1,"semiconjlogn");
		scheduler.Register(new ConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,1,0),1,"conjlogn");
		scheduler.Register(new ConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,0.1,0),1,"conjlogn");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");
		*/

		scheduler.Register(new ConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,1,10),1,"conjlogn");
		scheduler.Register(new ConjugateMove<ConjugateGamma,ConjugateLogNormalTreeProcess>(sigma,lognormaltree,0.1,10),1,"conjlogn");

		/*
		scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");
		*/

		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(beta,1),10,"alpha");
		scheduler.Register(new SimpleMove(beta,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(rate,1),10,"rate");
		scheduler.Register(new SimpleMove(rate,0.1),10,"rate");
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
		// chronogram->Sample();
		sigma->Sample();
		lognormaltree->Sample();
		alpha->Sample();
		beta->Sample();
		rate->Sample();
		stationary->Sample();
		relrate->Sample();
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
		return lognormaltree->GetTotalLength() * matrix->GetRate() * GetMeanRate();
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tsigma\tmeanrho\tvarrho\trootlog\talpha\tmeanrate\tvarrate\tstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << sigma->val();
		os << '\t' << GetMeanRho() << '\t' << GetVarRho();
		os << '\t' << *(lognormaltree->GetNodeVal(lognormaltree->GetRoot()->GetNode()));
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

	GTRLogNormalModel* model = new GTRLogNormalModel(datafile,treefile);

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
