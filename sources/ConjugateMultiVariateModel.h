#include "ProbModel.h"
#include "OneMatrixPhyloProcess.h"
#include "IID.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "Move.h"



class ConjugateExpoIntegrale : public Dvar<PosReal>{

	private:

		ConjugateMultiNormal* up;
		ConjugateMultiNormal* down;
		Dvar<PosReal>* time;
		int index;

	public:

	ConjugateExpoIntegrale(ConjugateMultiNormal* inup, ConjugateMultiNormal* indown, Dvar<PosReal>* intime, int inindex ){
		up =inup;
		down = indown;
		time = intime;
		index = inindex;
		Register(up);
		Register(down);
		Register(time);
	}

	~ConjugateExpoIntegrale(){};


	void specialUpdate(){
		setval( (exp(up->val()[index]) + exp(down->val()[index]))/2 * time->val()); //look
	}
};

class ConjugateBranchExpoIntegrale : public BranchValPtrTree<Dvar<PosReal> >{


	protected:


	ConjugateMultiVariateTreeProcess* mvtree;
	int index;


	public:

	ConjugateBranchExpoIntegrale(ConjugateMultiVariateTreeProcess* inmvtree, int inindex){
		index = inindex;
		mvtree = inmvtree;
		RecursiveCreate(GetRoot());
	}

	~ConjugateBranchExpoIntegrale(){
		RecursiveDelete(GetRoot());
	}


	Tree* GetTree(){
		return mvtree->GetTree();
	}

	Chronogram* GetLengthTree(){
		return mvtree->GetLengthTree();
	}

	ConjugateExpoIntegrale* CreateBranchVal(const Link* link){
		return new ConjugateExpoIntegrale(mvtree->GetMultiNormal(link), mvtree->GetMultiNormal(link->Out()), GetLengthTree()->GetBranchVal(link->GetBranch()), index);
	}
};


class SpecialMove : public MCUpdate {

	ConjugateInverseWishart* sigma;
	ConjugateMultiVariateTreeProcess* tree;
	RvarVec* priorOnSigmaZero;
	int nmove;
	double tuning;

	public:

	SpecialMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, RvarVec* inpriorOnSigmaZero, int innmove){
		sigma = insigma;
		tree = intree;
		nmove =innmove;
		priorOnSigmaZero = inpriorOnSigmaZero;
	}

	double Move(double tuning){
		sigma->ActivateSufficientStatistic();
		double r = 0;
		double p = 0;
		for(int i=0; i<nmove; i++){
			r += tree->Move(tuning);
			for(int j=0; j<10; j++){
				p += priorOnSigmaZero->Move(tuning);
			}
		}
		sigma->InactivateSufficientStatistic();
		return (r+p)/(11*nmove);
	}
};

class ConjugateMultiVariateModel : public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	MCScheduler scheduler;

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
	TaxonSet* taxonset;
	DataMatcher* otherdata;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 10 for amino-acids. 61 for codons)
	int Nstate;
	// Size of the multivariate vector
	int dim;


	// ---------
	// the random variables of the model
	// ---------

	// Ultrametric tree
	LengthTree* lengthtree;

	//Covariance prior and matrix:
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* sigma;

	// NodeMultiVariateTree
	ConjugateMultiVariateTreeProcess* multiVariateTreeProcess;


	// relative exchange rates of the matrix
	IIDExp* relrate;
	// equilibrium frequencies of the matrix
	Dirichlet* stationary;
	// substitution matrix is relrate * stationary
	GTRRandomSubMatrix* matrix;



	// usefull...
	Dvar<PosReal>* FixedUnity;


	OneMatrixPhyloProcess* phyloprocess;

	public:
	// constructor
	// this is where the entire graph structure of the model is created

	ConjugateMultiVariateModel(string datafile, string treefile, string otherdatafile) : scheduler(this)	{
		// fetch data from file
		data = new FileSequenceAlignment(datafile);


		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)

		taxonset = data->GetTaxonSet();


 		otherdata = new DataMatcher(otherdatafile);


		dim = otherdata->GetNdata() + 1;

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);


		cout << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------


		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1

		// Root of MCMC !!
		FixedUnity = new Const<PosReal>(1);


		// Create an ultrametric tree with branch lenght mu*time
		chronogram =  new Chronogram(tree, FixedUnity);

		Rvar<PosReal>** temp1  = new Rvar<PosReal>*[dim];
		for(int i=0; i<dim; i++){
			temp1[i] = new Exponential(FixedUnity, Exponential::MEAN);
			temp1[i]->SetName("priorOnS0");
// 			temp1[i]->setval(1);
//  			temp1[i]->Clamp();
		}
		priorOnSigmaZero = new RvarVec(temp1, dim);

		sigmaZero = new SigmaZero(priorOnSigmaZero);


		// Creation of the Wishart Matrix //loook P
		sigma = new ConjugateInverseWishart(sigmaZero, dim+1);
// 		sigma->SetAndClampAtSigmaZero();



		// Create an NodeMultiVariateTree
		multiVariateTreeProcess = new ConjugateMultiVariateTreeProcess(sigma, chronogram);

		ConjugateBranchExpoIntegrale* ratetree = new ConjugateBranchExpoIntegrale(multiVariateTreeProcess, 0);

		//Set the values in the file to the values of the process:
		otherdata->SetPosData(0);
		multiVariateTreeProcess->SetAndClamp(otherdata, 1);

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


		relrate->setall(1);
		relrate->Clamp();
		stationary->setuniform();
		stationary->Clamp();


		// a phylogenetic process
		// this is where everything comes together
		// the process is of type GTRRAS
		// which means that it gives the same GTR matrix for every site and every branch
		// but gives each site its own rate from the IID vector "rate"
		//
		phyloprocess = new OneMatrixPhyloProcess(ratetree,matrix,data);

		// phyloprocess = new GTRRASPhyloProcess(expotree,relrate,stationary,rate,data);
		cout << "unfold\n";
		// a phyloprocess needs to "unfold"
		// (recursive construction of all branch/site substitution processes)
		phyloprocess->Unfold();

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		cout << "register\n";
		RootRegister(FixedUnity);
		RootRegister(relrate);
		RootRegister(stationary);

		Register();

		cout << "model created\n";

		cout << "initialise -> ";
		Sample();
		cout << "ok\n";
		Trace(cout);

		cout << "update -> ";
		Update();
		cout << "ok\n";
		cout << "scheduler ->";
		MakeScheduler();
		cout << "ok\n";

		multiVariateTreeProcess->SetName("multiVariateTreeProcess");
		FixedUnity->SetName("FixedUnity");
		relrate->SetName("relrate");
		stationary->SetName("stationary");
		sigmaZero->SetName("sigmaZero");
		sigma->SetName("sigma");
		chronogram->SetName("chronogram");
		matrix->SetName("matrix");

		cout << "start sigma :" << * sigma << "\n";

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~ConjugateMultiVariateModel() {}

	// Update is called at the beginning
	// it samples a configuration from the prior
	// and updates all nodes (logprobs, backups, etc)

	double GetLength(){
		return chronogram->GetTotalLength();
	}



	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += chronogram->GetLogProb();
		total += priorOnSigmaZero->GetLogProb();
		total += sigma->GetLogProb();
		total += multiVariateTreeProcess->GetLogProb();
		//total += relrate->GetLogProb();
		//total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	//2 methods for Read() of Sample
	ConjugateInverseWishart* GetSigma(){
		return sigma;
	}

	ConjugateMultiVariateTreeProcess* GetMultiVariateTreeProcess(){
		return multiVariateTreeProcess;
	}



	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		//scheduler.Register(new SimpleMove(chronogram,0.01),3,"chronogram");
		scheduler.Register(new SpecialMove(sigma, multiVariateTreeProcess, priorOnSigmaZero, 10),1,"conjugate sigma");
		//scheduler.Register(new SimpleMove(stationary,0.1),10,"stationary");
		//scheduler.Register(new SimpleMove(relrate,1),10,"relrate");
	}


	// Metropolis Hastings, successively called on each component of the model

	double Move(double tuning)	{
		scheduler.Cycle(1, 1,false,false);
		return 1;
	}

	// Draw a sample from the prior

	void drawSample()	{
		stationary->Sample();
		relrate->Sample();
		priorOnSigmaZero->Sample();
		sigma->Sample();
		multiVariateTreeProcess->Sample();
		phyloprocess->Sample();
	}



	void Monitor(ostream& os)	{
		scheduler.ToStream(os);
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlikelihood\tsigma\tMVTP\tcov\tvar_taux\tvar_long\tdet\trate-root\tlong_root\tprior_rate\tprior_long";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior();
		os << '\t' << GetLogLikelihood();
		os << '\t' << sigma->GetLogProb();
		os << '\t' << multiVariateTreeProcess->GetLogProb();
		os << '\t' << sigma->GetMatrix()[0][1]; //5
		os << '\t' << sigma->GetMatrix()[0][0];
		os << '\t' << sigma->GetMatrix()[1][1];
		os << '\t' << sigma->GetDeterminant();
		os << '\t' << multiVariateTreeProcess->GetMultiNormal(multiVariateTreeProcess->GetRoot())->val()[0]; //9
		os << '\t' << multiVariateTreeProcess->GetMultiNormal(multiVariateTreeProcess->GetRoot())->val()[1];
		os << '\t' << priorOnSigmaZero->val(0)->val();
		os << '\t' << priorOnSigmaZero->val(1)->val();
		os << '\n';
	}

	void ToStream(ostream& os)	{
		os << *stationary << '\n';
		os << *relrate << '\n';
		os << *chronogram << '\n';
		os << *multiVariateTreeProcess << '\n';
		os << *sigma << '\n';
		for(int i=0; i<dim; i++){
			os << *priorOnSigmaZero->val(i) << '\n';
		}

	}
	void FromStream(istream& is)	{
		is >> *stationary;
		is >> *relrate;
		is >> *chronogram;
		is >> *multiVariateTreeProcess;
		is >> *sigma;
		for(int i=0; i<dim; i++){
			is >> *(priorOnSigmaZero->val(i));
		}
	}

	// various accessory functions to print tree with branch lengths
	// and make stochastic mapping of substitution events along the tree
	void PrintTree(ostream& os){
		os << 	multiVariateTreeProcess->Construct() + "\n";
	}

};
