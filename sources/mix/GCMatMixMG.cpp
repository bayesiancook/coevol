
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugateOneMatrixPhyloProcess.h"
#include "ConjugateGammaTree.h"

#include "Move.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"

#include "BranchMatrixMixture.h"
#include "GCProcess.h"
#include "PrecisionNormalTreeProcess.h"

class OmegaMixtureRandomMatrixTree : public MixtureRandomMatrixTree<PosReal>	{

	public:

	OmegaMixtureRandomMatrixTree(Var<PosReal>* inalpha, Var<PosReal>* inbeta, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* innucmatrixtree, bool innormalise = false)	{
		alpha = inalpha;
		beta = inbeta;
		statespace = instatespace;
		nucmatrixtree = innucmatrixtree;
		normalise = innormalise;
		CreateRandomVariable();
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return nucmatrixtree->GetTree();}

	Rvar<PosReal>* CreateRandomVariable()	{
		rvar = new Gamma(alpha,beta);
		rvar->SetName("rvar");
		return rvar;
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		RandomSubMatrix* matrix = new RandomMGOmegaCodonSubMatrix(statespace, nucmatrixtree->GetBranchVal(link->GetBranch()), rvar, normalise);
		matrix->SetName("matrix");
		return matrix;
	}
	
	private:
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	bool normalise;
};


class OmegaBranchMatrixFiniteMixture : public BranchMatrixFiniteMixture<PosReal>	{

	public:

	OmegaBranchMatrixFiniteMixture(int insize, int incomponentnumber, Var<PosReal>* inalpha, Var<PosReal>* inbeta, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* innucmatrixtree, bool innormalise = false) :
			BranchMatrixFiniteMixture<PosReal>(innucmatrixtree->GetTree(), insize, incomponentnumber)	{

		alpha = inalpha;
		beta = inbeta;
		statespace = instatespace;
		nucmatrixtree = innucmatrixtree;
		normalise = innormalise;
		Create();
	}

	Tree* GetTree() {return nucmatrixtree->GetTree();}

	protected:

	MixtureRandomMatrixTree<PosReal>* CreateComponent(int k)	{
		return new OmegaMixtureRandomMatrixTree(alpha,beta,statespace,nucmatrixtree,normalise);
	}

	private:
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	bool normalise;
};


class GCMixMGModel : public ProbModel	{

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
	Const<Real>* Zero;

	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	ConjugateGammaTree* gamtree;
	
	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* gctree;

	Dirichlet* relrate;
	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;

	Dvar<PosReal>* PriorOmega;
	Exponential* alpha;
	Exponential* beta;

	int K;
	
	BranchMatrixMixturePhyloProcess<PosReal>* phyloprocess;
	
	public:

	OmegaBranchMatrixFiniteMixture* omega;

	GCMixMGModel(string datafile, string treefile, int inK)	{
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
		
		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		// a tree with all branch lengths iid from an exponential distribution of rate lambda 
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new ConjugateGammaTree(tree,mu,lambda);

		// another log normal process for the variations of gc
		PriorTheta = new Const<PosReal>(1);		
		theta = new Gamma(One,PriorTheta);
		gctree = new LogNormalTreeProcess(gamtree,theta);
		
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stattree = new GCStatTree(gctree);
		nucmatrixtree = new NucMatrixTree(relrate,stattree);

		PriorOmega = new Const<PosReal>(1);
		alpha = new Exponential(PriorOmega, Exponential::MEAN);
		beta = new Exponential(PriorOmega, Exponential::MEAN);

		// number of components
		K = inK;
		// mixture of K components
		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();
		omega = new OmegaBranchMatrixFiniteMixture(Nsite,K,alpha,beta,codonstatespace,nucmatrixtree);

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixMixturePhyloProcess<PosReal>(gamtree,omega,codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "register\n";

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorTheta);
		RootRegister(PriorLambda);
		RootRegister(PriorOmega);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(omega->GetWeightVector());
		Register();

		cerr << "model created\n";

		cerr << "initialise\n";
		Sample();
		Update();
		Trace(cerr);

		cerr << "scheduler\n";
		MakeScheduler();

		phyloprocess->SetName("phyloprocess");

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GCMixMGModel() {}

	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();
		total += theta->GetLogProb();
		total += gctree->GetLogProb();
		total += relrate->GetLogProb();
		total += omega->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}


	// Metropolis Hastings, successively called on each component of the model
	
	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(beta,1),10,"beta");
		scheduler.Register(new SimpleMove(beta,0.1),10,"beta");
		scheduler.Register(new SimpleMove(omega,1),1,"omega");
		scheduler.Register(new SimpleMove(omega,0.1),1,"omega");

		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");

		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),100,"theta");
		scheduler.Register(new SimpleMove(gctree,1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.01),30,"gc");

		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}

	/*
	double Move(double tuning)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	// Draw a sample from the prior

	void drawSample()	{
		lambda->Sample();
		gamtree->Sample();

		alpha->Sample();
		beta->Sample();
		omega->Sample();

		relrate->Sample();
		theta->Sample();
		theta->setval(10);
		gctree->Sample();
		gctree->Reset();

		phyloprocess->Sample();
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetLength()	{
		return gamtree->GetTotalLength(); 
		// return gamtree->GetTotalLength() * matrix->GetRate(); 
	}

	double GetMeanOmega()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*omega)[i]->val();
		}
		mean /= Nsite;
		return mean;
	}

	double GetVarOmega()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nsite; i++)	{
			double tmp = (*omega)[i]->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= Nsite;
		var /= Nsite;
		var -= mean * mean;
		return var;
	}

	void TraceOmega(ostream& os)	{
		os << *omega;
		os.flush();
	}
		

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#t1\tt2\tt\tlogprior\tlnL\tmeano\tstderro\talpha\tbeta\tlength\tlambda\ttstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanOmega();
		os << '\t' << sqrt(GetVarOmega());
		os << '\t' << alpha->val();
		os << '\t' << beta->val();
		os << '\t' << GetLength();
		os << '\t' << lambda->val();
		os << '\t' << stattree->GetMeanGCContent();
		os << '\t' << stattree->GetVarGCContent();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *alpha << '\n';
		os << *beta << '\n';
		os << *omega << '\n';
		os << *gamtree << '\n';
		os << *theta << '\n';
		os << *gctree << '\n';
		os << *relrate << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *omega;
		is >> *gamtree;
		is >> *theta;
		is >> *gctree;
		is >> *relrate;
	}

	void Monitor(ostream& os)	{
		scheduler.ToStream(os);
	}


};


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	int K = atoi(argv[3]);
	string name = argv[4];

	// MGModel* model = new MGModel(datafile,treefile);
	GCMixMGModel* model = new GCMixMGModel(datafile,treefile,K);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	ofstream ros((name + ".omega").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Move(0.1);
		model->Move(0.01);
		model->Trace(os);
		model->omega->ToStream(ros);
		ros << '\n';
		ros.flush();
		ofstream mon_os((name + ".monitor").c_str());
		model->Monitor(mon_os);
	}
}
