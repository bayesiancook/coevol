
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugateOneMatrixPhyloProcess.h"
#include "ConjugateGammaTree.h"

#include "Move.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"

#include "MatrixMixture.h"

class OmegaMixtureRandomMatrix : public MixtureRandomMatrix<PosReal>	{

	public:

	OmegaMixtureRandomMatrix(Var<PosReal>* inalpha, Var<PosReal>* inbeta, CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, bool innormalise = false)	{
		alpha = inalpha;
		beta = inbeta;
		statespace = instatespace;
		nucmatrix = inmatrix;
		normalise = innormalise;
		CreateRandomVariable();
		CreateRandomSubMatrix();
	}

	Rvar<PosReal>* CreateRandomVariable()	{
		rvar = new Gamma(alpha,beta);
		rvar->SetName("rvar");
		return rvar;
	}

	RandomSubMatrix* CreateRandomSubMatrix()	{
		matrix = new RandomMGOmegaCodonSubMatrix(statespace, nucmatrix, rvar, normalise);
		matrix->SetName("matrix");
		return matrix;
	}
	
	private:
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
	bool normalise;
};


class OmegaMatrixInfiniteMixture : public MatrixInfiniteMixture<PosReal>	{

	public:

	OmegaMatrixInfiniteMixture(int insize, int incomponentnumber, Var<PosReal>* inalpha, Var<PosReal>* inbeta, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, bool innormalise = false) :
			MatrixInfiniteMixture<PosReal>(insize, incomponentnumber)	{

		alpha = inalpha;
		beta = inbeta;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		normalise = innormalise;
		Create();
	}

	protected:

	MixtureRandomMatrix<PosReal>* CreateComponent(int k)	{
		OmegaMixtureRandomMatrix* tmp = new OmegaMixtureRandomMatrix(alpha,beta,statespace,nucmatrix,normalise);
		tmp->GetRandomVariable()->Corrupt(true);
		tmp->GetRandomVariable()->Update();
		return tmp;
	}

	private:
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
	bool normalise;
};


class MixMGModel : public ProbModel	{

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

	double* siteomega;
	int size;

	// ---------
	// the random variables of the model
	// ---------

	// tree and branch lengths
	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	ConjugateGammaTree* gamtree;
	
	// nucleotide relative exchange rates
	IIDExp* relrate;
	Dirichlet* stationary;

	// nucleotide substitution matrix is relrate * stationary
	GTRRandomSubMatrix* matrix;
	
	Dvar<PosReal>* PriorOmega;
	Exponential* alpha;
	Exponential* beta;

	int K;
	
	MatrixInfiniteMixturePhyloProcess<PosReal>* phyloprocess;
	
	public:

	OmegaMatrixInfiniteMixture* omega;
	MixMGModel(string datafile, string treefile, int inK)	{
		// fetch data from file

		size = 0;
	
		nucdata = new FileSequenceAlignment(datafile);

		codondata = new CodonSequenceAlignment(nucdata, true);

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		siteomega = new double[Nsite];
		for (int i=0; i<Nsite; i++)	{
			siteomega[i] = 0;
		}

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
		gamtree = new ConjugateGammaTree(tree,mu,lambda);
		relrate = new IIDExp(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);

		matrix = new GTRRandomSubMatrix(relrate,stationary);

		PriorOmega = new Const<PosReal>(1);
		alpha = new Exponential(PriorOmega, Exponential::MEAN);
		alpha->setval(2.0);
		beta = new Exponential(PriorOmega, Exponential::MEAN);

		// initial number of components 
		K = inK;
		// infinite mixture
		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();
		omega = new OmegaMatrixInfiniteMixture(Nsite,K,alpha,beta,codonstatespace,matrix);

		cerr << "create phylo process\n";
		phyloprocess = new MatrixInfiniteMixturePhyloProcess<PosReal>(gamtree,omega,codondata);
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
		// RootRegister(omega->GetWeightVector());
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

	~MixMGModel() {}

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
		scheduler.Register(new SimpleMove(relrate,1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new ProfileMove(stationary,0.3,1),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.1,2),10,"stat4");
		scheduler.Register(new SimpleMove(stationary,0.1),10,"stat");
		scheduler.Register(new SimpleMove(stationary,0.01),10,"stat");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,1),10,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,0.1),10,"lengthrelrate");
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

		stationary->Sample();
		relrate->Sample();

		phyloprocess->Sample();
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetLength()	{
		return gamtree->GetTotalLength(); 
		// return gamtree->GetTotalLength() * matrix->GetRate(); 
	}

	double GetPosFrac()	{

		double count = 0;
		for (int i=0; i<Nsite; i++)	{
			if ((*omega)[i]->val() > 1)	{
				count++;
			}
		}
		return count / Nsite;
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

	/*
	void TraceOmega(ostream& os)	{
		for (int i=0; i<Nsite; i++)	{
			double tmp = (*omega)[i]->val();
			os << tmp << '\t';
		}
		os << '\n';
		os.flush();
	}
	*/	

	/*
	void TraceOmega(ostream& os)	{

		size++;
		double count = 0;
		for (int i=0; i<Nsite; i++)	{
			siteomega[i] += (*omega)[i]->val();
			if (siteomega[i] / size > 1)	{
				count++;
			}
		}
		os << count / Nsite;
		for (int i=0; i<Nsite; i++)	{
			os << '\t' << siteomega[i] / size;
		}
		os << '\n';
		os.flush();
	}
	*/

	void TraceOmega(ostream& os)	{

		size++;
		double count = 0;
		for (int i=0; i<Nsite; i++)	{
			if ((*omega)[i]->val() > 1)	{
				siteomega[i] ++;
				count++;
			}
		}
		os << count / Nsite;
		for (int i=0; i<Nsite; i++)	{
			os << '\t' << siteomega[i] / size;
		}
		os << '\n';
		os.flush();
	}
		

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tK\tposfrac\tmeano\tstderro\talpha\tbeta\tlength\tlambda\ttstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << omega->GetComponentNumber();
		os << '\t' << GetPosFrac();
		os << '\t' << GetMeanOmega();
		os << '\t' << sqrt(GetVarOmega());
		os << '\t' << alpha->val();
		os << '\t' << beta->val();
		os << '\t' << GetLength();
		os << '\t' << lambda->val();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *alpha << '\n';
		os << *beta << '\n';
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
	int K = atoi(argv[3]);
	int burnin = atoi(argv[4]);
	string name = argv[5];

	// MGModel* model = new MGModel(datafile,treefile);
	MixMGModel* model = new MixMGModel(datafile,treefile,K);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	ofstream ros((name + ".omega").c_str());
	model->TraceHeader(os);

	int size = 0;
	while (1)	{
		size++;
		model->Move(1);
		model->Move(0.1);
		// model->Move(0.01);
		model->Trace(os);
		if (size > burnin)	{
			model->TraceOmega(ros);
			ros.flush();
		}
		ofstream mon_os((name + ".monitor").c_str());
		ofstream dmon_os((name + ".detailedmonitor").c_str());
		model->Monitor(mon_os, dmon_os);
	}
}
