
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugateOneMatrixPhyloProcess.h"
#include "ConjugateGammaTree.h"
#include "Chrono.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "IID.h"


class SiteMGModel : public ProbModel	{

	public:
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

	int burnin;
	int size;

	double* siteomega;
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
	GammaIIDArray* omega;

	RandomSubMatrix** codonmatrix;

	SiteMatrixPhyloProcess* phyloprocess;

	Chrono ch1, ch2;
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	SiteMGModel(string datafile, string treefile, int inburnin)	{
		// fetch data from file

		burnin = inburnin;
		size = 0;

		nucdata = new FileSequenceAlignment(datafile);

		codondata = new CodonSequenceAlignment(nucdata, true);

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		cerr << Nsite << '\t' << Nstate << '\n';

		siteomega = new double[Nsite];
		for (int i=0; i<Nsite; i++)	{
			siteomega[i] = 0;
		}

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
		beta = new Exponential(PriorOmega, Exponential::MEAN);
		omega = new GammaIIDArray(Nsite,alpha,beta);


		// codonmatrix = new MGOmegaCodonRandomSubMatrix*[Nsite];
		codonmatrix = new RandomSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			codonmatrix[i] = new RandomMGOmegaCodonSubMatrix((CodonStateSpace*) codondata->GetStateSpace(),matrix,(*omega)[i]);
		}

		cerr << "create phylo process\n";
		phyloprocess = new SiteMatrixPhyloProcess(gamtree,codonmatrix,codondata);
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
		ch1.Reset();
		ch2.Reset();

		Trace(cerr);
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~SiteMGModel() {}

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

	double Move(double tuning)	{

		ch1.Start();
		for (int i=0; i<10; i++)	{
			for (int j=0; j<10; j++)	{
				lambda->Move(tuning);
				lambda->Move(0.1*tuning);
			}
			for (int j=0; j<10; j++)	{
				alpha->Move(tuning);
				beta->Move(tuning);

				alpha->Move(0.1*tuning);
				beta->Move(0.1*tuning);
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
			/*
			gamtree->Integrate();
			for (int j=0; j<10; j++)	{
				lambda->Move(tuning);
			}
			gamtree->Resample();
			*/
		}
		ch1.Stop();

		ch2.Start();
		phyloprocess->Move(tuning);
		ch2.Stop();
		return 1;
	}

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

	void PushOmega()	{
		for (int i=0; i<Nsite; i++)	{
			siteomega[i] += (*omega)[i]->val();
		}
		size++;
	}

	void WriteOmega(ostream& os)	{
		for (int i=0; i<Nsite; i++)	{
			os << siteomega[i] / size << '\t';
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << ch1.GetTime()/1000 << '\t' << ch2.GetTime()/1000 << '\t'<<  (ch1.GetTime() + ch2.GetTime()) / 1000 << '\t' << GetLogPrior() << '\t' << GetLogLikelihood();
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

	void MakeScheduler() {}

};

