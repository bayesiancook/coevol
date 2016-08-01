
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"

#include "OneMatrixPhyloProcess.h"

#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"

#include "ConjugateGammaTree.h"

#include "GTRSubMatrix.h"

#include "Chrono.h"
#include "Jeffreys.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "IID.h"


class ConjugateSiteOmegaMove : public MCUpdate, public Mnode {

	GammaIIDArray* omega;
	MGCodonSiteMatrixPathConjugateProcess* pathconjprocess;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

	public:

	int GetSize()	{
		return omega->GetSize();
	}

	ConjugateSiteOmegaMove(GammaIIDArray* inomega, MGCodonSiteMatrixPathConjugateProcess* inpathconjprocess)	{
		omega = inomega;
		pathconjprocess = inpathconjprocess;
		alpha = omega->GetAlpha();
		beta = omega->GetBeta();
		omega->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(false);
		pathconjprocess->ComputeTotSuffStat();
		ResampleOmega();
		Update();
		return 1;
	}

	void ResampleOmega()	{
		for (int i=0; i<GetSize(); i++)	{
			double currentomega = omega->GetVal(i)->val();
			double tmp = Random::Gamma(alpha->val() + pathconjprocess->GetTotNonSynCount(i), beta->val() + pathconjprocess->GetTotNonSynBeta(i) / currentomega);
			omega->GetVal(i)->setval(tmp);
		}
	}

};

class SiteOmegaModel : public ProbModel	{

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

	// tree and branch lengths
	Const<PosReal>* One;
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	ConjugateGammaTree* gamtree;

	// nucleotide relative exchange rates
	Dirichlet* relrate;
	Dirichlet* stationary;

	// nucleotide substitution matrix is relrate * stationary
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	Jeffreys* alpha;
	Jeffreys* beta;
	GammaIIDArray* omega;

	RandomSubMatrix** codonmatrix;

	PathConjugateProcess* pathconjprocess;
	MGCodonSiteMatrixPathConjugateProcess* mgpathconjprocess;
	PhyloProcess* phyloprocess;

	double mappingfreq;
	bool clamptree;
	int conjpath;
	int conjom;
	bool normalise;
	int nrep;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	SiteOmegaModel(string datafile, string treefile, int inconjpath, int inconjom, double inmappingfreq, bool inclamptree, bool innormalise, int innrep, bool sample, GeneticCodeType type=Universal)	{


		conjpath = inconjpath;
		conjom = inconjom;
		mappingfreq = inmappingfreq;
		if (mappingfreq == -1)	{
			/*
			if (omegaratiotree && (rawNstate == Nnuc))	{
				mappingfreq = 0.2;
			}
			else	{
				mappingfreq = 1;
			}
			*/
			mappingfreq = 1;
		}
		clamptree = inclamptree;
		normalise = innormalise;
		nrep = innrep;

		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true,type);

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

		One = new Const<PosReal>(1.0);
		PriorLambda = new Const<PosReal>(10);

		// a tree with all branch lengths iid from an exponential distribution of rate lambda
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new ConjugateGammaTree(tree,One,lambda);

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);

		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		double minjeff = 0.01;
		double maxjeff = 100;
		alpha = new Jeffreys(minjeff,maxjeff,One);
		beta = new Jeffreys(minjeff,maxjeff,One);
		omega = new GammaIIDArray(Nsite,alpha,beta);

		codonmatrix = new RandomSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			codonmatrix[i] = new RandomMGOmegaCodonSubMatrix((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,(*omega)[i]);
		}

		cerr << "create phylo process\n";
		mgpathconjprocess = 0;
		if (conjpath)	{
			cerr << "codon conjugate\n";
			mgpathconjprocess = new MGCodonSiteMatrixPathConjugateProcess(gamtree, codonmatrix, codondata);
			pathconjprocess = mgpathconjprocess;
			phyloprocess = new PathConjugatePhyloProcess(pathconjprocess);
			
		}
		/*
		if (conjpath)	{
			cerr << "conjugate\n";
			pathconjprocess = new SiteMatrixPathConjugateProcess(gamtree, codonmatrix, codondata);
			phyloprocess = new PathConjugatePhyloProcess(pathconjprocess);
			
		}
		*/
		else	{
			pathconjprocess = 0;
			phyloprocess = new SiteMatrixPhyloProcess(gamtree,codonmatrix,codondata);
		}
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorLambda);
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		cerr << "make scheduler\n";
		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 5 : 1;
		}
		MakeScheduler();
		if (sample)	{
			cerr << "update\n";
			Update();
		}
		Trace(cerr);
		cerr << "model created\n";
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~SiteOmegaModel() {}

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
		total += alpha->GetLogProb();
		total += beta->GetLogProb();
		total += omega->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{

		double ret = 0;
		if (pathconjprocess)	{
			ret = pathconjprocess->GetLogProb();
		}
		ret = phyloprocess->GetLogProb();
		return ret;
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

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
		os << "#logprior\tlnL\tmeano\tstderro\talpha\tbeta\tlength\tlambda\ttstatent\trrent\n";
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
		os << '\t' << stationary->val().GetEntropy();
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
		os << *stationary << '\n';
		os << *relrate << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *alpha;
		is >> *beta;
		is >> *omega;
		is >> *gamtree;
		is >> *stationary;
		is >> *relrate;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjprocess,mappingfreq),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		for (int i=0; i<nrep; i++)	{

			if (! clamptree)	{
				scheduler.Register(new SimpleMove(gamtree,1),10,"gamtree");
				scheduler.Register(new SimpleMove(gamtree,0.1),10,"gamtree");
				scheduler.Register(new SimpleMove(gamtree,0.01),10,"gamtree");
			}

			scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.01),10,"lambda");

			if (conjom)	{
				if (! mgpathconjprocess)	{
					cerr << "error in conjugate site omega move: pathconjprocess not allocated\n";
					exit(1);
				}
				scheduler.Register(new ConjugateSiteOmegaMove(omega,mgpathconjprocess),1,"conj om");
			}
			else	{
				scheduler.Register(new SimpleMove(omega,10),10,"om");
				scheduler.Register(new SimpleMove(omega,1),10,"om");
				scheduler.Register(new SimpleMove(omega,0.1),10,"om");
			}

			scheduler.Register(new SimpleMove(alpha,10),100,"omalpha");
			scheduler.Register(new SimpleMove(alpha,1),100,"omalpha");
			scheduler.Register(new SimpleMove(alpha,0.1),100,"omalpha");

			scheduler.Register(new SimpleMove(beta,10),100,"ombeta");
			scheduler.Register(new SimpleMove(beta,1),100,"ombeta");
			scheduler.Register(new SimpleMove(beta,0.1),100,"ombeta");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}
	}
};

