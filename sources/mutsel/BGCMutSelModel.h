
// a chronogram
//
// 2 independent log normal processes
// - log mutation rate (per unit of time)
// - log popsize (meant to be clamped at two leaf nodes)
// log(omega) = beta + alpha * log(popsize)

#ifndef POPSIZEMODEL_H
#define POPSIZEMODEL_H


#include "Random.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "Chronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ScatterMatrix.h"
#include "MGMixAAMutSelCodonSubMatrix.h"
#include "HKYRandomSubMatrix.h"

class InstantStat : public Dvar<Profile>	{

	public:

	InstantStat(Var<Real>* ingc)	{
		setval(Profile(Nnuc));
		bkvalue = Profile(Nnuc);
		gc = ingc;
		Register(gc);
	}

	void specialUpdate()	{
		double f = exp(gc->val()) / (1 + exp(gc->val()));
		(*this)[0] = (*this)[3] = 0.5 * (1 - f);
		(*this)[1] = (*this)[2] = 0.5 * f;
	}

	private:
	
	Var<Real>* gc;

};

class BGCMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	SequenceAlignment* codondata;
	ContinuousData* contdata;
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
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;
	
	// autocorrelated process
	Dvar<PosReal>* PriorSigma;
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;

	Normal* mutgc;
	InstantStat* stationary;
	
	GTRRandomSubMatrixWithNormRates* nucmatrix;
	
	IIDNormal* MeanAALogFitness;
	int P; // number of degrees of freedom
	Dirichlet* weight;
	NormalIIDArray* iidarray;
	Normal* bgc;

	RandomMGMixBGCAAMutSelCodonSubMatrix* codonmatrix;

	// phylo process
	OneMatrixRASPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	BGCMutSelModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal)	{
		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();

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

		// a tree with all branch lengths iid from an exponential distribution of mean meanlength
		// meanlength is itself endowed with an exponential prior of mean 0.1
		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu); 
		mu->ClampAt(1);
		chronogram = new Chronogram(tree,mu);

		// a log normal process on that tree for the variations of the mutation rate
		PriorSigma = new Const<PosReal>(1);		
		sigma = new Gamma(One,PriorSigma);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma);
		
		// another log normal process for the variations of pop size

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		mutgc = new Normal(Zero,One);
		stationary = new InstantStat(mutgc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);

		MeanAALogFitness = new IIDNormal(Naa,Zero,One);
		MeanAALogFitness->SetAtZero();
		MeanAALogFitness->Clamp();


		cerr << "iidarray\n";
		P = inP;
		weight = new Dirichlet(P);
		iidarray = new NormalIIDArray(Naa,P,Zero,One);
		for (int p=0; p<P; p++)	{
			for (int k=0; k<Naa; k++)	{
				(*iidarray->GetVal(p))[k] = 0.1 * Random::Uniform();
			}
		}

		bgc = new Normal(Zero,One);

		cerr << "codon matrix\n";
		cerr << ((CodonStateSpace*) codondata->GetStateSpace())->GetNstate() << '\n';
		codonmatrix = new RandomMGMixBGCAAMutSelCodonSubMatrix((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,MeanAALogFitness,iidarray,weight,bgc);

		cerr << "process\n";
		// a phylogenetic process
		phyloprocess = new OneMatrixRASPhyloProcess(lognormaltree, 0, codonmatrix, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(weight);
		Register();

		MakeScheduler();

		if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
			cerr << '\n';
			GetTree()->Print(cerr);
			cerr << '\n';
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~BGCMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LengthTree* GetChronogram() {return chronogram;}
	
	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += mu->GetLogProb();
		total += chronogram->GetLogProb();

		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();

		total += relrate->GetLogProb();
		total += mutgc->GetLogProb();

		total += MeanAALogFitness->GetLogProb();
		total += iidarray->GetLogProb();
		total += weight->GetLogProb();
		total += bgc->GetLogProb();
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
		scheduler.Register(new SimpleMove(chronogram,1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),30,"chrono");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),30,"lognormal");

		scheduler.Register(new ProfileMove(relrate,1,1),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.3,1),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.1,2),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,3),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.03),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),30,"relrates");

		scheduler.Register(new SimpleMove(mutgc,0.3),30,"mutgc");
		scheduler.Register(new SimpleMove(mutgc,0.1),30,"mutgc");
		scheduler.Register(new SimpleMove(mutgc,0.01),30,"mutgc");

		scheduler.Register(new SimpleMove(bgc,0.3),30,"bgc");
		scheduler.Register(new SimpleMove(bgc,0.1),30,"bgc");
		scheduler.Register(new SimpleMove(bgc,0.01),30,"bgc");

		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.3),30,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.1),30,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.01),30,"mutgc/bgc");

		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.1),100,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.05),100,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.02),100,"mutgc/bgc");

		scheduler.Register(new SimpleMove(MeanAALogFitness,1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.01),10,"mean");

		scheduler.Register(new ProfileMove(weight,1,1),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.3,2),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.3,1),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.1,2),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.1,4),20,"weight");
		scheduler.Register(new SimpleMove(weight,0.1),20,"weight");
		scheduler.Register(new SimpleMove(weight,0.01),20,"weight");

		scheduler.Register(new SimpleMove(iidarray,1),4,"iidarray/sigma");
		scheduler.Register(new SimpleMove(iidarray,0.1),4,"iidarray/sigma");
		scheduler.Register(new SimpleMove(iidarray,0.01),4,"iidarray/sigma");
	
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
		mu->Sample();
		// chronogram->Sample();
		sigma->Sample();
		sigma->setval(10);
		lognormaltree->Sample();
		mutgc->Sample();
		mutgc->setval(0);
		relrate->Sample();
		weight->Sample();
		// MeanAALogFitness->Sample();
		// iidarray->Sample();
		bgc->Sample();
		bgc->setval(0);
		codonmatrix->specialUpdate();
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

	double GetBGC()	{
		return bgc->val();
	}

	double GetMutGC()	{
		return mutgc->val();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tbgc\tmutgc\tgrandmean\tgrandvar\tweightent\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << GetBGC();
		os << '\t' << GetMutGC();
		os << '\t' << iidarray->GetGrandMean();
		os << '\t' << iidarray->GetGrandVar();
		os << '\t' << weight->GetEntropy();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *relrate << '\n';
		os << *mutgc << '\n';
		os << *MeanAALogFitness << '\n';
		os << *weight << '\n';
		os << *iidarray << '\n';
		os << *bgc << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *mutgc;
		is >> *MeanAALogFitness;
		is >> *weight;
		is >> *iidarray;
		is >> *bgc;
	}

};

#endif

