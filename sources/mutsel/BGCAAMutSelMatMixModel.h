
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
// #include "PhyloProcess.h"
#include "GTRSubMatrix.h"
#include "MGAAMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"

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

class BGCAAMutSelMixtureRandomMatrix : public MixtureRandomMatrix<RealVector>	{

	public:

	BGCAAMutSelMixtureRandomMatrix(Var<Real>* inmean, Var<PosReal>* invar, Var<Real>* inbgc, CodonStateSpace* instatespace, RandomSubMatrix* inmatrix)	{
		mean = inmean;
		var = invar;
		bgc = inbgc;
		statespace = instatespace;
		nucmatrix = inmatrix;
		CreateRandomVariable();
		CreateRandomSubMatrix();
	}

	Rvar<RealVector>* CreateRandomVariable()	{
		rvar = new IIDNormal(Naa,mean,var);
		rvar->SetName("rvar / iid normal");
		return rvar;
	}

	RandomSubMatrix* CreateRandomSubMatrix()	{
		matrix = new RandomMGBGCAAMutSelCodonSubMatrix(statespace, nucmatrix, rvar, bgc);
		matrix->SetName("aa mutsel matrix");
		return matrix;
	}
	
	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<Real>* bgc;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};


class BGCAAMutSelMatrixFiniteMixture : public MatrixFiniteMixture<RealVector>	{

	public:

	BGCAAMutSelMatrixFiniteMixture(int insize, int incomponentnumber, Var<Real>* inmean, Var<PosReal>* invar, Var<Real>* inbgc,  CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix) :
			MatrixFiniteMixture<RealVector>(insize, incomponentnumber)	{

		mean = inmean;
		var = invar;
		bgc = inbgc;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetIIDNormal(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	
	protected:

	IIDNormal* GetIIDNormal(int k)	{
		IIDNormal* tmp = dynamic_cast<IIDNormal*>(GetComponent(k));
		if (! tmp)	{
			cerr << "error in AAMutSelMatrixFiniteMixture: component is not iid normal\n";
			exit(1);
		}
		return tmp;
	}

	MixtureRandomMatrix<RealVector>* CreateComponent(int k)	{
		return new BGCAAMutSelMixtureRandomMatrix(mean,var,bgc,statespace,nucmatrix);
	}

	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<Real>* bgc;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};

class BGCAAMutSelMatMixValMove : public MCUpdate	{

	public:

	BGCAAMutSelMatMixValMove(BGCAAMutSelMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	BGCAAMutSelMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
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
	CodonStateSpace* codonstatespace;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;
	Const<Real>* Zero;

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

	Normal* bgc;
	
	int P; // number of degrees of freedom

	public :

	BGCAAMutSelMatrixFiniteMixture* aamutselmix;
	MatrixMixturePhyloProcess<RealVector>* phyloprocess;

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

		double* freq = new double[Nnuc];
		nucdata->GetEmpiricalFreq(freq);
		for (int k=0; k<Nnuc; k++)	{
			cerr << nucdata->GetStateSpace()->GetState(k) << '\t' << freq[k] << '\n';
		}
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
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		
		// another log normal process for the variations of pop size

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		mutgc = new Normal(Zero,One);
		stationary = new InstantStat(mutgc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		bgc = new Normal(Zero,One);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		aamutselmix = new BGCAAMutSelMatrixFiniteMixture(Nsite,P,Zero,One,bgc,codonstatespace,nucmatrix);

		cerr << "create phylo process\n";
		phyloprocess = new MatrixMixturePhyloProcess<RealVector>(lognormaltree,aamutselmix,codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(aamutselmix->GetWeightVector());
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


	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LengthTree* GetChronogram() {return chronogram;}
	
	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,false,true);
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

		total += bgc->GetLogProb();

		total += aamutselmix->GetLogProb();
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

		scheduler.Register(new MatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new MatMixValMove<RealVector>(aamutselmix,0.3,10),1,"aamutsel mix simple val");
		scheduler.Register(new MatMixValMove<RealVector>(aamutselmix,0.1,10),1,"aamutsel mix simple val");
		scheduler.Register(new BGCAAMutSelMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCAAMutSelMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCAAMutSelMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCAAMutSelMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");

		scheduler.Register(new SimpleMove(mutgc,0.3),10,"mutgc");
		scheduler.Register(new SimpleMove(mutgc,0.1),10,"mutgc");
		scheduler.Register(new SimpleMove(mutgc,0.01),10,"mutgc");

		scheduler.Register(new SimpleMove(bgc,0.3),10,"bgc");
		scheduler.Register(new SimpleMove(bgc,0.1),10,"bgc");
		scheduler.Register(new SimpleMove(bgc,0.01),10,"bgc");

		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.3),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(mutgc,bgc,0.01),10,"mutgc/bgc");

		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.05),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(mutgc,bgc,0.02),10,"mutgc/bgc");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.003),10,"relrates");

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
		cerr << "stat\n";
		mutgc->Sample();
		mutgc->setval(0);
		cerr << "relrate\n";
		relrate->Sample();
		bgc->Sample();
		bgc->setval(0);
		aamutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
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

	double GetGrandMeanLogFitness()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*aamutselmix)[i]->GetMean();
		}
		mean /= Nsite;
		return mean;
	}

	double GetGrandVarLogFitness()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*aamutselmix)[i]->GetVar();
		}
		mean /= Nsite;
		return mean;
	}

	double GetBGC()	{
		return bgc->val();
	}

	double GetMutGC()	{
		return mutgc->val();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tbgc\tmutgc\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << aamutselmix->GetEffSize();
		os << '\t' << GetGrandMeanLogFitness();
		os << '\t' << GetGrandVarLogFitness();
		os << '\t' << *mutgc;
		os << '\t' << *bgc;
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
		os << *bgc << '\n';
		os << *aamutselmix << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *mutgc;
		is >> *bgc;
		is >> *aamutselmix;
	}
};

#endif

