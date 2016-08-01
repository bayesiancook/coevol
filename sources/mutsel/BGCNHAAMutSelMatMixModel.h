
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
#include "BranchMatrixMixture.h"
#include "Normal.h"

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

class BGCAAMutSelMixtureRandomMatrixTree : public MixtureRandomMatrixTree<RealVector>	{

	public:

	BGCAAMutSelMixtureRandomMatrixTree(Tree* intree, Var<Real>* inmean, Var<PosReal>* invar, Var<Real>* inrootbgc, Var<Real>* intreebgc, CodonStateSpace* instatespace, RandomSubMatrix* inrootmatrix, RandomSubMatrix* intreematrix)	{
		mean = inmean;
		var = invar;
		rootbgc = inrootbgc;
		treebgc = intreebgc;
		statespace = instatespace;
		rootnucmatrix = inrootmatrix;
		treenucmatrix = intreematrix;
		tree = intree;
	//	SetWithRoot(true);
		CreateRandomVariable();
		rootmatrix = new RandomMGBGCAAMutSelCodonSubMatrix(statespace, rootnucmatrix, rvar, rootbgc);
		treematrix = new RandomMGBGCAAMutSelCodonSubMatrix(statespace, treenucmatrix, rvar, treebgc);
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return tree;}

	Rvar<RealVector>* CreateRandomVariable()	{
		rvar = new IIDNormal(Naa,mean,var);
		rvar->SetName("rvar / iid normal");
		return rvar;
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return rootmatrix;
		}
		return treematrix;
	}
	
	
	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<Real>* rootbgc;
	Var<Real>* treebgc;
	Tree* tree;
	CodonStateSpace* statespace;
	RandomSubMatrix* rootnucmatrix;
	RandomSubMatrix* treenucmatrix;
	RandomMGBGCAAMutSelCodonSubMatrix* rootmatrix;
	RandomMGBGCAAMutSelCodonSubMatrix* treematrix;
};


class BGCNHAAMutSelMatrixFiniteMixture : public BranchMatrixFiniteMixture<RealVector>	{

	public:

	BGCNHAAMutSelMatrixFiniteMixture(Tree* intree, int insize, int incomponentnumber, Var<Real>* inmean, Var<PosReal>* invar, Var<Real>* inrootbgc, Var<Real>* intreebgc, CodonStateSpace* instatespace, RandomSubMatrix* inrootnucmatrix, RandomSubMatrix* intreenucmatrix) :
			BranchMatrixFiniteMixture<RealVector>(intree, insize, incomponentnumber)	{

		mean = inmean;
		var = invar;
		rootbgc = inrootbgc;
		treebgc = intreebgc;
		statespace = instatespace;
		rootnucmatrix = inrootnucmatrix;
		treenucmatrix = intreenucmatrix;
		tree = intree;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetIIDNormal(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	Tree* GetTree() {return tree;}
	
	protected:

	IIDNormal* GetIIDNormal(int k)	{
		IIDNormal* tmp = dynamic_cast<IIDNormal*>(GetComponent(k));
		if (! tmp)	{
			cerr << "error in AAMutSelMatrixFiniteMixture: component is not iid normal\n";
			exit(1);
		}
		return tmp;
	}

	MixtureRandomMatrixTree<RealVector>* CreateComponent(int k)	{
		return new BGCAAMutSelMixtureRandomMatrixTree(tree,mean,var,rootbgc,treebgc,statespace,rootnucmatrix,treenucmatrix);
	}

	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<Real>* rootbgc;
	Var<Real>* treebgc;
	CodonStateSpace* statespace;
	RandomSubMatrix* rootnucmatrix;
	RandomSubMatrix* treenucmatrix;
	Tree* tree;
};

class BGCNHAAMutSelMatMixValMove : public MCUpdate	{

	public:

	BGCNHAAMutSelMatMixValMove(BGCNHAAMutSelMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	BGCNHAAMutSelMatrixFiniteMixture* mix;
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
	Normal* rootmutgc;
	Normal* treemutgc;
	InstantStat* rootstationary;
	InstantStat* treestationary;
	GTRRandomSubMatrixWithNormRates* rootnucmatrix;
	GTRRandomSubMatrixWithNormRates* treenucmatrix;

	Normal* rootbgc;
	Normal* treebgc;
	
	int P; // number of degrees of freedom

	public :

	BGCNHAAMutSelMatrixFiniteMixture* aamutselmix;
	BranchMatrixMixturePhyloProcess<RealVector>* phyloprocess;

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
		rootmutgc = new Normal(Zero,One);
		treemutgc = new Normal(Zero,One);
		rootstationary = new InstantStat(rootmutgc);
		treestationary = new InstantStat(treemutgc);
		rootnucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,rootstationary);
		treenucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,treestationary);

		rootbgc = new Normal(Zero,One);
		treebgc = new Normal(Zero,One);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		aamutselmix = new BGCNHAAMutSelMatrixFiniteMixture(tree,Nsite,P,Zero,One,rootbgc,treebgc,codonstatespace,rootnucmatrix,treenucmatrix);

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixMixturePhyloProcess<RealVector>(lognormaltree,aamutselmix,codondata);
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
		total += rootmutgc->GetLogProb();
		total += treemutgc->GetLogProb();

		total += rootbgc->GetLogProb();
		total += treebgc->GetLogProb();

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

		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.3,10),1,"aamutsel mix simple val");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.1,10),1,"aamutsel mix simple val");
		scheduler.Register(new BGCNHAAMutSelMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCNHAAMutSelMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCNHAAMutSelMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new BGCNHAAMutSelMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");

		scheduler.Register(new SimpleMove(rootmutgc,0.3),10,"mutgc");
		scheduler.Register(new SimpleMove(rootmutgc,0.1),10,"mutgc");
		scheduler.Register(new SimpleMove(rootmutgc,0.01),10,"mutgc");

		scheduler.Register(new SimpleMove(rootbgc,0.3),10,"bgc");
		scheduler.Register(new SimpleMove(rootbgc,0.1),10,"bgc");
		scheduler.Register(new SimpleMove(rootbgc,0.01),10,"bgc");

		scheduler.Register(new JointSimpleMove(rootmutgc,rootbgc,0.3),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(rootmutgc,rootbgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(rootmutgc,rootbgc,0.01),10,"mutgc/bgc");

		scheduler.Register(new AdditiveCompensatoryMove(rootmutgc,rootbgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(rootmutgc,rootbgc,0.05),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(rootmutgc,rootbgc,0.02),10,"mutgc/bgc");

		scheduler.Register(new SimpleMove(treemutgc,0.3),10,"mutgc");
		scheduler.Register(new SimpleMove(treemutgc,0.1),10,"mutgc");
		scheduler.Register(new SimpleMove(treemutgc,0.01),10,"mutgc");

		scheduler.Register(new SimpleMove(treebgc,0.3),10,"bgc");
		scheduler.Register(new SimpleMove(treebgc,0.1),10,"bgc");
		scheduler.Register(new SimpleMove(treebgc,0.01),10,"bgc");

		scheduler.Register(new JointSimpleMove(treemutgc,treebgc,0.3),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(treemutgc,treebgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new JointSimpleMove(treemutgc,treebgc,0.01),10,"mutgc/bgc");

		scheduler.Register(new AdditiveCompensatoryMove(treemutgc,treebgc,0.1),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(treemutgc,treebgc,0.05),10,"mutgc/bgc");
		scheduler.Register(new AdditiveCompensatoryMove(treemutgc,treebgc,0.02),10,"mutgc/bgc");

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
		cerr << "root\n";
		rootmutgc->Sample();
		rootmutgc->setval(0);
		cerr << "tree\n";
		treemutgc->Sample();
		treemutgc->setval(0);
		cerr << "relrate\n";
		relrate->Sample();
		rootbgc->Sample();
		rootbgc->setval(0);
		treebgc->Sample();
		treebgc->setval(0);
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

	double GetRootBGC()	{
		return rootbgc->val();
	}

	double GetRootMutGC()	{
		return rootmutgc->val();
	}

	double GetTreeBGC()	{
		return treebgc->val();
	}

	double GetTreeMutGC()	{
		return treemutgc->val();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tbgc\tmutgc\trootbgc\trootmutgc\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << aamutselmix->GetEffSize();
		os << '\t' << GetGrandMeanLogFitness();
		os << '\t' << GetGrandVarLogFitness();
		os << '\t' << *treemutgc;
		os << '\t' << *treebgc;
		os << '\t' << *rootmutgc;
		os << '\t' << *rootbgc;
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
		os << *rootmutgc << '\n';
		os << *rootbgc << '\n';
		os << *treemutgc << '\n';
		os << *treebgc << '\n';
		os << *aamutselmix << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *rootmutgc;
		is >> *rootbgc;
		is >> *treemutgc;
		is >> *treebgc;
		is >> *aamutselmix;
	}
};

#endif

