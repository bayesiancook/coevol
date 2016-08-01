
// a chronogram
// 2 independent log normal processes
// - log Ks
// - log omega = log Kn/Ks

#ifndef BRANCHOMEGA_H
#define BRANCHOMEGA_H


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
#include "ConjugateMultiVariateTreeProcess.h"


class InstantStat : public Dvar<Profile>	{

	public:

	InstantStat(Var<RealVector>* inup, Var<RealVector>* indown, int inindex = 0)	{
		setval(Profile(Nnuc));
		bkvalue = Profile(Nnuc);
		up = inup;
		down = indown;
		Register(up);
		Register(down);
		index = inindex;
	}

	void specialUpdate()	{
		double total = 0;
		double gcup = exp((*up)[i+index]) / (1 + exp((*up)[i+index]));
		double gcdown = exp((*down)[i+index]) / (1 + exp((*down)[i+index]));
		gc = 0.5 * (gcup + gcdown);
		(*this)[0] = (*this)[3] = 0.5 * (1 - gc);
		(*this)[1] = (*this)[2] = 0.5 * gc;
	}

	private:
	
	int index;
	Var<RealVector>* up;
	Var<RealVector>* down;
	double gc;

};

class StatTree : public BranchValPtrTree<InstantStat>	{

	public:

	StatTree(Var<Profile>* instat0, MultiVariateTreeProcess* inprocess, int inindex)	{
		SetWithRoot(true);
		index = inindex;
		stat0 = instat0;
		process = inprocess;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	protected:

	InstantStat* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new InstantStat(stat0,process->GetNodeVal(link->GetNode()),process->GetNodeVal(link->GetNode()));
		}
		return new InstantStat(stat0,process->GetNodeVal(link->GetNode()),process->GetNodeVal(link->Out()->GetNode()));
	}

	private:
	MultiVariateTreeProcess* process;
	int index;
	Var<Profile>* stat0;

};


class NucMatrixTree : public BranchValPtrTree<GTRRandomSubMatrixWithNormRates>	{


	public:

	NucMatrixTree(Var<Profile>* inrelrate, StatTree* instattree) {
		SetWithRoot(true);
		stattree = instattree;
		relrate = inrelrate;
		RecursiveCreate(GetRoot());
	}

	protected:

	GTRRandomSubMatrixWithNormRates* CreateBranchVal(const Link* link)	{
		return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	StatTree* stattree;
	Var<Profile>* relrate;

};

class MatrixTree : public BranchValPtrTree<RandomMGOmegaCodonSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* inomegatree) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = new Const<PosReal>(1.0);
		RecursiveCreate(GetRoot());
	}

	Var<PosReal>* GetRootOmega() {return rootomega;}

	protected:

	RandomMGOmegaCodonSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	NucMatrixTree* nucmatrixtree;
	Const<PosReal>* rootomega;

};

class BranchMatrixPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchMatrixPhyloProcess(LengthTree* intree, BranchValPtrTree<RandomMGOmegaCodonSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:
	BranchValPtrTree<RandomMGOmegaCodonSubMatrix>* matrixtree;
};

class BranchOmegaGCModel : public ProbModel {

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
	int Ncont;

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
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;

	StatTree* stattree;
	NucMatrixTree* nucmatrixtree;

	Gamma* one;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* Sigma;
	
	// NodeMultiVariateTree
	ConjugateMultiVariateTreeProcess* process;

	
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* omegatree;

	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* gentimetree;

	MatrixTree* matrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	BranchOmegaGCModel(string datafile, string treefile, string contdatafile, bool sample=true, GeneticCodeType type=Universal)	{
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

		cerr << contdatafile << '\n';
		if (contdatafile != "None")	{
			cerr << "read cont data \n";
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

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

		// a log normal process on that tree for the variations of the mutation rate
		PriorSigma = new Const<PosReal>(1);		
		sigma = new Gamma(One,PriorSigma);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma);
		
		// another log normal process for the variations of pop size
		PriorTau = new Const<PosReal>(1);		
		tau = new Gamma(One,PriorTau);
		omegatree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		
		PriorTheta = new Const<PosReal>(1);		
		theta = new Gamma(One,PriorTheta);
		gentimetree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		if (contdata)	{
			gentimetree->ClampAt(contdata,0);
		}
		
		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1

		diag = new Rvar<PosReal>*[Nnuc];
		one = new Gamma(One,One);
		one->ClampAt(1);
		for (int k=0; k<Nnuc; k++)	{
			diag[k] = one;
		}
		priorOnSigmaZero = new RvarVec(diag,Nnuc);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		Sigma = new ConjugateInverseWishart(sigmaZero,Nnuc+4);
		process = new ConjugateMultiVariateTreeProcess(Sigma,chronogram);
		process->ClampRoot();
	
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		stattree = new StatTree(stationary,process,0);
		nucmatrixtree = new NucMatrixTree(relrate,stattree);
		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree);

		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(lognormaltree, matrixtree, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(PriorTau);
		RootRegister(PriorTheta);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(omegatree->GetRootRate());
		RootRegister(gentimetree->GetRootRate());
		RootRegister(matrixtree->GetRootOmega());
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(one);
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

	~BranchOmegaGCModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetOmegaTree() {return omegatree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetGenTimeTree() {return gentimetree;}
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

		total += tau->GetLogProb();
		total += omegatree->GetLogProb();

		total += theta->GetLogProb();
		total += gentimetree->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		total += Sigma->GetLogProb();
		total += process->GetLogProb();

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

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(omegatree,1),10,"omega");
		scheduler.Register(new SimpleMove(omegatree,0.1),10,"omega");
		scheduler.Register(new SimpleMove(omegatree,0.01),10,"omega");

		scheduler.Register(new SimpleMove(theta,10),100,"theta");
		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),10,"theta");
		scheduler.Register(new SimpleMove(gentimetree,10),10,"gentime");
		scheduler.Register(new SimpleMove(gentimetree,1),10,"gentime");
		scheduler.Register(new SimpleMove(gentimetree,0.1),10,"gentime");
		scheduler.Register(new SimpleMove(gentimetree,0.01),10,"gentime");

		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,10,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.01,10),1,"conjugate sigma - process");

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
		sigma->setval(10);
		lognormaltree->Sample();

		tau->Sample();
		tau->setval(10);
		omegatree->Sample();

		theta->Sample();
		theta->setval(10);
		gentimetree->Sample();

		stationary->Sample();
		relrate->Sample();

		Sigma->Sample();
		Sigma->SetIdentity();
		process->Sample();
		process->Reset();

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

	double GetMeanLogGenTime()	{
		return gentimetree->GetMeanLogRate();
	}
	
	double GetVarLogGenTime()	{
		return gentimetree->GetVarLogRate();
	}

	double GetMeanOmega()	{
		return omegatree->GetMeanRate();
	}
	
	double GetVarOmega()	{
		return omegatree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tomega\tstderr\tsigma\tmeanrho\tvarrho\tmeangentime\tvargentime\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << GetMeanOmega();
		os << '\t' << sqrt(GetVarOmega());
		os << '\t' << sigma->val();
		os << '\t' << GetMeanRho() << '\t' << GetVarRho();
		os << '\t' << GetMeanLogGenTime() << '\t' << GetVarLogGenTime();
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
		os << *tau << '\n';
		os << *omegatree << '\n';
		os << *theta << '\n';
		os << *gentimetree << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
		os << *Sigma << '\n';
		os << *process << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *tau;
		is >> *omegatree;
		is >> *theta;
		is >> *gentimetree;
		is >> *relrate;
		is >> *stationary;
		is >> *Sigma;
		is >> *process;
	}

};

#endif

