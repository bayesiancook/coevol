
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
#include "MeanExpTree.h"
#include "IncompleteGamma.h"

class LogMutSelOmega : public Dvar<Real>	{

	public:

	LogMutSelOmega(Var<PosReal>* inalpha, Var<PosReal>* inbeta, Var<Real>* inlogpopeff, int inN = 1000, double inmax = 100)	{
		alpha = inalpha;
		beta = inbeta;
		logpopeff = inlogpopeff;
		N = inN;
		max = inmax;
		CreateQuadrature();
		Register(alpha);
		Register(beta);
		Register(logpopeff);
	}

	protected:

	void specialUpdate()	{
		UpdateQuadrature();
		double p = AverageFixProbRatio(0);
		if ((! p>0) || (isinf(p)) || (isnan(p)))	{
			cerr << "error in MutSelOmega::specialUpdate : " << p << '\n';
		}
		setval(log(p));
	}

	private:

	void CreateQuadrature()	{

		x = new double[N+1];
		w = new double[N+1];
		y = new double[N];
		v = new double[N];
	}

	void UpdateQuadrature()	{

		double lnga1 = Random::logGamma(alpha->val());
		double Beta = beta->val() / exp(logpopeff->val());
		x[0] = 0;
		w[0] = 0;
		for (int i=1; i<=N; i++) {
			x[i] = max * exp(log(((double) i) / N) / alpha->val());
			w[i] = IncompleteGamma(x[i]*Beta,alpha->val(),lnga1);
		}
		for (int i=0; i<N; i++) {
			y[i] = (x[i+1] + x[i]) / 2;
			v[i] = w[i+1] - w[i];
		}
	}
	 
	double FixProbRatio(double x)	{

		if (fabs(x)<1e-10)	{
			return 1;
		}
		if (x > 100)	{
			return 0;
		}
		if (x < -100)	{
			return -x;
		}
		return x / (exp(x) -1);
	}


	double AverageFixProbRatio(double s0)	{

		double total = 0;
		for (int i=0; i<N; i++) {
			total += v[i] * FixProbRatio(y[i]+s0);
		}
		if (isinf(total))	{
			cerr << "inf prob in MGOmegaGCCodoSubMatrix\n";
			cerr << "rescaled bgc : " << s0 << '\n';
			cerr << '\n';
			for (int i=0; i<N; i++) {
				cerr << v[i] << '\t' << y[i] << '\t' << FixProbRatio(y[i]+s0) << '\n';
			}
			exit(1);
		}
		return total;
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	Var<Real>* logpopeff;

	int N;
	double max;
	double* x;
	double* w;
	double* y;
	double* v;

};

class LogOmegaTree : public NodeValPtrTree<Dvar<Real> >	{


	public:

	LogOmegaTree(NodeVarTree<Real>* inpopefftree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, int inN = 1000, double inmax = 100) {
		popefftree = inpopefftree;
		N = inN;
		max = inmax;
		alpha = inalpha;
		beta = inbeta;
		RecursiveCreate(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	protected:

	void specialUpdate(Link* from)	{
		GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
		
	Dvar<Real>* CreateNodeVal(const Link* link)	{
		return new LogMutSelOmega(alpha, beta, popefftree->GetNodeVal(link->GetNode()),N,max);
	}

	Tree* GetTree() {return popefftree->GetTree();}

	private:

	int N;
	double max;

	NodeVarTree<Real>* popefftree;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};

class MatrixTree : public BranchValPtrTree<RandomMGOmegaCodonSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree) {
		SetWithRoot(true);
		nucmatrix = innucmatrix;
		omegatree = inomegatree;
		statespace = instatespace;
		rootomega = new Const<PosReal>(1.0);
		RecursiveCreate(GetRoot());
	}

	Var<PosReal>* GetRootOmega() {return rootomega;}

	protected:

	RandomMGOmegaCodonSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
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

class PopSizeModel : public ProbModel {

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
	GTRRandomSubMatrixWithNormRates* nucmatrix;
	
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* popsizetree;

	Gamma* alpha;
	Gamma* beta;

	LogOmegaTree* logomegatree;
	MeanExpTree* omegatree;
	
	MatrixTree* matrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	PopSizeModel(string datafile, string treefile, string contdatafile, bool sample=true, GeneticCodeType type=Universal)	{
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
		popsizetree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		alpha = new Gamma(One,One);
		beta = new Gamma(One,One);

		if (contdata)	{
			popsizetree->ClampAt(contdata,0);
		}
		else	{
			alpha->ClampAt(0.18);
			beta->ClampAt(0.00012);
		}

		logomegatree = new LogOmegaTree(popsizetree,alpha,beta);
		omegatree = new MeanExpTree(logomegatree,chronogram,MEAN);

		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree);

		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(lognormaltree, matrixtree, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(PriorTau);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(popsizetree->GetRootRate());
		RootRegister(matrixtree->GetRootOmega());
		RootRegister(relrate);
		RootRegister(stationary);
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

	~PopSizeModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetPopSizeTree() {return popsizetree;}
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

		total += tau->GetLogProb();
		total += popsizetree->GetLogProb();

		total += alpha->GetLogProb();
		total += beta->GetLogProb();

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
		scheduler.Register(new SimpleMove(chronogram,1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),30,"chrono");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),30,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(popsizetree,10),100,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,1),100,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.1),100,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.01),100,"popsize");

		scheduler.Register(new SimpleMove(alpha,1),30,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),30,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.01),30,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.001),30,"alpha");

		scheduler.Register(new SimpleMove(beta,1),30,"beta");
		scheduler.Register(new SimpleMove(beta,0.1),30,"beta");
		scheduler.Register(new SimpleMove(beta,0.01),30,"beta");
		scheduler.Register(new SimpleMove(beta,0.001),30,"beta");

		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

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

		tau->Sample();
		tau->setval(10);
		popsizetree->Sample();

		alpha->Sample();
		beta->Sample();
		if (contdata)	{
			alpha->setval(0);
			beta->setval(0);
		}

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

	double GetMeanLogPopSize()	{
		return popsizetree->GetMeanLogRate();
	}
	
	double GetVarLogPopSize()	{
		return popsizetree->GetVarLogRate();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}
	
	/*
	double GetVarOmega()	{
		return omegatree->GetVarRate();
	}
	*/

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetAlpha()	{
		return alpha->val();
	}

	double GetBeta()	{
		return  beta->val();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tomega\tmeanlogpopsize\tvarlogpopsize\talpha\tbeta\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << GetMeanOmega();
		os << '\t' << GetMeanLogPopSize() << '\t' << GetVarLogPopSize();
		os << '\t' << *alpha << '\t' << *beta;
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
		os << *popsizetree << '\n';
		os << *alpha << '\n';
		os << *beta << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *tau;
		is >> *popsizetree;
		is >> *alpha;
		is >> *beta;
		is >> *relrate;
		is >> *stationary;
	}

};

#endif

