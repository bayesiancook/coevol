
// a chronogram
//
// 2 independent log normal processes
// - log mutation rate (per unit of time)
// - log popsize (meant to be clamped at two leaf nodes)
// log(omega) = beta + alpha * log(popsize)

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
#include "GTRSubMatrix.h"
#include "MGAAMutSelCodonSubMatrix.h"

#include "BranchMatrixMixture.h"
#include "PrecisionNormalTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"

#include "ACGTProcess.h"
#include "Normal.h"

class AAMutSelMixtureRandomMatrixTree : public MixtureRandomMatrixTree<RealVector>	{

	public:

	AAMutSelMixtureRandomMatrixTree(Var<Real>* inmean, Var<PosReal>* invar, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* inmatrixtree, LengthTree* inpopsizetree)	{
		mean = inmean;
		var = invar;
		statespace = instatespace;
		nucmatrixtree = inmatrixtree;
		popsizetree = inpopsizetree;
		CreateRandomVariable();
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return nucmatrixtree->GetTree();}

	Rvar<RealVector>* CreateRandomVariable()	{
		rvar = new IIDNormal(Naa,mean,var);
		rvar->SetName("rvar / iid normal");
		return rvar;
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		RandomSubMatrix* matrix = 0;
		if (popsizetree)	{
			matrix = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrixtree->GetBranchVal(link->GetBranch()), rvar,0,0,popsizetree->GetBranchVal(link->GetBranch()));
		}
		else	{
			matrix = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrixtree->GetBranchVal(link->GetBranch()), rvar,0,0,0);
		}
		matrix->SetName("matrix");
		return matrix;
	}
	
	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	LengthTree* popsizetree;
};


class AAMutSelBranchMatrixFiniteMixture : public BranchMatrixFiniteMixture<RealVector>	{

	public:

	AAMutSelBranchMatrixFiniteMixture(int insize, int incomponentnumber, Var<Real>* inmean, Var<PosReal>* invar, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* innucmatrixtree, LengthTree* inpopsizetree) :
			BranchMatrixFiniteMixture<RealVector>(innucmatrixtree->GetTree(), insize, incomponentnumber)	{

		mean = inmean;
		var = invar;
		statespace = instatespace;
		nucmatrixtree = innucmatrixtree;
		popsizetree = inpopsizetree;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetIIDNormal(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	Tree* GetTree() {return nucmatrixtree->GetTree();}
	
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
		return new AAMutSelMixtureRandomMatrixTree(mean,var,statespace,nucmatrixtree,popsizetree);
	}

	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	LengthTree* popsizetree;
};


class AAMutSelBranchMatMixValMove : public MCUpdate	{

	public:

	AAMutSelBranchMatMixValMove(AAMutSelBranchMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	AAMutSelBranchMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};

class ACGTPopSizeAAMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	CodonSequenceAlignment* datacopy;
	ContinuousData* contdata;
	TaxonSet* taxonset;
	CodonStateSpace* codonstatespace;

	bool withgc;
	bool withpopsize;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ntaxa;

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

	// pop size
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* popsizetree;

	// parameters for the multivariate tree process describing variations of ATGC
	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* Sigma;
	
	// NodeMultiVariateTree
	ConjugateMultiVariateTreeProcess* process;

	Dirichlet* relrate;
	StatTree* stattree;
	BranchStatTree* branchstattree;
	NucMatrixTree* nucmatrixtree;
	
	int P; // number of degrees of freedom

	public :

	AAMutSelBranchMatrixFiniteMixture* aamutselmix;
	BranchMatrixMixturePhyloProcess<RealVector>* phyloprocess;

	// constructor
	// this is where the entire graph structure of the model is created

	ACGTPopSizeAAMutSelModel(string datafile, string treefile, int inP, bool inwithgc, bool inwithpopsize, bool sample=true, GeneticCodeType type=Universal)	{

		bool dc = false;

		withgc = inwithgc;
		withpopsize = inwithpopsize;

		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		if (dc)	{
			cerr << "before deleting constant sites : " << codondata->GetNsite() << '\n';
			codondata->DeleteAAConstantSites();
			ofstream os("woconst.ali"); 
			codondata->ToStream(os);
			os.close();
			cerr << "after deleting constant sites : " << codondata->GetNsite() << '\n';
			exit(1);
		}
		datacopy = 0;
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();
		Ntaxa = taxonset->GetNtaxa();

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
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		
		// another log normal process for the variations of pop size
		PriorTau = new Const<PosReal>(1);		
		tau = new Gamma(One,PriorTau);
		popsizetree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		popsizetree->GetRootRate()->ClampAt(0);
		
		GammaDiag = new Gamma*[Nnuc];
		diag = new Rvar<PosReal>*[Nnuc];
		for (int k=0; k<Nnuc; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(1);
		}
		priorOnSigmaZero = new RvarVec(diag,Nnuc);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		Sigma = new ConjugateInverseWishart(sigmaZero,Nnuc + 3);

		process = new ConjugateMultiVariateTreeProcess(Sigma,chronogram);
		
		stattree = new StatTree(process,Nnuc,0);
		branchstattree = new BranchStatTree(stattree);
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		nucmatrixtree = new NucMatrixTree(relrate,branchstattree);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		if (withpopsize)	{
			aamutselmix = new AAMutSelBranchMatrixFiniteMixture(Nsite,P,Zero,One,codonstatespace,nucmatrixtree,popsizetree);
		}
		else	{
			aamutselmix = new AAMutSelBranchMatrixFiniteMixture(Nsite,P,Zero,One,codonstatespace,nucmatrixtree,0);
		}

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixMixturePhyloProcess<RealVector>(lognormaltree,aamutselmix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorTau);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(aamutselmix->GetWeightVector());
		RootRegister(aamutselmix->GetRoot());
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

	~ACGTPopSizeAAMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetPopSizeTree() {return popsizetree;}
	BranchStatTree* GetBranchStatTree() {return branchstattree;}
	StatTree* GetStatTree() {return stattree;}
	LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return Ntaxa;}
	int GetNstate() {return Nstate;}
	CodonSequenceAlignment* GetCodonData() {return codondata;}
	
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

		for (int k=0; k<Nnuc; k++)	{
			total += GammaDiag[k]->GetLogProb();
		}
		total += Sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();

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
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(popsizetree,1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.01),10,"popsize");

		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,10,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.1,10),1,"conjugate sigma - process");

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
		tau->Sample();
		tau->setval(10);
		popsizetree->Sample();
		cerr << "relrate\n";
		relrate->Sample();
		cerr << "gc\n";
		cerr << "multivariate process\n";
		for (int k=0; k<Nnuc; k++)	{
			GammaDiag[k]->Sample();
			GammaDiag[k]->setval(1);
		}
		Sigma->Sample();
		Sigma->SetIdentity();
		process->Sample();
		process->Reset();

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

	double GetMeanLogPopSize()	{
		return popsizetree->GetMeanLogRate();
	}
	
	double GetVarLogPopSize()	{
		return popsizetree->GetVarLogRate();
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

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tacgt\tvar\tlogpopsize\tvar\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << aamutselmix->GetEffSize();
		os << '\t' << GetGrandMeanLogFitness();
		os << '\t' << GetGrandVarLogFitness();
		// os << '\t' << aamutselmix->GetWeightEntropy();
		os << '\t' << stattree->GetMeanEntropy();
		os << '\t' << stattree->GetVarEntropy();
		os << '\t' << GetMeanLogPopSize();
		os << '\t' << GetVarLogPopSize();
		os << '\t' << relrate->val().GetEntropy();
		for (int k=0; k<Nnuc; k++)	{
			for (int l=k+1; l<Nnuc; l++)	{
				os << '\t' << (*Sigma)[k][l];
			}
		}
		for (int k=0; k<Nnuc; k++)	{
			os << '\t' << (*Sigma)[k][k];
		}
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
		for (int k=0; k<Nnuc; k++)	{
			os << *GammaDiag[k] << '\n';
		}
		os << *Sigma << '\n';
		os << *process << '\n';
		os << *relrate << '\n';
		os << *aamutselmix << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *tau;
		is >> *popsizetree;
		for (int k=0; k<Nnuc; k++)	{
			is >> *GammaDiag[k];
		}
		is >> *Sigma;
		is >> *process;
		is >> *relrate;
		is >> *aamutselmix;
	}

	void MakeDataCopy()	{
		datacopy = new CodonSequenceAlignment(codondata);
	}

	void PostPredSample()	{
		phyloprocess->PostPredSample();
		if (! datacopy)	{
			MakeDataCopy();
		}
		phyloprocess->GetLeafData(datacopy);
	}

	CodonSequenceAlignment* GetPostPredCodonData() {return datacopy;}
};

