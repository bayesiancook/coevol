
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
#include "ProfileMGAAMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"
#include "Normal.h"

class AAMutSelMixtureRandomMatrix : public MixtureRandomMatrix<Profile>	{

	public:

	AAMutSelMixtureRandomMatrix(Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* inmatrix)	{
		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = inmatrix;
		CreateRandomVariable();
		CreateRandomSubMatrix();
	}

	Rvar<Profile>* CreateRandomVariable()	{
		rvar = new Dirichlet(center,concentration);
		rvar->SetName("rvar / dirichlet");
		return rvar;
	}

	RandomSubMatrix* CreateRandomSubMatrix()	{
		matrix = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrix, rvar);
		matrix->SetName("aa mutsel matrix");
		return matrix;
	}
	
	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};


class AAMutSelMatrixFiniteMixture : public MatrixFiniteMixture<Profile>	{

	public:

	AAMutSelMatrixFiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix) :
			MatrixFiniteMixture<Profile>(insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	void setuniform()	{
		for (int k=0; k<GetComponentNumber(); k++)	{
			GetDirichlet(k)->setuniform();
		}
	}
	
	double* GetSiteStat(int i)	{
		return GetRandomVariable(i)->GetArray();
	}

	protected:

	Dirichlet* GetDirichlet(int k)	{
		Dirichlet* tmp = dynamic_cast<Dirichlet*>(GetComponent(k));
		if (! tmp)	{
			cerr << "error in AAMutSelMatrixFiniteMixture: component is not dirichlet\n";
			exit(1);
		}
		return tmp;
	}

	MixtureRandomMatrix<Profile>* CreateComponent(int k)	{
		return new AAMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};

class AAMutSelMatMixValMove : public MCUpdate	{

	public:

	AAMutSelMatMixValMove(AAMutSelMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	AAMutSelMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};


class AAMutSelModel : public ProbModel {

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

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ntaxa;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;
	Const<PosReal>* NAA;

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
	
	int P; // number of degrees of freedom

	Dirichlet* center;
	Exponential* concentration;

	public :

	AAMutSelMatrixFiniteMixture* aamutselmix;
	MatrixMixturePhyloProcess<Profile>* phyloprocess;

	// constructor
	// this is where the entire graph structure of the model is created

	AAMutSelModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal)	{
		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
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
		
		One = new Const<PosReal>(1.0);
		NAA = new Const<PosReal>(20.0);

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
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		concentration = new Exponential(NAA,Exponential::MEAN);
		center = new Dirichlet(Naa);

		aamutselmix = new AAMutSelMatrixFiniteMixture(Nsite,P,center,concentration,codonstatespace,nucmatrix);

		cerr << "create phylo process\n";
		phyloprocess = new MatrixMixturePhyloProcess<Profile>(lognormaltree,aamutselmix,codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(NAA);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(relrate);
		RootRegister(center);
		RootRegister(stationary);
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

	~AAMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return  Ntaxa;}
	int GetNstate() {return Nstate;}
	int GetNsite() {return Nsite;}
	CodonSequenceAlignment* GetCodonData() {return codondata;}
	double* GetSiteStat(int site) {return aamutselmix->GetSiteStat(site);}
	
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
		total += stationary->GetLogProb();

		total += center->GetLogProb();
		total += concentration->GetLogProb();

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

		scheduler.Register(new MatMixWeightAllocMove<Profile>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new MatMixValMove<Profile>(aamutselmix,0.3,10),1,"aamutsel mix simple val");
		scheduler.Register(new MatMixValMove<Profile>(aamutselmix,0.1,10),1,"aamutsel mix simple val");
		scheduler.Register(new MatMixValMove<Profile>(aamutselmix,0.01,10),1,"aamutsel mix simple val");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new AAMutSelMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new AAMutSelMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelMatMixValMove(aamutselmix,0.1,5,1),10,"aamutsel mix subset");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aamutselmix,1),1,"aamutsel mix weight alloc");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(concentration,1),10,"concentration");
		scheduler.Register(new SimpleMove(concentration,0.1),10,"concentration");

		scheduler.Register(new ProfileMove(center,0.1,1),10,"center");
		scheduler.Register(new ProfileMove(center,0.03,2),10,"center");
		scheduler.Register(new SimpleMove(center,0.01),10,"center");
		scheduler.Register(new SimpleMove(center,0.003),10,"center");

		scheduler.Register(new ProfileMove(relrate,1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.3,2),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.1,3),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.03),10,"relrates");
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
		// sigma->setval(10);
		lognormaltree->Sample();
		// lognormaltree->Reset();
		cerr << "stat\n";
		stationary->Sample();
		cerr << "relrate\n";
		relrate->Sample();

		center->Sample();
		concentration->Sample();

		aamutselmix->Sample();
		aamutselmix->setuniform();
	
		cerr << "iid\n";
		// iidarray->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
		cerr << "number of children for sigma: " << sigma->down.size() << '\n';
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

	double GetMeanEntropy()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*aamutselmix)[i]->GetEntropy();
		}
		mean /= Nsite;
		return mean;
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tmeanent\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << aamutselmix->GetEffSize();
		os << '\t' << GetMeanEntropy();
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
		os << *stationary << '\n';
		os << *center << '\n';
		os << *concentration << '\n';
		os << *aamutselmix << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *stationary;
		is >> *center;
		is >> *concentration;
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

