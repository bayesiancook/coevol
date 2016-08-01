
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
#include "Mixture.h"
#include "Normal.h"


class DirichletMixture : public FiniteMixture<Profile>	{

	public:

	DirichletMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration) : FiniteMixture<Profile>(insize, incomponentnumber)	{
		center = incenter;
		concentration = inconcentration;
		Create();
	}

	protected:

	Rvar<Profile>* CreateComponent(int k)	{
		return new Dirichlet(center,concentration);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;

};

class ProfileMixValMove : public MCUpdate	{

	public:

	ProfileMixValMove(DirichletMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	DirichletMixture* mix;
	double tuning;
	int N;
	int nrep;
};

class DirichletArray : public ValPtrArray< Rvar<Profile> >	{

	public: 

	DirichletArray(DirichletMixture* inmix, Var<PosReal>* inconcentration) : ValPtrArray<Rvar<Profile> >(inmix->GetSize()) {
		mix = inmix;
		concentration = inconcentration;
		Create();
	}

	Dirichlet* GetDirichlet(int site)	{
		Dirichlet* tmp = dynamic_cast<Dirichlet*>(GetVal(site));
		if (! tmp)	{
			cerr << "error in DirichletArray: invalid cast\n";
			exit(1);
		}
		return tmp;
	}

	double Move(double tuning)	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += this->GetVal(i)->Move(tuning);
		}
		return total / this->GetSize();
	}

	double Move(double tuning,n)	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += this->GetDirichlet(i)->Move(tuning,n);
		}
		return total / this->GetSize();
	}

	void drawSample()	{
		for (int i=0; i<this->GetSize(); i++)	{
			this->GetVal(i)->Sample();
		}
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<this->GetSize(); i++)	{
			total += this->GetVal(i)->GetLogProb();
		}
		return total ;
	}

	protected:

	Rvar<Profile>* CreateVal(int site)	{
		return new Dirichlet(mix->GetVal(site),concentration);
	}
};


class DirichletArrayMove : public MCUpdate	{

	public:

	DirichletArrayMove(DirichletArray* inarray, double intuning, int inN, int innrep) : array(inarray), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += array->Move(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	DirichletArray* array;
	double tuning;
	int N;
	int nrep;
};

class SiteAAMutSelModel : public ProbModel {

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

	public :

	Dirichlet* hypercenter;
	Exponential* hyperconcentration;
	Exponential* concentration;

	DirichletMixture* profilemix;
	DirichletArray* profilearray;
	RandomSubMatrix** matrixarray;
	SiteMatrixPhyloProcess* phyloprocess;

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
		
		One = new Const<PosReal>(1);
		NAA = new Const<PosReal>(20);

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

		hypercenter = new Dirichlet(Naa);
		hyperconcentration = new Exponential(NAA,Exponential::MEAN);
		concentration = new Exponential(NAA,Exponential::MEAN);
		profilemix = new DirichletMixture(Nsite,P,hypercenter,hyperconcentration);
		profilearray = new DirichletArray(profilemix,concentration);
		matrixarray = new RandomSubMatrix*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			matrixarray[i] = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrix, profilearray->GetDirichlet(i),0,0,0);
		}
		
		phyloprocess = new SiteMatrixPhyloProcess(lognormaltree,matrixarray,codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(NAA);
		RootRegister(hypercenter);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(profilemix->GetWeightVector());
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
	CodonSequenceAlignment* GetCodonData() {return codondata;}
	
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}

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

		total += profilemix->GetLogProb();
		total += profilearray->GetLogProb();
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

		scheduler.Register(new ProfileMove(hypercenter,0.01,2),10,"hypercenter");
		scheduler.Register(new ProfileMove(hypercenter,0.03,2),10,"hypercenter");
		scheduler.Register(new ProfileMove(hypercenter,0.01,5),10,"hypercenter");
		scheduler.Register(new SimpleMove(hypercenter,0.001),10,"hypercenter");

		scheduler.Register(new SimpleMove(hyperconcentration,1),100,"hyperconcentration");
		scheduler.Register(new SimpleMove(hyperconcentration,0.1),100,"hyperconcentration");

		scheduler.Register(new SimpleMove(concentration,1),100,"concentration");
		scheduler.Register(new SimpleMove(concentration,0.1),100,"concentration");

		scheduler.Register(new MixWeightAllocMove<Profile>(profilemix,1),1,"profile mix weight alloc");
		scheduler.Register(new MixValMove<Profile>(profilemix,0.3,10),1,"profile mix simple val");
		scheduler.Register(new MixValMove<Profile>(profilemix,0.1,10),1,"profile mix simple val");
		scheduler.Register(new ProfileMixValMove(profilemix,1,1,10),1,"profile mix subset");
		scheduler.Register(new ProfileMixValMove(profilemix,1,2,10),1,"profile mix subset");
		scheduler.Register(new ProfileMixValMove(profilemix,0.3,3,10),1,"profile mix subset");
		scheduler.Register(new ProfileMixValMove(profilemix,0.3,4,10),1,"profile mix subset");

		scheduler.Register(new DirichletArrayMove(profilearray,1,1,10),1,"profile array subset");
		scheduler.Register(new DirichletArrayMove(profilearray,1,2,10),1,"profile array subset");
		scheduler.Register(new DirichletArrayMove(profilearray,0.3,3,10),1,"profile array subset");
		scheduler.Register(new DirichletArrayMove(profilearray,0.3,4,10),1,"profile array subset");

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

		hypercenter->Sample();
		hyperconcentration->Sample();
		concentration->Sample();

		profilemix->Sample();
		profilearray->Sample();
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

	double GetMeanPriorEntropy()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*profilemix)[i]->GetEntropy();
		}
		mean /= Nsite;
		return mean;
	}

	double GetMeanEntropy()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += profilearray->GetVal(i)->GetEntropy();
		}
		mean /= Nsite;
		return mean;
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tmeanpriorent\tmeanent\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << profilemix->GetEffSize();
		os << '\t' << GetMeanPriorEntropy();
		os << '\t' << GetMeanEntropy();
		os << '\t' << stationary->GetMeanEntropy();
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
		os << *hypercenter << '\n';
		os << *hyperconcentration << '\n';
		os << *concentration << '\n';
		os << *profilemix << '\n';
		os << *profilearray << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *stationary;
		is >> *hypercenter;
		is >> *hyperconcentration;
		is >> *concentration;
		is >> *profilemix;
		is >> *profilearray;
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

