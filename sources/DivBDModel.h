
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "Jeffreys.h"
#include "DivBDCalibratedChronogram.h"
#include "GTRSubMatrix.h"
#include "OneMatrixPhyloProcess.h"
#include "GeneralConjugatePath.h"

#include "AIS.h"


class PosSum : public Dvar<PosReal>	{

	public:

	PosSum(Var<PosReal>* ina, Var<PosReal>* inb)	{
		a = ina;
		b = inb;
		Register(a);
		Register(b);
		specialUpdate();
	}

	protected:

	void specialUpdate()	{
		setval(a->val() + b->val());
	}

	private:

	Var<PosReal>* a;
	Var<PosReal>* b;
};

class DivBDModel: public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
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

	Jeffreys* lambda;
	// Jeffreys* lambdamu;
	// Gamma* lambdamu;
	Jeffreys* mu;
	// PosSum* lambda;

	int divmodel;
	// 0 : mu = 0;
	// 1 : lambda mu free

	// chronogram
	DivBDCalibratedChronogram* chronogram;

	// autocorrelated process
	Jeffreys* sigma;
	LogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* matrix;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;
	// OneMatrixPhyloProcess* phyloprocess;
	bool pathconjugate;

	AISController* ais;
	int burnin;
	int every;
	double initlambda;
	double finallambda;
	double initmu;
	double finalmu;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	DivBDModel(string datafile, string treefile, string calibfile, int indivmodel, double ininitlambda, double infinallambda, double ininitmu, double infinalmu, bool inpathconjugate, bool sample)	{

		double minjeff = 1e-9;
		double maxjeff = 1e9;

		pathconjugate = inpathconjugate;
		divmodel = indivmodel;
		initlambda = ininitlambda;
		finallambda = infinallambda;
		initmu = ininitmu;
		finalmu = infinalmu;

		// fetch data from file
		data = new FileSequenceAlignment(datafile);
		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)

		taxonset = data->GetTaxonSet();

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

		lambda = new Jeffreys(minjeff,maxjeff,Zero);
		// lambdamu = new Jeffreys(minjeff,maxjeff,Zero);
		// lambdamu = new Gamma(One,One);
		mu = new Jeffreys(minjeff,maxjeff,Zero);
		lambda->setval(0.02);
		// lambdamu->setval(0.001);
		mu->setval(0.001);
		if (! divmodel)	{
			mu->ClampAt(0);
		}
		// lambda = new PosSum(mu,lambdamu);

		CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);
		chronogram = new DivBDCalibratedChronogram(tree,One,lambda,mu,calibset);

		// a log normal process on that tree
		sigma = new Jeffreys(minjeff,maxjeff,Zero);
		sigma->setval(1.0);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		// lognormaltree->Reset();

		// substitution matrix
		relrate = new Dirichlet(Nstate*(Nstate-1)/2);
		stationary = new Dirichlet(Nstate);
		relrate->setuniform();
		stationary->setuniform();
		// relrate->Clamp();
		// stationary->Clamp();
		matrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		// a phylogenetic process
		if (pathconjugate)	{
			pathconjtree = new OneMatrixPathConjugateTree(lognormaltree,matrix,data);
			phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
		}
		else	{
			pathconjtree = 0;
			phyloprocess = new OneMatrixPhyloProcess(lognormaltree,matrix,data);
		}
		cerr << "unfold\n";
		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(stationary);

		Register();

		MakeScheduler();

		Update();
		TraceHeader(cerr);
		Trace(cerr);

		if (IsAIS())	{
			MakeAISController();
		}

		cerr << "model created\n";

	}

	bool IsAIS()	{
		return ((initlambda != -1) && (finallambda != -1));
	}

	bool FixedLambda()	{
		return (initlambda != -1);
	}

	bool FixedMu()	{
		return (initmu != -1);
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~DivBDModel() {}


	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,false,false);
		return 1;
	}

	int GetNtaxa()	{
		return taxonset->GetNtaxa();
	}

	LogNormalTreeProcess* GetLogNormalProcess()	{
		return lognormaltree;
	}

	Tree* GetTree()	{
		return tree;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	Chronogram* GetChronogram()	{
		return chronogram;
	}

	DivBDCalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<DivBDCalibratedChronogram*>(chronogram);
	}

	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		// total += lambdamu->GetLogProb();
		if (divmodel)	{
			total += mu->GetLogProb();
		}
		total += chronogram->GetLogProb();
		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();
		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	void MakeAISController()	{
		ais = new AISController();
		ais->AISRegister(lambda,initlambda,finallambda);
		ais->AISRegister(mu,initmu,finalmu);
	}

	double AISSet(double p)	{
		double before = lambda->GetLogProb() + mu->GetLogProb();
		double ret = ais->AISSet(p);
		double after = lambda->GetLogProb() + mu->GetLogProb();
		return ret - (after - before);
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		if (pathconjugate)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		/*
		scheduler.Register(new SimpleMove(lambdamu,1),10,"lambda-mu");
		scheduler.Register(new SimpleMove(lambdamu,0.1),10,"lambda-mu");
		*/

		if (! FixedLambda())	{
			scheduler.Register(new SimpleMove(lambda,1),10,"lambda-mu");
			scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda-mu");
		}

		if (! FixedMu())	{
			scheduler.Register(new SimpleMove(mu,1),10,"mu");
			scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		}

		scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
		scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
		scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");

		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.01),10,"sigma");

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

	}

	// Metropolis Hastings, successively called on each component of the model
	// old fashioned move
	//

	// Draw a sample from the prior

	void drawSample()	{
		lambda->Sample();
		// lambdamu->Sample();
		if (divmodel)	{
			mu->Sample();
		}
		chronogram->Sample();
		sigma->Sample();
		lognormaltree->Sample();
		stationary->Sample();
		relrate->Sample();
		phyloprocess->Sample();
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetLambda()	{
		return lambda->val();
	}

	double GetMu()	{
		return mu->val();
	}

	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}

	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetRootAge()	{
		return GetCalibratedChronogram()->GetRootAge();
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlambda\tmu\t(lambda-mu)T\trootage\tsigma\tmeanrho\tvarrho\tcheckbounds";
		os << '\n';

	}
	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		// os << '\t' << GetCalibratedChronogram()->GetLogProb() << '\t' << GetCalibratedChronogram()->GetYuleLogProb();
		os << '\t' << lambda->val();
		os << '\t' << mu->val();
		os << '\t' << (lambda->val() - mu->val()) * GetRootAge();
		os << '\t' << GetRootAge();
		os << '\t' << sigma->val();
		os << '\t' << GetMeanRho() << '\t' << GetVarRho();
		os << '\t' << GetCalibratedChronogram()->CheckBounds();
		os << '\n';
		os.flush();
	}

	void WriteSuffStat(ostream& os)	{

		os << lambda->val() << '\t' << mu->val();
		list<double> agelist;
		GetCalibratedChronogram()->GetAgeList(agelist);
		for (list<double>::iterator i=agelist.begin(); i!=agelist.end(); i++)	{
			os << '\t' << (*i) * GetRootAge();
		}
		os << '\n';
	}

	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		// os << *lambdamu << '\n';
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *chronogram->GetScale() << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		// is >> *lambdamu;
		is >> *mu;
		is >> *chronogram;
		is >> *chronogram->GetScale();
		is >> *sigma;
		is >> *lognormaltree;
		is >> *relrate;
		is >> *stationary;
	}

};
