
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "DolloSubMatrix.h"
#include "PrecisionNormalTreeProcess.h"
#include "CalibratedChronogram.h"
#include "GeneralConjugatePath.h"
#include "Jeffreys.h"
#include "DolloSubMatrix.h"
#include "DolloPhyloProcess.h"

/*
A relatively simple model of gene loss:
- you first read a data matrix (presence/absence matrix) and a tree
- a chronogram (possibly, with fossil calibrations, but not mandatory)
- along this chronogram a log normal Brownian process describing the variation in the rate of loss over time (between lineages). These are the variations of the relative rate of loss.
- now, we also have an aboslute rate of loss (lossrate).
- and a DolloSubMatrix, which takes this rate of loss as its unique argument.
- finally, a "phyloprocess" is created: this is a simple i.i.d process: all genes (whose presence/absence profile is recorded in the data matrix) are lost at the same rates
*/

class DolloLogNormalModel : public ProbModel {

	// a fixed tree (read from file)
	Tree* tree;

	// a dataset (binary characters in that case)
	// state 1 : present
	// state 0 : absent
	SequenceAlignment* data;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (2 here)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	// constants
	Const<Real>* Zero;
	Const<PosReal>* One;

	// chronogram
	Chronogram* chronogram;

	// min and max for truncated jeffreys priors
	double minjeff;
	double maxjeff;

	// sigma is the rate of variation per unit of time of a lognormal Brownian process
	Jeffreys* sigma;
	// the log normal brownian process
	LogNormalTreeProcess* lognormaltree;

	// rate of loss
	Jeffreys* lossrate;
	// a rate matrix for a Markov process
	// modeling gene loss
	RandomDolloSubMatrix* matrix;

	// phylo process
	PhyloProcess* phyloprocess;

	bool iscalib;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	bool isCalibrated()	{
		return iscalib;
	}

	DolloLogNormalModel(string datafile, string treefile, string calibfile, double rootage, double rootstdev, bool sample)	{

		// fetch data from file
		data = new FileSequenceAlignment(datafile);
		cerr << "data ok\n";
		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();
		if (Nstate != 2)	{
			cerr << "error : Dollo works only on binary characters\n";
			exit(1);
		}

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

		if (calibfile != "None")	{
			iscalib = true;

			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
		}
		else	{
			iscalib = false;
			chronogram = new Chronogram(tree,One);
		}

		minjeff = 1e-3;
		maxjeff = 1e3;
		// a log normal process on that tree
		sigma = new Jeffreys(minjeff,maxjeff,One);
		// why?
		// sigma->ClampAt(1);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		lognormaltree->Reset();
		lognormaltree->Clamp();

		lossrate = new Jeffreys(minjeff,maxjeff,One);
		matrix = new RandomDolloSubMatrix(lossrate,false);

		// a phylogenetic process
		phyloprocess = new DolloOneMatrixPhyloProcess(lognormaltree,matrix,data);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		// register all the roots of the graphical model
		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(lognormaltree->GetRootRate());

		// recursive traversal of the model, starting from the roots
		// to register all the random variables
		Register();

		// the MCMC scheduler
		MakeScheduler();

		// update the model
		Update();
		TraceHeader(cerr);
		Trace(cerr);

		cerr << "model created\n";

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~DolloLogNormalModel() {}

	// various accessors
	LogNormalTreeProcess* GetLogNormalProcess()	{
		return lognormaltree;
	}

	Tree* GetTree()	{
		return tree;
	}

	Chronogram* GetChronogram()	{
		return chronogram;
	}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	Var<PosReal>* GetScale()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	// define here the series of MCMC updates, and their tuning parameters
	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");

		if (isCalibrated())	{
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
		}

		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.01),10,"sigma");

		scheduler.Register(new SimpleMove(lossrate,1),10,"sigma");
		scheduler.Register(new SimpleMove(lossrate,0.1),10,"sigma");
		scheduler.Register(new SimpleMove(lossrate,0.01),10,"sigma");

	}

	// Draw a sample from the prior
	void drawSample()	{
		cerr << "sample\n";
		chronogram->Sample();
		sigma->Sample();
		lognormaltree->Sample();
		lossrate->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;
		total += chronogram->GetLogProb();
		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();
		total += lossrate->GetLogProb();
		return total;
	}

	// in fact, returns the log prob path
	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// various summary statistics
	// used to check mcmc convergence
	double GetMeanRate()	{
		return lognormaltree->GetMeanRate();
	}

	double GetVarRate()	{
		return lognormaltree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlossrate\trootage\tsigma\tmeanrate\tvarrate";
		os << '\n';

	}
	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << lossrate->val();
		os << '\t' << sigma->val();
		os << '\t' << GetRootAge();
		os << '\t' << GetMeanRate() << '\t' << GetVarRate();
		os << '\n';
		os.flush();
	}

	// saving the parameter configuration
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *lossrate << '\n';
	}

	// getting the parameter configuration from file
	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		is >> *sigma;
		is >> *lognormaltree;
		is >> *lossrate;
	}

};
