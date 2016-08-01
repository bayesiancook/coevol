
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "IID.h"
#include "ScalableNormalTreeProcess.h"
#include "BDCalibratedChronogram.h"

#include "GeneralConjugatePath.h"

template <class R, class D> class DSemiConjugateMove : public MCUpdate	{

	public:

	DSemiConjugateMove(R* inrandom, D* indsemi, double intuning, int inn) : random(inrandom), dsemi(indsemi), tuning(intuning), n(inn) {}

	double Move(double tuning_modulator = 1)	{
		dsemi->ActivateSufficientStatistic();
		double total = 0;
		for (int i=0; i<n; i++)	{
			total += random->Move(tuning* tuning_modulator);
		}
		total /= n;
		dsemi->InactivateSufficientStatistic();
		return total;
	}

	protected:

	R* random;
	D* dsemi;
	double tuning;
	int n;
};

class DSemiConjugateMappingMove : public MCUpdate	{

	public:

	DSemiConjugateMappingMove(PhyloProcess* inprocess, PathConjugateTree* inpathconjtree) : process(inprocess), pathconjtree(inpathconjtree) {}

	double Move(double tuning_modulator=1)	{
		pathconjtree->InactivateSufficientStatistic();
		process->Move(1);
		pathconjtree->ActivateSufficientStatistic();
		return 1;
	}

	protected:

	PhyloProcess* process;
	PathConjugateTree* pathconjtree;
};

class GTRLogNormalModel : public ProbModel {

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

	Const<PosReal>* One;

	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	double meanchi;
	double meanchi2;

	// chronogram
	Chronogram* chronogram;

	// autocorrelated process
	Gamma* sigma;
	Gamma* nu;
	LogNormalTreeProcess* lognormaltree;

	// substitution matrix is relrate * stationary
//	Dirichlet* relrate;
	IIDExp* relrate;
	Dirichlet* stationary;
	// GTRRandomSubMatrixWithNormRates* matrix;
	GTRRandomSubMatrix* matrix;

	// phylo process
	// OneMatrixPhyloProcess* phyloprocess;
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	bool conjpath;
	int chronoprior;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRLogNormalModel(string datafile, string treefile, string calibfile, double inrootage, double inrootstdev, int inchronoprior)	{

		conjpath = true;
		chronoprior = inchronoprior;
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

		One = new Const<PosReal>(1);


		double rootage = inrootage;
		double rootstdev = inrootstdev;

		double a = rootage * rootage / rootstdev / rootstdev;
		double b = rootage / rootstdev / rootstdev;
		RootAlpha = new Const<PosReal>(a);
		RootBeta = new Const<PosReal>(b);
		CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

		if (chronoprior)	{
			cerr << "BD\n";
			meanchi = 1e-3;
			meanchi2 = 1e-3;
			MeanChi = new Const<PosReal>(meanchi);
			MeanChi2 = new Const<PosReal>(meanchi2);
			Chi = new Exponential(MeanChi,Exponential::MEAN);
			Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
			Chi->ClampAt(1.5);
			Chi2->ClampAt(1.5);
			chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,RootAlpha,RootBeta,calibset);
		}
		else	{
			chronogram = new CalibratedChronogram(tree,One,RootAlpha,RootBeta,calibset);
		}

		// a log normal process on that tree
		sigma = new Gamma(One,One);
		nu = new Gamma(One,One);
		cerr << "log normal tree\n";
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,nu);
		cerr << "ok\n";

		// substitution matrix
		// relrate = new Dirichlet(Nstate*(Nstate-1)/2);
		relrate = new IIDExp(Nstate*(Nstate-1)/2);
		for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
			(*relrate)[i] = 1;
		}
		relrate->Clamp();
		stationary = new Dirichlet(Nstate);
		// matrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);
		matrix = new GTRRandomSubMatrix(relrate,stationary,false);

		// a phylogenetic process
		cerr << "phyloprocess\n";
		if (conjpath)	{
			pathconjtree = new OneMatrixPathConjugateTree(lognormaltree,matrix,data);
			phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
		}
		else	{
			pathconjtree = 0;
			phyloprocess = new OneMatrixPhyloProcess(lognormaltree,matrix,data);
		}
		cerr << "unfold\n";
		phyloprocess->Unfold();
		phyloprocess->Sample();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(RootAlpha);
		RootRegister(RootBeta);
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		MakeScheduler();

		Update();
		TraceHeader(cerr);
		Trace(cerr);

		cerr << "model created\n";

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GTRLogNormalModel() {}

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
		if (chronoprior)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();
		total += nu->GetLogProb();
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

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	BDCalibratedChronogram* GetBDCalibratedChronogram()	{
		return dynamic_cast<BDCalibratedChronogram*>(chronogram);
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{
		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		int nrep = conjpath ? 30 : 1;
		for (int i=0; i<nrep; i++)	{
			scheduler.Register(new SimpleMove(nu,1),10,"nu");
			scheduler.Register(new SimpleMove(nu,0.1),10,"nu");
			scheduler.Register(new SimpleMove(nu,0.01),10,"nu");
			scheduler.Register(new SimpleMove(sigma,10),10,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");

			if (chronoprior)	{
				scheduler.Register(new SimpleMove(Chi,10),10,"p1");
				scheduler.Register(new SimpleMove(Chi,1),10,"p1");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"p1");
				scheduler.Register(new SimpleMove(Chi2,10),10,"p2");
				scheduler.Register(new SimpleMove(Chi2,1),10,"p2");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"p2");
			}

			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");

			scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
			scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
			scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

			/*
			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
			*/
			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}
	}

	// Metropolis Hastings, successively called on each component of the model
	// old fashioned move
	//

	// Draw a sample from the prior

	void drawSample()	{
		if (chronoprior)	{
			Chi->Sample();
			Chi2->Sample();
		}
		chronogram->Sample();
		sigma->Sample();
		nu->Sample();
		lognormaltree->Sample();
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

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetRootAge()	{
		return GetCalibratedChronogram()->GetScale()->val();
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tscale";
		if (chronoprior)	{
			os << "\tp1\tp2";
		}
		os << "\tsigma\tmeanrho\tvarrho\n";
	}

	double GetLogBDPrior()	{
		return GetBDCalibratedChronogram()->GetLogProb() - GetBDCalibratedChronogram()->logJacobian();
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetRootAge();
		if (chronoprior)	{
			os << '\t' << *Chi << '\t' << *Chi2;
			os << '\t' << GetBDCalibratedChronogram()->GetTotalLogGG() + GetBDCalibratedChronogram()->logAbsBDUncalibratedPrior();
		}
		os << '\t' << sigma->val();
		os << '\t' << nu->val();
		os << '\t' << GetMeanRho() << '\t' << GetVarRho();
		os << '\t' << *lognormaltree->GetRootRate();
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
	}

	void FromStream(istream& is)	{
	}

};


int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string calibfile = argv[3];
	double rootage = atof(argv[4]);
	double rootstdev = atof(argv[5]);
	string name = argv[6];
	int chronoprior = 0;
	string prior = argv[7];
	if (prior == "bd")	{
		chronoprior = 1;
	}
	else if (prior == "uni")	{
		chronoprior = 0;
	}
	else	{
		cerr << "error : does not recognise " << argv[7] << '\n';
		exit(1);
	}

	if (chronoprior)	{
		cerr << "birth death prior\n";
	}
	else	{
		cerr << "uniform prior\n";
	}

	GTRLogNormalModel* model = new GTRLogNormalModel(datafile,treefile,calibfile,rootage,rootstdev,chronoprior);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Trace(os);
		os.flush();
		ofstream mos((name + ".monitor").c_str());
		model->Monitor(mos);
	}
}
