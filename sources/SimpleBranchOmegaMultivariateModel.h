
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "CalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "GeneralConjugatePath.h"
#include "Jeffreys.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "MeanChronogram.h"
#include "AuxCoevol.h"


class BranchOmegaMultivariateModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
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

	Const<Real>* Zero;
	Const<PosReal>* One;

	Chronogram* chronogram;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* sigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	ConjugateMultiVariateTreeProcess* process;

	MeanExpTreeFromMultiVariate* synratetree;
	MeanExpTreeFromMultiVariate* omegatree;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamptree;
	bool meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int nrep;

	int df;

	bool iscalib;

	public:


	BranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, double priorsigma, int indf, int contdatatype, bool inclamptree, bool inmeanexp, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		// if we want the divergence times to be fixed
		clamptree = inclamptree;

		// arithmetic / geodesic averages (see Poujol et al)
		meanexp = inmeanexp;

		// number of components of the Brownian process corresponding to rates
		// here L = 2: dS and dN/dS
		L = 2;

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		// make codon alignment 
		codondata = new CodonSequenceAlignment(nucdata, true, type);

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		// number of MCMC cycles per point
		nrep = innrep;
		if (nrep == 0)	{
			nrep = 30;
		}

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		cerr << "tree and data ok\n";
		cerr << '\n';

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		cerr << "new chrono\n";
		if (calibfile != "None")	{
			iscalib = true;
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);
			chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
		}
		else	{
			iscalib = false;
			chronogram = new Chronogram(tree,One);
		}
		if (clamptree)	{
			chronogram->Clamp();
		}
		cerr << "ok\n";

		// Ncont : number of quantitative traits
		// L : number of substitution parameters coevolving with traits (typically, 2: dS and dN/dS).
		// create an array of positive variables kappa_i, i=1..Ncont + L
		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		// create a diagonal matrix, with the kappa_i along the diagonal
		sigmaZero = new SigmaZero(DiagArray);

		// create covariance matrix
		// from an inverse wishart of parameter sigma0
		// and df degrees of freedom
		df = Ncont + L + indf;
		cerr << "sigma\n";
		sigma = new ConjugateInverseWishart(sigmaZero, df);

		// create a multivariate brownian process (of dimension Ncont + L)
		// along the chronogram, and with covariance matrix sigma
		cerr << "process\n";
		process = new ConjugateMultiVariateTreeProcess(sigma,chronogram);

		// condition the multivariate process
		// on the matrix of quantitative traits.
		// note the offset here : first trait corresponds to entry L+1 of the process, etc.
		// this is because the first L entries of the process correspond to the substitution variables (dS, dN/dS)
		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,L+i,i,contdatatype);
		}

		// just for numerical stability of the starting point
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		// create the branch lengths resulting from combining
		// the times given by the chronogram with the rate (first entry of the multivariate process)
		cerr << "syn and omega\n";
		synratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);

		// create the dN/dS on each branch, nased on the second entry of the multivariate process
		omegatree = new MeanExpTreeFromMultiVariate(process,1,MEAN,false,meanexp);
		cerr << "matrix\n";

		// create a GTR nucleotide matrix
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);

		// create a series of codon matrices, one for each branch of the chronogram
		// each codon matrix is parameterized by the global nucleotide matrix, combined with the branch-specific dN/dS
		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

		// create a phylo process based on this array of branch specific matrices
		// and condition it on the multiple alignment codondata
		pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);


		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "sample\n";
		if (sample)	{
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		cerr << "register\n";
		// register model
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BranchOmegaMultivariateModel() {}


	// accessors
	Tree* GetTree() {return tree;}

	SequenceAlignment* GetData()	{
		return codondata;
	}

	MeanExpTreeFromMultiVariate* GetSynRateTree() {return synratetree;}
	MeanExpTreeFromMultiVariate* GetOmegaTree() {return omegatree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}

	Chronogram* GetChronogram() {return chronogram;}

	bool isCalibrated()	{
		return iscalib;
	}

	CalibratedChronogram* GetCalibratedChronogram()	{
		if (! iscalib)	{
			cerr << "error : calibrated chronogram does not exist under uncalibrated model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	Var<PosReal>* GetScale()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	// probability computation

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += chronogram->GetLogProb();
		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = pathconjtree->GetLogProb();
		return ret;
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	// MCMC schedule
	virtual void MakeScheduler()	{

		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<nrep; i++)	{

			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			if (isCalibrated())	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");
			}

			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.01,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.001,10),1,"conjugate sigma - process");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"kappa");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"kappa");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"kappa");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}
	}

	void drawSample()	{
		chronogram->Sample();
		DiagArray->Sample();
		sigma->Sample();
		// sigma->SetIdentity();
		process->Sample();
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		relrate->Sample();
		stationary->Sample();
		phyloprocess->Sample();
	}

	// summary statistics
	double GetTotalLength()	{
		return GetSynRateTree()->GetTotal();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tdS\tomega";
		if (isCalibrated())	{
			os << "\trootage";
		}

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		os << "\tstatent";
		os << "\trrent";
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetTotalLength();
		os << '\t' << GetMeanOmega();
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *relrate;
		is >> *stationary;
	}
};

#endif