
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
	SequenceAlignment* data1;
	SequenceAlignment* data2;
	TaxonSet* taxonset;

	// number of columns
	int Nsite1;
	int Nsite2;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate1;
	int Nstate2;

	ContinuousData* contdata;
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

	MeanExpTreeFromMultiVariate* synratetree1;
	MeanExpTreeFromMultiVariate* synratetree2;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate1;
	Dirichlet* stationary1;
	GTRRandomSubMatrixWithNormRates* nucmatrix1;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate2;
	Dirichlet* stationary2;
	GTRRandomSubMatrixWithNormRates* nucmatrix2;

	// phylo process
	PathConjugateTree* pathconjtree1;
	PhyloProcess* phyloprocess1;
	PathConjugateTree* pathconjtree2;
	PhyloProcess* phyloprocess2;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamptree;
	bool meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int nrep;

	int df;

	bool iscalib;

	double** cov;

	public:


	BranchOmegaMultivariateModel(string datafile1, string datafile2, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, double priorsigma, int indf, int contdatatype, bool inclamptree, bool inmeanexp, int innrep, bool sample=true)	{

		clamptree = inclamptree;
		meanexp = inmeanexp;
		L = 2;

		// get data from file
		data1 = new FileSequenceAlignment(datafile1);
		Nsite1 = data1->GetNsite();
		Nstate1 = data1->GetNstate();

		// get data from file
		data2 = new FileSequenceAlignment(datafile2);
		Nsite2 = data2->GetNsite();
		Nstate2 = data2->GetNstate();

		nrep = innrep;
		if (nrep == 0)	{
			nrep = 30;
		}

		taxonset = data1->GetTaxonSet();

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

		df = Ncont + L + indf;
		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		cerr << "new chrono\n";
		if (calibfile != "None")	{
			iscalib = true;
			cerr << "calibrated chronogram deactivated\n";
			exit(1);
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
		synratetree1 = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
		synratetree2 = new MeanExpTreeFromMultiVariate(process,1,INTEGRAL,false,meanexp);

		// create a GTR nucleotide matrix
		relrate1 = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary1 = new Dirichlet(Nnuc);
		nucmatrix1 = new GTRRandomSubMatrixWithNormRates(relrate1,stationary1,true);

		// create a GTR nucleotide matrix
		relrate2 = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary2 = new Dirichlet(Nnuc);
		nucmatrix2 = new GTRRandomSubMatrixWithNormRates(relrate2,stationary2,true);

		pathconjtree1 = new OneMatrixPathConjugateTree(synratetree1,nucmatrix1,data1);
		phyloprocess1 = new PathConjugatePhyloProcess(pathconjtree1);

		pathconjtree2 = new OneMatrixPathConjugateTree(synratetree2,nucmatrix2,data2);
		phyloprocess2 = new PathConjugatePhyloProcess(pathconjtree2);

		cerr << "unfold\n";
		phyloprocess1->Unfold();
		phyloprocess2->Unfold();
		cerr << "sample\n";
		if (sample)	{
			if (phyloprocess1)	{
				phyloprocess1->Sample();
				phyloprocess2->Sample();
			}
		}

		cerr << "register\n";
		// register model
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(relrate1);
		RootRegister(stationary1);
		RootRegister(relrate2);
		RootRegister(stationary2);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}

		cov = new double*[L + Ncont];
		for (int i=0; i<L+Ncont; i++)	{
			cov[i] = new double[L+Ncont];
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BranchOmegaMultivariateModel() {}


	// accessors
	Tree* GetTree() {return tree;}

	MeanExpTreeFromMultiVariate* GetSynRateTree1() {return synratetree1;}
	MeanExpTreeFromMultiVariate* GetSynRateTree2() {return synratetree2;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	Chronogram* GetChronogram() {return chronogram;}

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

		total += relrate1->GetLogProb();
		total += stationary1->GetLogProb();
		total += relrate2->GetLogProb();
		total += stationary2->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		ret += phyloprocess1->GetLogProb();
		ret += phyloprocess2->GetLogProb();
		/*
		ret += pathconjtree1->GetLogProb();
		ret += pathconjtree2->GetLogProb();
		*/
		return ret;
	}

	// MCMC schedule
	virtual void MakeScheduler()	{

		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess1,pathconjtree1),1,"mapping + sufficient stat");
		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess2,pathconjtree2),1,"mapping + sufficient stat");

		for (int i=0; i<nrep; i++)	{

			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.01,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.001,10),1,"conjugate sigma - process");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			scheduler.Register(new ProfileMove(relrate1,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate1,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate1,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary1,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary1,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary1,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary1,0.001),10,"stat");

			scheduler.Register(new ProfileMove(relrate2,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate2,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate2,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary2,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary2,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary2,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary2,0.001),10,"stat");
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
		relrate1->Sample();
		stationary1->Sample();
		phyloprocess1->Sample();
		relrate2->Sample();
		stationary2->Sample();
		phyloprocess2->Sample();
	}

	// summary statistics
	double GetTotalLength1()	{
		return GetSynRateTree1()->GetTotal();
	}

	double GetTotalLength2()	{
		return GetSynRateTree2()->GetTotal();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tdS1\tds2";

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "cov_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "cov_" << k << '_' << k;
		}
		/*
		os << "\tstatent1";
		os << "\tstatent2";
		os << "\trrent1";
		os << "\trrent2";
		*/
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetTotalLength1();
		os << '\t' << GetTotalLength2();
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		process->GetEmpiricalCovariance(cov);
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << cov[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << cov[k][k];
		}
		/*
		os << '\t' << stationary1->val().GetEntropy();
		os << '\t' << stationary2->val().GetEntropy();
		os << '\t' << relrate1->val().GetEntropy();
		os << '\t' << relrate2->val().GetEntropy();
		*/
		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *relrate1 << '\n';
		os << *stationary1 << '\n';
		os << *relrate2 << '\n';
		os << *stationary2 << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *relrate1;
		is >> *stationary1;
		is >> *relrate2;
		is >> *stationary2;
	}
};

#endif
