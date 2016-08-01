
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "GeneralConjugatePath.h"
#include "Jeffreys.h"
#include "PrecisionNormalTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "MeanChronogram.h"
#include "SumConstrained.h"

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(Var<Profile>* inrelrate, BranchVarTree<Profile>* instattree, Var<Profile>* inrootstat) {
		SetWithRoot(true);
		stattree = instattree;
		relrate = inrelrate;
		rootstat = inrootstat;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new GTRRandomSubMatrixWithNormRates(relrate,rootstat,false);
		}
		return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),true);
		// return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),false);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;
	Var<Profile>* relrate;
	Var<Profile>* rootstat;

};

class CompoModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
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

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	int chronoprior;
	double softa;

	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;
	// 2 : bd with cauchy proper lower bounds

	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;
	GammaTree* syngammatree;
	LengthTree* lengthtree;
	LengthTree* synratetree;
	MeanExpTreeFromMultiVariate* meanexpsynratetree;
	Jeffreys* synsigma;
	LogNormalTreeProcess* lognormalsyntree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	// ConjugateInverseWishart* sigma;
	Rvar<CovMatrix>* sigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* process;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;

	Dirichlet* rootstationary;

	SumConstrainedMapping* changeofbasis;
	SumConstrainedStatTree* stattree;

	MatrixTree* matrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	int clampdiag;
	int clamptree;
	int meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int nrep;

	int L;
	int df;
	int jeffdim;

	int Ksyn;
	int Kcomp;
	int Kcont;

	bool iscalib;

	public:

	ConjugateInverseWishart* GetConjugateInverseWishart() {
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigma);
		if (! tmp)	{
			cerr << "error : dynamic castt of conjugate inverse wishart : " << sigma << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConjugateMultiVariateTreeProcess* GetConjugateMultiVariateTreeProcess() {
		ConjugateMultiVariateTreeProcess* tmp = dynamic_cast<ConjugateMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}


	CompoModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double insofta,  double inmeanchi, double inmeanchi2, int indf, int inclampdiag, int contdatatype, string rrtype, int inclamptree, int inmeanexp, int innormalise, int innrep, string inmix, string rootfile, int inseparatesyn, bool sample=true)	{

		clamptree = inclamptree;
		meanexp = inmeanexp;
		clampdiag = inclampdiag;

		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		cerr << "ali\n";
		// get data from file
		data = new FileSequenceAlignment(datafile);

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		nrep = innrep;
		if (nrep == 0)	{
			nrep = 30;
		}

		taxonset = data->GetTaxonSet();

		cerr << "tree\n";
		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		cerr << "cont data\n";
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

		double mindiag = 0.001;
		double maxdiag = 1000;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		iscalib = false;

		syngammatree = 0;
		lognormalsyntree = 0;
		synratetree = 0;

		L = Nstate;

		if (calibfile == "Unconstrained")	{
			L--;
			PriorMu = new Const<PosReal>(0.1);
			mu = new Gamma(One,PriorMu);
			syngammatree = new GammaTree(tree,One,mu);
			lengthtree = syngammatree;
			synratetree = syngammatree;
		}
		else	{
			PriorMu = new Const<PosReal>(1);
			mu = new Gamma(One,PriorMu);
			mu->ClampAt(1);

			if (calibfile != "None")	{
				iscalib = true;

				double a = rootage * rootage / rootstdev / rootstdev;
				double b = rootage / rootstdev / rootstdev;
				if (rootage == -1)	{
					a = b = -1;
				}
				CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

				if (chronoprior == 0)	{
					chronogram = new CalibratedChronogram(tree,mu,a,b,calibset);
				}
				else {
					MeanChi = new Const<PosReal>(meanchi);
					MeanChi2 = new Const<PosReal>(meanchi2);
					Chi = new Exponential(MeanChi,Exponential::MEAN);
					Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
					chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,a,b,calibset,chronoprior,softa);
				}
			}
			else	{
				chronogram = new Chronogram(tree,mu);
			}

			if (clamptree)	{
				chronogram->Clamp();
			}

			lengthtree = chronogram;

			if (inseparatesyn)	{
				synsigma = new Jeffreys(mindiag,maxdiag,Zero);
				synsigma->setval(1);
				lognormalsyntree = new LogNormalTreeProcess(lengthtree,synsigma,INTEGRAL);
				synratetree = lognormalsyntree;
				L--;
			}
		}

		/*
		if (mix == "branch")	{
			MixAlpha = new Jeffreys(0.1,10);
			MixAlpha->setval(1.0);
			gammatree = new GammaTree(lengthtree->GetTree(),MixAlpha,MixAlpha);
			// gammatree->GetBranchVal(gammatree->GetRoot()->Next()->Out()->GetBranch())->ClampAt(1.0);
			partition = new BranchPartition(tree);
			gammamixtree = 0;
		}
		else if (mix != "None")	{
			partition = new BranchPartition(tree,mix);
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
			// gammamixtree = new GammaMixTree(partition,One,One,true); // clamp first component to 1
			gammatree = 0;
		}
		else	{
			partition = new BranchPartition(tree);
			gammamixtree = 0;
			gammatree = 0;
		}

		Nmat = GetPartition()->GetNComponent();
		*/

		if (L == Nstate)	{
			Ksyn = 0;
			Kcomp = 1;
			Kcont = Nstate;
		}
		else	{
			Ksyn = -1;
			Kcomp = 0;
			Kcont = Nstate-1;
		}

		cerr << "Kcont : " << Kcont << '\n';

		df = Ncont + L + indf;

		// create an array of positive variables kappa_i, i=1..Ncont + L
		jeffdim = Ncont + L - Nstate + 2;
		DiagArray = new JeffreysIIDArray(jeffdim,mindiag,maxdiag,Zero);
		for (int i=0; i<jeffdim; i++)	{
			// DiagArray->GetVal(i)->ClampAt(1.0);
			DiagArray->GetVal(i)->setval(1.0 + 0.01 * (Random::Uniform() - 0.5));
		}
		Var<PosReal>** diagarray = new Var<PosReal>*[Ncont + L];
		int index = 0;
		int k = 0;
		if (Ksyn != -1)	{
			diagarray[0] = DiagArray->GetVal(0);
			index++;
			k++;
		}
		for (int i=0; i<Nstate-1; i++)	{
			diagarray[index] = DiagArray->GetVal(k);
			index++;
		}
		k++;
		for (int i=0; i<Ncont; i++)	{
			diagarray[index] = DiagArray->GetVal(k);
			index++;
			k++;
		}
		if (index != (Ncont + L))	{
			cerr << "error when constructing diag array: index\n";
			exit(1);
		}
		if (k != jeffdim)	{
			cerr << "error when constructing diag array: k\n";
			exit(1);
		}

		cerr << L << '\t' << Ncont << '\t' << jeffdim << '\t' << Nstate << '\n';
		cerr << Ksyn << '\t' << Kcomp << '\t' << Kcont << '\n';

		// create a diagonal matrix, with the kappa_i along the diagonal
		sigmaZero = new SigmaZero(diagarray, Ncont+L);

		// create covariance matrix
		// from an inverse wishart of parameter sigma0
		cerr << "sigma\n";
		if (clampdiag)	{
			sigma = new DiagonalCovMatrix(sigmaZero, df);
		}
		else	{
			sigma = new ConjugateInverseWishart(sigmaZero, df);
		}
		// sigma->SetIdentity();

		// create a multivariate brownian process (of dimension Ncont + L)
		// along the chronogram, and with covariance matrix sigma
		cerr << "process\n";

		if (clampdiag)	{
			process = new MultiVariateTreeProcess(sigma,lengthtree);
		}
		else	{
			process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),lengthtree);
		}

		// process->Reset();

		// condition the multivariate process
		// on the matrix of quantitative traits.
		// note the offset here : first trait corresponds to entry L+1 of the process, etc.
		// this is because the first L entries of the process correspond to the substitution variables (dS, dN/dS)
		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,Kcont+i,i,contdatatype);
			(*process->GetRootVal())[Kcont+i] = contdata->GetMeanLog(i);
		}

		// just for numerical stability of the starting point
		for (int i=0; i<L; i++)	{
			process->CutOff(1,i);
		}

		if (Ksyn != -1)	{
			// create the branch lengths resulting from combining
			// the times given by the chronogram with the rate (first entry of the multivariate process)
			cerr << "syn tree\n";
			meanexpsynratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
			synratetree = meanexpsynratetree;
		}

		cerr << "matrix\n";

		// create a GTR nucleotide matrix
		relrate = new Dirichlet(Nstate*(Nstate-1)/2);
		if ((Nstate == Naa) && (rrtype == "lg"))	{
			double total = 0;
			for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
				(*relrate)[i] = LG_RR[i];
				total += (*relrate)[i];
			}
			for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
				(*relrate)[i] /= total;
			}
			relrate->Clamp();
		}

		rootstationary = new Dirichlet(Nstate);

		changeofbasis = new SumConstrainedMapping(Nstate);

		/*
		cerr << '\n';
		changeofbasis->ToStream(cerr);
		cerr << '\n';
		*/

		cerr << Kcomp << '\n';
		stattree = new SumConstrainedStatTree(process,Kcomp,changeofbasis);

		matrixtree = new MatrixTree(relrate,stattree,rootstationary);

		// create a phylo process based on this array of branch specific matrices
		// and condition it on the multiple alignment codondata
		pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, data);
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
		RootRegister(PriorMu);
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		RootRegister(relrate);
		RootRegister(rootstationary);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~CompoModel() {}

	bool Unconstrained()	{
		return syngammatree;
	}

	bool SeparateSyn()	{
		return lognormalsyntree;
	}

	bool isCalibrated()	{
		return iscalib;
	}

	// accessors
	Tree* GetTree() {return tree;}

	SequenceAlignment* GetData()	{
		return data;
	}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}

	MeanExpTreeFromMultiVariate* GetSynRateTree() {return meanexpsynratetree;}

	LengthTree* GetLengthTree() {return lengthtree;}

	Chronogram* GetChronogram() {
		if (Unconstrained())	{
			cerr << "error : chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return chronogram;
	}
	CalibratedChronogram* GetCalibratedChronogram()	{
		if (Unconstrained())	{
			cerr << "error : calibrated chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
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

		if (Unconstrained())	{
			total += mu->GetLogProb();
			total += syngammatree->GetLogProb();
		}
		else	{
			if (chronoprior)	{
				total += Chi->GetLogProb();
				total += Chi2->GetLogProb();
			}
			total += chronogram->GetLogProb();
		}

		if (SeparateSyn())	{
			total += synsigma->GetLogProb();
			total += lognormalsyntree->GetLogProb();
		}

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();
		total += rootstationary->GetLogProb();

		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		// double ret = pathconjtree->GetLogProb();
		return ret;
	}

	// MCMC schedule
	virtual void MakeScheduler()	{

		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<nrep; i++)	{

			if (Unconstrained())	{
				scheduler.Register(new SimpleMove(mu,1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(mu,0.1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(syngammatree,1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.01),10,"syngamtree");
			}
			else if (! clamptree)	{
				if (chronoprior)	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.001),10,"chrono");
			}

			if (isCalibrated())	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");
			}

			if (SeparateSyn())	{
				scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");
				scheduler.Register(new SimpleMove(lognormalsyntree,10),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,1),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.1),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.01),10,"log normal syn");
			}

			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.3),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.03),10,"multinormal");

			if (! clampdiag)	{
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,1),1,"conjugate sigma - process");
			}
			else	{
				scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
			}

			/*
			if (! clampdiag)	{
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.3,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.03,10),1,"conjugate sigma - process");
			}
			else	{
				scheduler.Register(new SimpleMove(process,1),10,"multinormal");
				scheduler.Register(new SimpleMove(process,0.3),10,"multinormal");
				scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
				scheduler.Register(new SimpleMove(process,0.03),10,"multinormal");

				scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
				scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
			}
			*/

			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,10,Kcont,1),10,"conjugate sigma - process - piecewise translation");
			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,1,Kcont,1),10,"conjugate sigma - process - piecewise translation");
			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,0.1,Kcont,1),10,"conjugate sigma - process - piecewise translation");
			// scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(sigma,process,1,20,1,30),1,"conjugate sigma - process - piecewise translation");
			// scheduler.Register(new MultiVariatePiecewiseTranslationMove(sigma,process,0.1,20,1,30),1,"conjugate sigma - process - piecewise translation");


			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(rootstationary,0.1,2),10,"stat");
			scheduler.Register(new ProfileMove(rootstationary,0.03,2),10,"stat");
			scheduler.Register(new ProfileMove(rootstationary,0.03,5),10,"stat");
			scheduler.Register(new SimpleMove(rootstationary,0.01),10,"stat");
		}
	}

	void drawSample()	{
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	// summary statistics

	double GetMeanSynRate()	{
		if (Unconstrained())	{
			return syngammatree->GetMean();
		}
		if (SeparateSyn())	{
			return lognormalsyntree->GetMeanRate();
		}
		return GetSynRateTree()->GetTotal();
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

	double GetMeanEntropy()	{
		return stattree->GetMeanEntropy();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\trate\tstatent";
		for (int k=0; k<Ncont; k++)	{
			os << "\trootcont" << k;
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		os << "\tlogdet";
		/*
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		*/
		os << "\tstatent";
		os << "\trrent";
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		os << '\t' << GetMeanEntropy();
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << (*process->GetRootVal())[Kcont + k];
		}
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		os << '\t' << sigma->GetLogDeterminant();
		/*
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		*/
		os << '\t' << rootstationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *mu << '\n';
		if (Unconstrained())	{
			os << *syngammatree << '\n';
		}
		else	{
			os << *chronogram << '\n';
			if (isCalibrated())	{
				os << *GetCalibratedChronogram()->GetScale() << '\n';
			}
			if (chronoprior)	{
				os << *Chi << '\t' << *Chi2 << '\n';
			}
		}
		if (SeparateSyn())	{
			os << *synsigma << '\n';
			os << *lognormalsyntree << '\n';
		}
		changeofbasis->ToStream(os);
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *relrate << '\n';
		os << *rootstationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		if (Unconstrained())	{
			is >> *syngammatree;
		}
		else	{
			is >> *chronogram;
			if (isCalibrated())	{
				is >> *GetCalibratedChronogram()->GetScale();
			}
			if (chronoprior)	{
				is >> *Chi >> *Chi2;
			}
		}
		if (SeparateSyn())	{
			is >> *synsigma;
			is >> *lognormalsyntree;
		}
		changeofbasis->FromStream(is);
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *relrate;
		is >> *rootstationary;
	}
};

#endif
