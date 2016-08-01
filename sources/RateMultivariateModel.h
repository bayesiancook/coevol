
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"

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


class BranchMatrixPhyloProcess : public PhyloProcess	{


	protected:

	public:

	BranchMatrixPhyloProcess(LengthTree* intree, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:
	BranchValPtrTree<RandomSubMatrix>* matrixtree;
};


class RateMultivariateModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
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
	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	int chronoprior;
	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;

	// chronogram

	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;
	Gamma* phi;
	IIDUniform* mean;

	MultiVariateTreeProcess* process;

	MeanExpTreeFromMultiVariate* ratetree;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;

	// for homogeneous model
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	// for gc model
	MeanLogitTreeFromMultiVariate* gctree;
	Beta* rootgc;
	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;

	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	// if true: OUP else BP
	bool autoregressive;

	bool clamproot;

	bool meanexp;

	// if true: gc model
	bool gc;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	bool conjpath;

	public:

	RateMultivariateModel() {
	}

	RateMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, bool ingc, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, bool inclamproot, bool inmeanexp, bool sample=true)	{

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		gc = ingc;
		L = gc ? 2 : 1;

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);
		Nsite = nucdata->GetNsite();	// # columns
		Nstate = nucdata->GetNstate();	// # states (20 for amino acids)

		if (inconjpath == -1)	{
			/*
			if (Nsite > 5000)	{
				conjpath = true;
			}
			else	{
				conjpath = false;
			}
			*/
			conjpath = true;
		}
		else	{
			conjpath = inconjpath;
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

		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu);
		mu->ClampAt(1);

		RootAlpha = 0;
		RootBeta = 0;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		if (calibfile != "None")	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			cerr << "gamma params : " << a << '\t' << b << '\n';
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,mu,RootAlpha,RootBeta,calibset);
			}
			else if (chronoprior == 1)	{
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,RootAlpha,RootBeta,calibset);
			}
			else	{
				cerr << "error : chronoprior\n";
				exit(1);
			}
		}
		else	{
			chronogram = new Chronogram(tree,mu);
		}

		GammaDiag = new Gamma*[Ncont+L];
		diag = new Rvar<PosReal>*[Ncont+L];
		for (int k=0; k<Ncont+L; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(priorsigma);
		}
		priorOnSigmaZero = new RvarVec(diag,Ncont+L);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		if (clampdiag)	{
			sigma = new DiagonalCovMatrix(sigmaZero, Ncont+L+1);
		}
		else	{
			sigma = new InverseWishartMatrix(sigmaZero, Ncont+L+1);
		}
		sigma->SetIdentity();

		if (autoregressive)	{
			phi = new Gamma(One,One);
			mean = new IIDUniform(Ncont+L,100);
			process = new AutoRegressiveMultiVariateTreeProcess(sigma,mean,phi,chronogram);
		}
		else	{
			phi = 0;
			mean = 0;
			process = new MultiVariateTreeProcess(sigma,chronogram);
		}

		// process->ClampTransversal(1,0);

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,L+i,i,contdatatype);
			}
		}

		if (clamproot)	{
			process->ClampRoot();
		}

		// cut off to avoid numerical errors
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		ratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);

		// make rate matrices
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (gc)	{
			gctree = new MeanLogitTreeFromMultiVariate(process,1,MEAN,false);
			rootgc = new Beta(One,One);
			stattree = new GCStatTree(gctree,rootgc);
			nucmatrixtree = new NucMatrixTree(relrate,stattree);
			nucmatrix = 0;
			stationary = 0;
		}
		else	{
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);
			gctree = 0;
			rootgc = 0;
			stattree = 0;
			nucmatrixtree = 0;
		}

		// make substitution mappings
		if (conjpath)	{
			if (gc)	{
				pathconjtree = new BranchMatrixPathConjugateTree(ratetree, nucmatrixtree, nucdata);
			}
			else	{
				pathconjtree = new OneMatrixPathConjugateTree(ratetree,nucmatrix,nucdata);
			}
			phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
		}
		else	{
			pathconjtree = 0;
			if (gc)	{
				phyloprocess = new BranchMatrixPhyloProcess(ratetree, nucmatrixtree, nucdata);
			}
			else	{
				phyloprocess = new OneMatrixPhyloProcess(ratetree, nucmatrix, nucdata);
			}
		}

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		// register model
		RootRegister(PriorMu);
		RootRegister(One);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		if (chronoprior == 1)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (autoregressive)	{
			RootRegister(Zero);
		}
		RootRegister(relrate);
		if (!gc)	{
			RootRegister(stationary);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();

			if (RootAlpha)	{
				cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
				cerr << "scale progeny : " << GetCalibratedChronogram()->GetScale()->down.size() << '\n';
			}
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~RateMultivariateModel() {}

	Tree* GetTree() {return tree;}

	MeanExpTreeFromMultiVariate* GetRateTree() {return ratetree;}
	MeanLogitTreeFromMultiVariate* GetGCTree() {return gctree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	Chronogram* GetChronogram() {return chronogram;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double Update(bool check = false)	{
		double ret = ProbModel::Update();
		phyloprocess->Sample();
		ret = ProbModel::Update();
		return ret;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if (chronoprior == 1)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();
		for (int k=0; k<Ncont + L; k++)	{
			total += GammaDiag[k]->GetLogProb();
		}

		total += sigma->GetLogProb();
		if (autoregressive)	{
			total += phi->GetLogProb();
			total += mean->GetLogProb();
		}
		total += process->GetLogProb();

		total += relrate->GetLogProb();

		if (gc)	{
			total += rootgc->GetLogProb();
		}
		else	{
			total += stationary->GetLogProb();
		}
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		int nrep = conjpath ? 30 : 1;
		for (int i=0; i<nrep; i++)	{
			if (chronoprior == 1)	{
				scheduler.Register(new SimpleMove(Chi,1),10,"chrono");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(Chi2,1),10,"chrono");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"chrono");
			}
			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			if (RootAlpha)	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
			}

			/*
			for (int k=0; k<Ncont+L; k++)	{
				scheduler.Register(new SimpleMove(GammaDiag[k],10),10,"theta");
				scheduler.Register(new SimpleMove(GammaDiag[k],1),10,"theta");
				scheduler.Register(new SimpleMove(GammaDiag[k],0.1),10,"theta");
			}
			*/

			scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");

			scheduler.Register(new SimpleMove(process,10),40,"multinormal");
			scheduler.Register(new SimpleMove(process,1),40,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),40,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),40,"multinormal");

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

				scheduler.Register(new SimpleMove(mean,1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
			}

			/*
			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,1,0,1),4,"translation multinormal");
			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,0.1,0,1),4,"translation multinormal");
			scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,0.01,0,1),4,"translation multinormal");
			*/

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			if (gc)	{
				scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");
			}
			else	{
				scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
			}

			/*
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,1,0),5,"root compensation");
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.1,0),5,"root compensation");
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.01,0),5,"root compensation");
			*/

		}

		// scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";

		if (chronoprior == 1)	{
			Chi->Sample();
			Chi2->Sample();
		}
		chronogram->Sample();
		cerr << "gamma diag\n";
		for (int k=0; k<Ncont+L; k++)	{
			GammaDiag[k]->Sample();
		}

		sigma->Sample();
		if (autoregressive)	{
			phi->Sample();
			mean->Sample();
		}

		process->Sample();

		process->CutOff(2,0);
		ratetree->specialUpdate();

		relrate->Sample();
		if (gc)	{
			rootgc->Sample();
		}
		else	{
			stationary->Sample();
		}

		phyloprocess->Sample();

		cerr << "ok\n";
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	double GetMaxRate()	{
		return ratetree->GetMax();
	}

	double GetLength()	{
		return ratetree->GetTotal();
	}

	double GetMeanGC()	{
		if (! gc)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree->GetMean();
	}

	double GetRootGC()	{
		if (! gc)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return rootgc->val();
	}

	Var<RealVector>* GetMeanVector()	{
		return mean;
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength";

		// os << "\trootleft\trootright";

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}

		if (RootAlpha)	{
			os << "\trootage";
		}
		if (chronoprior == 1)	{
			os << "\tp1\tp2";
		}

		if (autoregressive)	{
			os << '\t' << "theta";
			os << '\t' << "dim";
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << "mean" << k;
			}
		}

		if (gc)	{
			os << "\tmeangc\trootgc";
		}
		else	{
			os << "\tstatent";
		}
		os << "\trrent";

		os << '\n';
	}

	void PrintEntries(ostream& os)	{

		os << "rate\n";
		if (gc == 3)	{
			os << "gc1\n";
			os << "gc2\n";
			os << "gc3\n";
		}
		else if (gc == 1)	{
			os << "gc\n";
		}
		for (int k=0; k<Ncont; k++)	{
			os << "character " << k + 1 << '\n';
		}
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();

		/*
		os << '\t' << *chronogram->GetBranchVal(tree->GetRoot()->Next()->GetBranch());
		os << '\t' << *chronogram->GetBranchVal(tree->GetRoot()->Next()->Next()->GetBranch());
		*/

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}

		if (RootAlpha)	{
			os << '\t' << GetRootAge();
		}
		if (chronoprior == 1)	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}

		if (autoregressive)	{
			os << '\t' << phi->val();
			os << '\t' << *mean;
		}

		if (gc)	{
			os << '\t' << GetMeanGC();
			os << '\t' << GetRootGC();
		}
		else	{
			os << '\t' << stationary->val().GetEntropy();
		}

		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *chronogram << '\n';
		if (RootAlpha)	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		if (chronoprior == 1)	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << *GammaDiag[k] << '\n';
		}
		os << *sigma << '\n';
		if (autoregressive)	{
			os << *phi << '\n';
			os << *mean << '\n';
		}
		os << '\n';
		os << *process << '\n';
		os << *relrate << '\n';
		if (gc)	{
			os << *rootgc << '\n';
		}
		else	{
			os << *stationary << '\n';
		}
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		if (RootAlpha)	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		if (chronoprior == 1)	{
			is >> *Chi >> *Chi2;
		}
		for (int k=0; k<Ncont+L; k++)	{
			is >> *GammaDiag[k];
		}
		is >> *sigma;
		if (autoregressive)	{
			is >> *phi;
			is >> *mean;
		}
		is >> *process;
		is >> *relrate;
		if (gc)	{
			is >> *rootgc;
		}
		else	{
			is >> *stationary;
		}
	}
};

#endif
