
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"

#include "BDChronogram.h"
#include "BDCalibratedChronogram.h"

#include "BranchProcess.h"
#include "MatrixTree.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"

#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"


#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"

#include "Jeffreys.h"

#include "Partition.h"
#include "SplitPartition.h"

#include "SplitLengthTree.h"

#include "SplitMultiVariateMove.h"
#include "MultiVariatePropagateMove.h"

#include "PartitionMultiVariateTreeProcess.h"

#include "WhiteNoise.h"
#include "TimeLineMultiVariateTreeProcess.h"
#include "MeanChronogram.h"

#include "WishartArray.h"

#include "ProteinSequenceAlignment.h"
#include "AminoAcidOmegaSubMatrix.h"
#include "SimilarityMatrix.h"
#include "AminoAcidMatrixTree.h"

#include "PosRealVectorMove.h"
#include "RenormalizedPosRealVector.h"
#include "MultiVariateRelRateCompensatoryMove.h"
#include "MultiVariateUgamCompMove.h"
#include "ChronoCompensatoryMove.h"

class BranchOmegaMultivariateModel : public ProbModel {

	public:

	// data fields

	double lnL;

	int nprocs;
	int myid;
	int* sitemin;
	int* sitemax;

	int sample;

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
	// number of states (4 for nucleic acids, 20 for amino-acids. !!! 4 for codons)
	int rawNstate;

	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Rvar<PosReal>* Chi;
	Rvar<PosReal>* Chi2;

	int Ninterpol;

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
	bool iscalib;

	GammaTree* syngammatree;

	LengthTree* lengthtree;

	// if different branches have different scaling factors
	BranchPartition* partition;
	SplitBranchPartition* splitpartition;
	GammaMixTree* gammamixtree;
	GammaTree* gammatree;
	Jeffreys* MixAlpha;

	JeffreysIIDArray* TimeLineDiagArray;
	SigmaZero* TimeLineSigmaZero;
	Rvar<CovMatrix>* timelinesigma;

	TimeIntervals* timeintervals;
	TimeLine* timeline;

	string mix;

	int Nmat;

	JeffreysIIDArray* DiagArray;
	JeffreysIIDArray* DiagArray0;
	SigmaZero* sigmaZero;
	SigmaZero* sigmaZero0;
	IIDArray<CovMatrix>* sigmaarray;
	Rvar<CovMatrix>* sigma;

	RandomVarArray<RealVector>* driftarray;
	Var<RealVector>* drift;

	RandomVarArray<PosReal>* driftphiarray;
	Var<PosReal>* driftphi;

	RandomVarArray<RealVector>* driftarray2;
	Var<RealVector>* drift2;

	RandomVarArray<PosReal>* driftphiarray2;
	Var<PosReal>* driftphi2;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	Gamma* phi;
	IIDUniform* mean;

	Jeffreys* synsigma;
	Jeffreys* omegasigma;
	LogNormalTreeProcess* lognormalsyntree;
	LogNormalTreeProcess* lognormalomegatree;

	MultiVariateTreeProcess* process;

	LengthTree* synratetree;
	// MeanExpTreeFromMultiVariate* synratetree;
	MeanExpTreeFromMultiVariate* nonsynratetree;
	MeanExpTreeFromMultiVariate* kappatree;
	MeanExpTreeFromMultiVariate* kappatree1;
	MeanExpTreeFromMultiVariate* kappatree2;
	MeanExpTreeFromMultiVariate* kappatree3;
	BranchValPtrTree<Dvar<PosReal> >* omegatree;
	BranchValPtrTree<Dvar<PosReal> >* omegatvgctree;
	BranchValPtrTree<Dvar<PosReal> >* omegatv0tree;
	BranchValPtrTree<Dvar<PosReal> >* omegatstree;
	MeanExpTreeFromMultiVariate* synratetvgctree;
	MeanExpTreeFromMultiVariate* synratetstree;

	MeanExpTreeFromMultiVariate* tstvtree;
	MeanExpTreeFromMultiVariate* tvgctree;

	// nucleotide mutation matrix is relrate * stationary
	IIDExp* exprelrate;
	RenormalizedPosRealVector* renormrelrate;
	// Dirichlet* relrate;
	Dirichlet* tsrelrate;
	Dirichlet* tvrelrate;

	// for homogeneous model
	Dirichlet* freestationary;
	InstantStat* gcstationary;
	Var<Profile>* stationary;
	int gcstat;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	// for gc model

	int jitter;
	int contjitter;
	Normal*** leafstates;
	Gamma* leafvar;

	Gamma* ugam;
	GammaTree* ugamtree;
	GammaWhiteNoiseProcess*wngamtree;

	int whitenoise;
	WhiteNoiseProcess* wntree;
	Gamma* wnvar;
	BranchVarTree<UnitReal>* gctree;
	Beta* rootgc;
	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;

	// for gc model
	MeanLogitTreeFromMultiVariate* gctree1;
	MeanLogitTreeFromMultiVariate* gctree2;
	MeanLogitTreeFromMultiVariate* gctree3;
	Beta* rootgc1;
	Beta* rootgc2;
	Beta* rootgc3;
	GCStatTree* stattree1;
	GCStatTree* stattree2;
	GCStatTree* stattree3;
	NucMatrixTree* nucmatrixtree1;
	NucMatrixTree* nucmatrixtree2;
	NucMatrixTree* nucmatrixtree3;

	// for both
	BranchValPtrTree<RandomSubMatrix>* matrixtree;
	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	// if true: OUP else BP
	bool autoregressive;

	bool clamproot;
	bool clamptree;

	bool meanexp;

	int omegaratiotree;
	// 0 : rate model
	// 1 : dS omega
	// 2 : dS dN
	// 3 : dS omegatsgc omegats0 omegatv

	bool separatesyn;
	bool separateomega;

	// if true: gc model
	int gc;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int matrixtype;
	// 0 : nucleotide matrix or amino acid matrix
	// 1 : MG Codon matrix
	// 2 : MG3 Codon matrix (3 different nucleotide processes at the three positions
	int  conjpath;
	bool priorsampling;
	double mappingfreq;

	bool normalise;

	int ncycle;
	int nrep;

	int mutmodel;
	// 0 : gtr
	// 1 : hky with variations of Ts/Tv along the tree
	// 2 3 : variants of hky with more complex underlying (constant) relative exchangeabilities
	// 4 : variations of Ts/Tv along the tree, with an underlying gtr matrix
	// 5 : variations of Ts/Tv0 and TvGC/Tv0 along the tree, with an underlying gtr matrix

	int dsindex;
	int gcindex;
	int tstvindex;
	int tvgcindex;
	int omegaindex;
	int omegatsindex;
	int omegatv0index;
	int omegatvgcindex;

	int df;

	double maxdrift;
	double maxphi;
	int uniformprior;

	string bounds;
	bool withdrift;
	bool withexpdrift;
	bool withreldrift;

	int clampsuffstat;
	string suffstatfile;

	bool withtimeline;

	int krkctype;
	int splitkrkctype;
	GeneticCodeType codetype;
	CodonStateSpace* codonstatespace;

	SimilarityMatrix * aasimilarityMatrix ;
	SplitAAMatrix* aasplitMatrix;

    int TAG1 = 0;

	public:

	TimeLineMultiVariateTreeProcess* GetTimeLineMultiVariateTreeProcess() {
		TimeLineMultiVariateTreeProcess* tmp = dynamic_cast<TimeLineMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : in dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	SequenceAlignment* GetData()	{
		if (codondata)	{
			return codondata;
		}
		return nucdata;
	}

	int GetNtaxa()	{
		return taxonset->GetNtaxa();
	}

	int GetNstate()	{
		return GetData()->GetNstate();
	}

	SplitTree* GetSplitTree()	{
		SplitTree* tmp = dynamic_cast<SplitTree*>(splittree);
		if (! tmp)	{
			cerr << "error in GetSplitTree : null pointer\n";
			cerr << splittree << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	bool Split()	{
		return Ninterpol != 1;
	}

	bool Unconstrained()	{
		return (syngammatree != 0);
	}

	bool isCalibrated()	{
		return iscalib;
	}

	BranchOmegaMultivariateModel() {
	}

	BranchPartition* GetPartition()	{
		if (partition && Split())	{
			return splitpartition;
		}
		return partition;
	}

	bool isJittered()	{
		return jitter;
	}

	bool isContJittered()	{
		return contjitter;
	}


	BranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double insofta,  double inmeanchi, double inmeanchi2, double priorsigma, string priorsigmafile, int indf, int inmutmodel, int ingc, bool inclampdiag, bool inautoregressive, int inconjpath, double inmappingfreq, int contdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inncycle, string inbounds, string inmix, int inNinterpol, int inwithdrift, int inuniformprior, string rootfile,  string insuffstatfile, bool intimeline, bool inseparatesyn, bool inseparateomega, int inkrkctype, int injitter, int inmyid, int innprocs, int insample, GeneticCodeType type)	{

		sample = insample;

		nprocs = innprocs;
		myid = inmyid;

		dsindex = -1;
		whitenoise = false;

		contjitter = 0;
		if (injitter > 2)	{
			injitter -= 2;
			contjitter = 1;
		}
		jitter = injitter;

		// splitkrkctype == 1: 0 rates for non nearest neighbors
		// splitkrkctype == 2: only one omega
		krkctype = inkrkctype;
		splitkrkctype = 0;
		if (krkctype >= 8)	{
		}
		else if (krkctype >= 4)	{
			splitkrkctype = 1;
			krkctype -= 4;
			krkctype = inkrkctype;
		}
		codetype = type;
		codonstatespace = new CodonStateSpace(codetype);

		separatesyn = inseparatesyn;
		separateomega = inseparateomega;

		withtimeline = intimeline;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		maxphi = 10.0;
		maxdrift = 100.0;
		uniformprior = inuniformprior;

		if (inwithdrift == 4)	{
			withdrift = false;
			withexpdrift = false;
			withreldrift = true;
		}
		else if (inwithdrift == 3)	{
			withdrift = true;
			withexpdrift = true;
			withreldrift = true;
		}
		else if (inwithdrift == 2)	{
			withdrift = true;
			withexpdrift = true;
			withreldrift = false;
		}
		else if (inwithdrift == 1)	{
			withdrift = true;
			withexpdrift = false;
			withreldrift = false;
		}
		else {
			withdrift = false;
			withexpdrift = false;
			withreldrift = false;
		}

		Ninterpol = inNinterpol;
		bounds = inbounds;
		mix = inmix;

		df = indf;
		mutmodel = inmutmodel;

		normalise = innormalise;

		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;


		// A FROZEN ACCIDENT...
		if (inomegaratiotree == 3)	{
			omegaratiotree = 4;
		}
		else if (inomegaratiotree == 4)	{
			omegaratiotree = 3;
		}
		else	{
			omegaratiotree = inomegaratiotree;
		}


		if (withtimeline && ! clamptree)	{
			cerr << "error : should clamp the tree to use timeline\n";
			exit(1);
		}

		gc = ingc;
		if (gc == 2)	{
			gc = 1;
			whitenoise = true;
		}

		if (gc == 3)	{
			if (! omegaratiotree)	{
				cerr << "error : can apply gc3 formalism only under codon model\n";
				cerr << "to activate codon model, use -dSdN or -dSomega options\n";
				exit(1);
			}
		}
		gcstat = 0;
		if (gc == -1)	{
			gc = 0;
			gcstat = 1;
		}

		// 4 : 3 independent processes : ts tv0 tvgc
		// OLD CODE
		// 5 : 2 independent processes : mixture of the log
		// 6 : 2 independent processes : mixture of the exp

		// NEW CODE
		// 5 : ts and tv0 but only one omega
		if (omegaratiotree == 5) 	{
			L = 3 + gc;
		}
		else if (omegaratiotree == 4) 	{
			L = 6 + gc;
		}
		else if (omegaratiotree == 3)	{
			// transition transversion
			L = 4 + gc;
		}
		else if (omegaratiotree)	{
			L = 2 + gc;
		}
		else	{
			L = 1 + gc;
		}

		if (mutmodel == 5)	{
			L += 2;
		}
		else if (mutmodel == 4)	{
			L += 1;
		}
		else if (mutmodel == 3)	{
			L += 3;
		}
		else if (mutmodel)	{
			L++;
		}

		if (separatesyn)	{
			L--;
		}
		if (separateomega)	{
			L--;
		}

		priorsampling = false;

		if (inconjpath == -1)	{
			conjpath = 2;
		}
		else if (inconjpath == 3)	{
			conjpath = 0;
			priorsampling = true;
		}
		else	{
			conjpath = inconjpath;
		}
		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 30 : 1;
		}
		ncycle = inncycle;

		matrixtype = 0;
		mappingfreq = inmappingfreq;

		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		rawNstate = nucdata->GetNstate();

		if (rawNstate == Naa)	{
			if ((! omegaratiotree) || (krkctype == -1))	{
				cerr << "error : with amino acid data, should specify a krkc model\n";
				exit(1);
			}
		}

		if (omegaratiotree && (rawNstate == Nnuc))	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();

		if (GetNprocs() > 1)	{
			MakeMPI();
		}

		taxonset = nucdata->GetTaxonSet();

		if (GetMyid())	{
			nucdata->SubSelect(GetSiteMin(),GetSiteMax());
			if (codondata)	{
				codondata->SubSelect(GetSiteMin(),GetSiteMax());
			}
		}

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		if (Split())	{
			if (withexpdrift || withreldrift)	{
				cerr << "error : split and exp drift not yet compatible\n";
				exit(1);
			}
			cerr << "subdivide\n";
			// tree->Subdivide(tree->GetRoot(),Ninterpol);
			splittree = new SplitTree(tree,Ninterpol);
		}
		else	{
			splittree = tree;
		}

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

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		iscalib = false;

		syngammatree = 0;
		if (calibfile == "Unconstrained")	{
			if (withexpdrift || withreldrift)	{
				cerr << "error : unconstrained and exp drift not compatible\n";
				exit(1);
			}
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
					if (meanchi != -1)	{
						MeanChi = new Const<PosReal>(meanchi);
						MeanChi2 = new Const<PosReal>(meanchi2);
						Chi = new Exponential(MeanChi,Exponential::MEAN);
						Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
					}
					else	{
						double min = 1e-6;
						double max = 1e6;
						Chi = new Jeffreys(min,max,Zero);
						Chi2 = new Jeffreys(min,max,Zero);
					}
					chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,a,b,calibset,chronoprior,softa);
				}
			}
			else	{
				if (chronoprior == 0)	{
					chronogram = new Chronogram(tree,mu);
				}
				else {
					if (meanchi != -1)	{
						MeanChi = new Const<PosReal>(meanchi);
						MeanChi2 = new Const<PosReal>(meanchi2);
						Chi = new Exponential(MeanChi,Exponential::MEAN);
						Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
					}
					else	{
						double min = 1e-6;
						double max = 1e6;
						Chi = new Jeffreys(min,max,Zero);
						Chi2 = new Jeffreys(min,max,Zero);
					}
					chronogram = new BDChronogram(tree,mu,Chi,Chi2);
				}
			}

			if (clamptree)	{
				chronogram->Clamp();
			}

			if (Split())	{
				lengthtree = new SplitLengthTree(chronogram,GetSplitTree());
			}
			else	{
				lengthtree = chronogram;
			}

		}

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
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
			gammatree = 0;
		}

		Nmat = GetPartition()->GetNComponent();

		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -2)	{
			ifstream is(priorsigmafile.c_str());
			for (int i=0; i<Ncont+L; i++)	{
				double tmp;
				is >> tmp;
				DiagArray->ClampAt(tmp,i);
			}
		}
		else if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		DiagArray0 = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -2)	{
			ifstream is(priorsigmafile.c_str());
			for (int i=0; i<Ncont+L; i++)	{
				double tmp;
				is >> tmp;
				DiagArray0->ClampAt(tmp,i);
			}
		}
		else if (priorsigma == -1)	{
			DiagArray0->ClampAt(1.0);
		}
		else	{
			DiagArray0->ClampAt(priorsigma);
		}

		sigmaZero = new SigmaZero(DiagArray);
		sigmaZero0 = new SigmaZero(DiagArray0);

		if (clampdiag)	{
			sigmaarray = new DiagonalWishartArray(Nmat, df, sigmaZero, sigmaZero0);
		}
		else	{
			sigmaarray = new InverseWishartArray(Nmat, df, sigmaZero, sigmaZero0);
		}

		sigma = sigmaarray->GetVal(0);

		// driftarray = new MultiVarArray(Nmat, Zero, ContDiagArray, ContDiagArray0);
		if (uniformprior)	{
			MultiUniArray* tmparray = new MultiUniArray(Nmat, Ncont+L, Zero, maxdrift);
			driftarray = tmparray;
			if (! withdrift)	{
				tmparray->ClampAtZero();
			}
		}
		else	{
			MultiVarArray* tmparray = new MultiVarArray(Nmat, Zero, DiagArray, DiagArray0);
			driftarray = tmparray;
			if (! withdrift)	{
				tmparray->ClampAtZero();
			}
		}
		drift = driftarray->GetVal(0);

		driftphiarray = 0;
		driftphi = 0;
		driftarray2 = 0;
		drift2 = 0;
		driftphiarray2 = 0;
		driftphi2 = 0;

		if (withexpdrift)	{
			if (uniformprior)	{
				driftphiarray = new PosUniIIDArray(Nmat,One,maxphi);
			}
			else	{
				driftphiarray = new GammaIIDArray(Nmat,One,One);
			}
			driftphi = driftphiarray->GetVal(0);
		}

		if (withreldrift)	{
			if (uniformprior)	{
				driftarray2 = new MultiUniArray(Nmat, Ncont+L, Zero, maxdrift);
				driftphiarray2 = new PosUniIIDArray(Nmat,One,maxphi);
			}
			else	{
				driftarray2 = new MultiVarArray(Nmat, Zero, DiagArray, DiagArray0);
				driftphiarray2 = new GammaIIDArray(Nmat,One,One);
			}
			drift2 = driftarray2->GetVal(0);
			driftphi2 = driftphiarray2->GetVal(0);
		}


		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(Ncont+L));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont+L));
			ifstream is(rootfile.c_str());
			if (! is)	{
				cerr << "error: cannot open root file " << rootfile << '\n';
				exit(1);
			}
			cerr << "root constraints \n";
			for (int i=0; i<Ncont + L; i++)	{
				(*rootmean)[i] = 0;
				(*rootvar)[i] = 0;
			}
			int  n;
			is >> n;
			if (n != Ncont)	{
				cerr << "error : number of entries in root calibration file " << rootfile << " does not match number of quantitative traits (" << Ncont << ")\n";
				exit(1);
			}
			for (int i=0; i<Ncont; i++)	{
				double mean, stdev;
				is >> mean >> stdev;
				(*rootmean)[i+L] = mean;
				(*rootvar)[i+L] = stdev * stdev;
				cerr << GetContinuousData()->GetCharacterName(i) << '\t' << mean << '\t' << stdev << '\n';
			}
		}

		if (withtimeline)	{
			if (Unconstrained())	{
				cerr << "error : timeline and unconstrained incompatible\n";
				exit(1);
			}
			TimeLineDiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
			TimeLineDiagArray->ClampAt(1.0);
			TimeLineSigmaZero = new SigmaZero(TimeLineDiagArray);
			timelinesigma = new DiagonalCovMatrix(TimeLineSigmaZero, Ncont+L+df);
			timelinesigma->SetIdentity();

			timeintervals = new TimeIntervals(chronogram);
			timeline = new TimeLine(chronogram,timeintervals, timelinesigma);
		}

		ugam = 0;
		ugamtree = 0;
		wngamtree = 0;
		if (jitter == 1)	{
			ugam = new Gamma(One,One);
			ugam->setval(10.0);
			ugamtree = new GammaTree(tree,ugam,ugam);
		}
		if (jitter == 2)	{
			ugam = new Gamma(One,One);
			ugam->setval(10.0);
			wngamtree = new GammaWhiteNoiseProcess(lengthtree,ugam);
			wngamtree->Reset();
		}

		synsigma = 0;
		omegasigma = 0;
		lognormalsyntree = 0;
		lognormalomegatree = 0;
		if (separatesyn)	{
			synsigma = new Jeffreys(mindiag,maxdiag,Zero);
			synsigma->setval(1);
			lognormalsyntree = new LogNormalTreeProcess(lengthtree,synsigma,INTEGRAL);
		}
		if (separateomega)	{
			omegasigma = new Jeffreys(mindiag,maxdiag,Zero);
			omegasigma->setval(1);
			lognormalomegatree = new LogNormalTreeProcess(lengthtree,omegasigma,MEAN);
		}

		if (autoregressive)	{
			phi = new Gamma(One,One);
			mean = new IIDUniform(Zero,Ncont+L,100);
			if (gammamixtree || gammatree)	{
				cerr << "error: cannot use partitions for autocorrelated processes\n";
				exit(1);
			}
			process = new AutoRegressiveMultiVariateTreeProcess(sigma,mean,phi,lengthtree);
		}
		else	{
			phi = 0;
			mean = 0;
			if (gammatree)	{
				// process = new PartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammatree,driftarray,driftphiarray, chronogram);
				process = new PartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammatree,driftarray,driftphiarray,chronogram, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);
			}
			else if (withtimeline)	{
				process = new TimeLineMultiVariateTreeProcess(sigma,timeline,chronogram);
			}
			else	{
				///process = new PartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammamixtree,driftarray,driftphiarray, chronogram);
				process = new PartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammamixtree,driftarray,driftphiarray,chronogram, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);
			}
		}

		if (whitenoise)	{
			wnvar = new Gamma(One,One);
			wntree = new WhiteNoiseProcess(lengthtree,Zero,wnvar);
		}
		else	{
			wnvar = 0;
			wntree = 0;
		}

		cerr << "process RESET\n";
		process->Reset();

		if (clamproot)	{
			process->ClampRoot();
		}

		// process->ClampTransversal(1,0);

		leafstates = 0;
		leafvar = 0;
		if (contjitter && (! contdata))	{
			cerr << "error: should have continuous data with contjitter option\n";
			exit(1);
		}

		if (contdata)	{
			if (contjitter)	{
				if (contdatatype)	{
					cerr << "error: model with leaf polymorphism only with log transformation\n";
					exit(1);
				}
				leafvar = new Gamma(One,One);
				leafvar->setval(0.1);
				leafstates = new Normal**[GetNtaxa()];
				for (int i=0; i<GetNtaxa(); i++)	{
					leafstates[i] = new Normal*[Ncont];
					for (int j=0; j<Ncont; j++)	{
						leafstates[i][j] = 0;
					}
				}
				for (int i=0; i<Ncont; i++)	{
					process->SetLeafStates(leafstates,contdata,leafvar,L+i,i);
				}
			}
			else	{
				for (int i=0; i<Ncont; i++)	{
					process->SetAndClamp(contdata,L+i,i,contdatatype);
				}
			}
		}

		if (clamproot)	{
			process->ClampRoot();
		}

		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		if (bounds != "None")	{
			FileBoundSet boundset(bounds,process->GetTree());
			process->SetBounds(boundset,L);
		}

		CreateSubstitutionProcess(sample);

		if (phyloprocess)	{
			phyloprocess->Unfold();
		}
		if (sample > 0)	{
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		// register model
		RootRegister(PriorMu);
		RootRegister(Zero);
		RootRegister(One);
		if (MeanChi)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (exprelrate)	{
			RootRegister(exprelrate);
		}
		/*
		if (relrate)	{
			RootRegister(relrate);
		}
		*/
		else if (tsrelrate)	{
			RootRegister(tsrelrate);
			RootRegister(tvrelrate);
		}
		if (freestationary)	{
			RootRegister(freestationary);
		}
		if (gammatree)	{
			RootRegister(MixAlpha);
		}
		if (rootmean)	{
			RootRegister(rootmean);
			RootRegister(rootvar);
		}
		Register();

		MakeScheduler();
		PreUpdate();

	}

	void PreUpdate()	{
		if (! GetMyid())	{
			if (sample > 0)	{
				Update();
			}
			if (GetNprocs() > 1)	{
				GlobalUpdateParameters();
			}
		}
		else	{
			SlaveUpdateParameters();
			if (sample > 0)	{
				phyloprocess->Sample();
				Update();
			}
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BranchOmegaMultivariateModel() {}

	void CreateSubstitutionProcess(int sample)	{

		int index = 0;

		dsindex = -1;
		gcindex = -1;
		tstvindex = -1;
		tvgcindex = -1;
		omegaindex = -1;
		omegatsindex = -1;
		omegatv0index = -1;
		omegatvgcindex = -1;

		omegatree = 0;
		omegatstree = 0;
		omegatvgctree = 0;
		omegatv0tree = 0;
		synratetstree = 0;
		synratetvgctree = 0;
        gctree = 0;
        gctree1 = 0;
        gctree2 = 0;
        gctree3 = 0;

		MultiNormal* rootval = process->GetMultiNormal(process->GetTree()->GetRoot());

		if (! Unconstrained())	{
			if (separatesyn)	{
				synratetree = lognormalsyntree;
			}
			else	{
				synratetree = new MeanExpTreeFromMultiVariate(process,index,INTEGRAL,false,meanexp);
				dsindex = index;
				index++;
			}
			nonsynratetree = 0;
		}

		if (omegaratiotree == 0)	{
			omegatree = 0;
		}
		else if (omegaratiotree == 1)	{
			if (separateomega)	{
				omegatree = lognormalomegatree;
			}
			else	{
				if (wngamtree)	{
					omegatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp,wngamtree);
				}
				else {
					omegatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp,ugamtree);
				}
				omegaindex = index;
				index++;
			}
		}
		else if (omegaratiotree == 3)	{
			synratetstree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			rootval->ClampAt(0,index);
			tstvindex = index;
			index++;
			omegatstree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatsindex = index;
			index++;
			omegatv0tree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatv0index = index;
			omegaindex = index;
			index++;
			omegatree = omegatv0tree;
		}
		else if (omegaratiotree == 4)	{
			synratetstree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			rootval->ClampAt(0,index);
			tstvindex = index;
			index++;
			synratetvgctree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			rootval->ClampAt(0,index);
			tvgcindex = index;
			index++;

			omegatstree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatsindex = index;
			index++;
			omegatv0tree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatv0index = index;
			omegaindex = index;
			index++;
			omegatvgctree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatvgcindex = index;
			index++;
			omegatree = omegatv0tree;
		}
		else if (omegaratiotree == 5)	{
			synratetstree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			rootval->ClampAt(0,index);
			tstvindex = index;
			index++;
			omegatstree = 0;
			omegatsindex = -1;
			omegatv0tree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
			omegatv0index = index;
			omegaindex = index;
			index++;
			omegatree = omegatv0tree;
		}
		else	{
			nonsynratetree = new MeanExpTreeFromMultiVariate(process,index,INTEGRAL,false,meanexp);
			omegaindex = index;
			index++;
			omegatree = new RatioTree(nonsynratetree,synratetree);
		}

		// relrate = 0;
		tsrelrate = 0;
		tvrelrate = 0;
		exprelrate = 0;
		renormrelrate = 0;

		stationary = 0;
		freestationary = 0;
		gcstationary = 0;
		rootgc = 0;

		if (rawNstate == 20)	{
			exprelrate = new IIDExp(Naa*(Naa-1)/2);
			renormrelrate = new RenormalizedPosRealVector(exprelrate);
			if (omegaratiotree == 3)	{
				// rootval->ClampAt(0,dsindex + 1);
				rootval->ClampAt(0,omegatsindex);
				rootval->ClampAt(0,omegatv0index);
			}
			else if (omegaratiotree == 5)	{
				rootval->ClampAt(0,omegatv0index);
			}
			else	{
				rootval->ClampAt(0,omegaindex);
			}
		}
		else if ((mutmodel == 0) || (mutmodel >= 4))	{
			exprelrate = new IIDExp(Nnuc*(Nnuc-1)/2);
			renormrelrate = new RenormalizedPosRealVector(exprelrate);
			// relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		}
		else if (mutmodel == 2)	{
			tsrelrate = new Dirichlet(2);
			tvrelrate = new Dirichlet(4);
		}

		if (! omegaratiotree)	{
			if (mutmodel >= 4)	{
				if (gc)	{
					if (tstvindex != -1)	{
						cerr << "error in create subs process: tstvindex\n";
						exit(1);
					}
					tstvindex = index;
					tstvtree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					rootval->ClampAt(0,index);
					index++;
					if (mutmodel == 5)	{
						tvgctree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						rootval->ClampAt(0,index);
						tvgcindex = index;
						index++;
					}

					gcindex = index;
					if (whitenoise)	{
						gctree = new MeanJitteredLogitTree(process,index,wntree,false);
						index++;
					}
					else	{
						gctree = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
					}
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					if (mutmodel == 5)	{
						nucmatrixtree = new GCNuc3X3MatrixTree(renormrelrate,stattree,tstvtree,tvgctree,normalise);
						// nucmatrixtree = new GCNuc3X3MatrixTree(relrate,stattree,tstvtree,tvgctree,normalise);
					}
					else	{
						nucmatrixtree = new GCNuc2X2MatrixTree(renormrelrate,stattree,tstvtree,normalise);
						// nucmatrixtree = new GCNuc2X2MatrixTree(relrate,stattree,tstvtree,normalise);
					}
					nucmatrix = 0;
				}
				else	{
					nucmatrix = 0;
					gctree = 0;
					rootgc = 0;
					stattree = 0;
					if (gcstat)	{
						rootgc = new Beta(One,One);
						gcstationary = new InstantStat(rootgc);
						stationary = gcstationary;
					}
					else	{
						freestationary = new Dirichlet(Nnuc);
						stationary = freestationary;
					}
					tstvindex = index;
					tstvtree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					rootval->ClampAt(0,index);
					index++;
					if (mutmodel == 5)	{
						tvgctree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						rootval->ClampAt(0,index);
						index++;
					}

					if (mutmodel == 5)	{
						nucmatrixtree = new Nuc3X3MatrixTree(renormrelrate,stationary,tstvtree,tvgctree,normalise);
						// nucmatrixtree = new Nuc3X3MatrixTree(relrate,stationary,tstvtree,tvgctree,normalise);
					}
					else	{
						nucmatrixtree = new Nuc2X2MatrixTree(renormrelrate,stationary,tstvtree,normalise);
						// nucmatrixtree = new Nuc2X2MatrixTree(relrate,stationary,tstvtree,normalise);
					}
				}
			}
			else if (mutmodel)	{
				cerr << "kappa model deprecated\n";
				exit(1);
				if (gc)	{
					gcindex = index;
					if (whitenoise)	{
						gctree = new MeanJitteredLogitTree(process,index,wntree,false);
						index++;
					}
					else	{
						gctree = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
					}
					tstvindex = index;
					kappatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					index++;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new HKYGCNucMatrixTree(kappatree,stattree,tsrelrate,tvrelrate,normalise);
					nucmatrix = 0;
					stationary = 0;
				}
				else	{
					nucmatrix = 0;
					gctree = 0;
					rootgc = 0;
					stattree = 0;
					if (gcstat)	{
						rootgc = new Beta(One,One);
						gcstationary = new InstantStat(rootgc);
						stationary = gcstationary;
					}
					else	{
						freestationary = new Dirichlet(Nnuc);
						stationary = freestationary;
					}
					tstvindex = index;
					kappatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					index++;
					nucmatrixtree = new HKYNucMatrixTree(kappatree,stationary,tsrelrate,tvrelrate,normalise);
				}
			}
			else	{
				if (gc)	{
					nucmatrix = 0;
					stationary = 0;
					gcindex = index;
					if (whitenoise)	{
						gctree = new MeanJitteredLogitTree(process,index,wntree,false);
						index++;
					}
					else	{
						gctree = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
					}
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new GTRGCNucMatrixTree(renormrelrate,stattree,normalise);
					// nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
				}
				else	{
					gctree = 0;
					rootgc = 0;
					stattree = 0;
					nucmatrixtree = 0;
					if (gcstat)	{
						rootgc = new Beta(One,One);
						gcstationary = new InstantStat(rootgc);
						stationary = gcstationary;
					}
					else	{
						freestationary = new Dirichlet(Nnuc);
						stationary = freestationary;
					}
					nucmatrix = new GTRRandomSubMatrixWithNormRates(renormrelrate,stationary,normalise);
					// nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
				}
			}

			// make substitution mappings
			pathconjtree = 0;
			phyloprocess = 0;
			if (sample != -1)	{
			if (conjpath)	{
				if (gc || mutmodel)	{
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, nucmatrixtree, nucdata);
				}
				else	{
					pathconjtree = new OneMatrixPathConjugateTree(synratetree,nucmatrix,nucdata);
				}

				phyloprocess = 0;
				if (clampsuffstat)	{
					pathconjtree->ReadFromFile(suffstatfile);
				}
				else	{
					if ((GetNprocs() == 1) || (GetMyid() > 0))	{
						phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
					}
				}
			}
			else	{
				if (GetNprocs() > 1)	{
					cerr << "error: mpi only in path conjugate mode\n";
					exit(1);
				}
				pathconjtree = 0;
				if (priorsampling)	{
					phyloprocess = 0;
				}
				else if (gc || mutmodel)	{
					phyloprocess = new BranchMatrixPhyloProcess(synratetree, nucmatrixtree, nucdata);
				}
				else	{
					phyloprocess = new OneMatrixPhyloProcess(synratetree, nucmatrix, nucdata);
				}
			}
			}
		}
		else if (rawNstate == Naa)	{

			cerr << "amino acid\n";
			kappatree = 0;
			nucmatrix = 0;

			cerr << "Using the model based on ";
			switch(krkctype){
				case 0 : aasimilarityMatrix = new PolarityBasedMatrix(); cerr << "POLARITY\n"; break;
				case 1 : aasimilarityMatrix = new VolumeBasedMatrix(); cerr << "VOLUME\n"; break;
				case 2 : aasimilarityMatrix = new ChargeBasedMatrix(); cerr << "CHARGE\n"; break;
				case 3 : aasimilarityMatrix = new PolarityAndVolumeBasedMatrix(); cerr << "POLARITY and VOLUME combined\n"; break;
				case -1 : cerr<< " NO Model Specfified" << endl; exit(0);
			}

			aasimilarityMatrix->Affiche();
			cerr << "matrix\n";
			freestationary = new Dirichlet(Naa);
			stationary = freestationary;
			aasplitMatrix = 0;
			if ((omegaratiotree == 3) || (omegaratiotree == 5))	{
				aasplitMatrix = new SplitAAMatrix(codonstatespace);
				// aasplitMatrix->Print(cerr, aasimilarityMatrix);
				matrixtree = new SplitAminoAcidOmegaMatrixTree(aasplitMatrix, splitkrkctype, aasimilarityMatrix,renormrelrate, stationary, synratetstree, omegatstree, omegatv0tree, One, One, One);
			}
			else	{
				matrixtree = new AminoAcidOmegaMatrixTree(aasimilarityMatrix,renormrelrate, stationary, omegatree, One);
			}

			nucmatrixtree = 0;

			gctree1 = 0;
			rootgc1 = 0;
			stattree1 = 0;
			nucmatrixtree1 = 0;

			gctree2 = 0;
			rootgc2 = 0;
			stattree2 = 0;
			nucmatrixtree2 = 0;

			gctree3 = 0;
			rootgc3 = 0;
			stattree3 = 0;
			nucmatrixtree3 = 0;

			gctree = 0;
			rootgc = 0;
			stattree = 0;
			nucmatrixtree = 0;

			// make substitution mappings
			pathconjtree = 0;
			phyloprocess = 0;
			if (sample != -1)	{
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree (synratetree, matrixtree, nucdata);

				phyloprocess = 0;
				if (clampsuffstat)	{
					pathconjtree->ReadFromFile(suffstatfile);
				}
				else	{
					if ((GetNprocs() == 1) || (GetMyid() > 0))	{
						phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
					}
				}
			}
			else	{
				if (GetNprocs() > 1)	{
					cerr << "error: mpi only in path conjugate mode\n";
					exit(1);
				}
				pathconjtree = 0;
				if (priorsampling)	{
					phyloprocess = 0;
				}
				phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, nucdata);
			}
			}
		}
		else	{

			if (mutmodel >= 4)	{
				cerr << "error : tstv or tstvgc not compatible with dsom2 or dsom3\n";
				exit(1);
			}
			matrixtype = 1;
			if (mutmodel)	{
				if (gc == 3)	{
					matrixtype = 2;
					if (mutmodel == 3)	{
						gctree = 0;
						rootgc = 0;
						stattree = 0;
						nucmatrixtree = 0;
						nucmatrix = 0;
						stationary = 0;
						gcindex = index;
						gctree1 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						gctree2 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						gctree3 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						tstvindex = index;
						kappatree1 = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						index++;
						kappatree2 = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						index++;
						kappatree3 = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						index++;
						kappatree = 0;
						rootgc1 = new Beta(One,One);
						rootgc2 = new Beta(One,One);
						rootgc3 = new Beta(One,One);
						stattree1 = new GCStatTree(gctree1,rootgc1);
						stattree2 = new GCStatTree(gctree2,rootgc2);
						stattree3 = new GCStatTree(gctree3,rootgc3);
						nucmatrixtree1 = new HKYGCNucMatrixTree(kappatree1,stattree1,tsrelrate,tvrelrate,normalise);
						nucmatrixtree2 = new HKYGCNucMatrixTree(kappatree2,stattree2,tsrelrate,tvrelrate,normalise);
						nucmatrixtree3 = new HKYGCNucMatrixTree(kappatree3,stattree3,tsrelrate,tvrelrate,normalise);
						if (omegaratiotree >= 3)	{
							cerr << "multiple omegas with gc3 not implemented yet\n";
							exit(1);
						}
						else	{
							matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);
						}

					}
					else	{
						gctree = 0;
						rootgc = 0;
						stattree = 0;
						nucmatrixtree = 0;
						nucmatrix = 0;
						stationary = 0;
						gcindex = index;
						gctree1 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						gctree2 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						gctree3 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
						index++;
						tstvindex = index;
						kappatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
						index++;
						rootgc1 = new Beta(One,One);
						rootgc2 = new Beta(One,One);
						rootgc3 = new Beta(One,One);
						stattree1 = new GCStatTree(gctree1,rootgc1);
						stattree2 = new GCStatTree(gctree2,rootgc2);
						stattree3 = new GCStatTree(gctree3,rootgc3);
						nucmatrixtree1 = new HKYGCNucMatrixTree(kappatree,stattree1,tsrelrate,tvrelrate,normalise);
						nucmatrixtree2 = new HKYGCNucMatrixTree(kappatree,stattree2,tsrelrate,tvrelrate,normalise);
						nucmatrixtree3 = new HKYGCNucMatrixTree(kappatree,stattree3,tsrelrate,tvrelrate,normalise);
						if (omegaratiotree >= 3)	{
							cerr << "multiple omegas with gc3 not implemented yet\n";
							exit(1);
						}
						matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);

					}
				}
				else if (gc == 1)	{
					gctree1 = 0;
					rootgc1 = 0;
					stattree1 = 0;
					nucmatrixtree1 = 0;

					gctree2 = 0;
					rootgc2 = 0;
					stattree2 = 0;
					nucmatrixtree2 = 0;

					gctree3 = 0;
					rootgc3 = 0;
					stattree3 = 0;
					nucmatrixtree3 = 0;

					nucmatrix = 0;
					stationary = 0;
					gcindex = index;
					gctree = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
					index++;
					tstvindex = index;
					kappatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					index++;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new HKYGCNucMatrixTree(kappatree,stattree,tsrelrate,tvrelrate,normalise);
					if (omegaratiotree == 4)	{
						matrixtree = new GCOmega3X3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, synratetvgctree, omegatstree, omegatv0tree, omegatvgctree, One);
						// matrixtree = new GCOmega3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatstree, omegatv0tree, omegatvgctree, One);
					}
					else if (omegaratiotree == 3)	{
						matrixtree = new GCOmega2X2MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, omegatstree, omegatv0tree, One);
					}
					else	{
						matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
					}


				}
				else	{
					gctree1 = 0;
					rootgc1 = 0;
					stattree1 = 0;
					nucmatrixtree1 = 0;

					gctree2 = 0;
					rootgc2 = 0;
					stattree2 = 0;
					nucmatrixtree2 = 0;

					gctree3 = 0;
					rootgc3 = 0;
					stattree3 = 0;
					nucmatrixtree3 = 0;

					gctree = 0;
					rootgc = 0;
					stattree = 0;
					if (gcstat)	{
						rootgc = new Beta(One,One);
						gcstationary = new InstantStat(rootgc);
						stationary = gcstationary;
					}
					else	{
						freestationary = new Dirichlet(Nnuc);
						stationary = freestationary;
					}
					nucmatrix = 0;
					tstvindex = index;
					kappatree = new MeanExpTreeFromMultiVariate(process,index,MEAN,false,meanexp);
					index++;
					nucmatrixtree = new HKYNucMatrixTree(kappatree,stationary,tsrelrate,tvrelrate,normalise);
					if (omegaratiotree == 4)	{
						matrixtree = new GCOmega3X3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, synratetvgctree, omegatstree, omegatv0tree, omegatvgctree, One);
						// matrixtree = new GCOmega3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatstree, omegatv0tree, omegatvgctree, One);
					}
					else if (omegaratiotree == 3)	{
						matrixtree = new GCOmega2X2MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, omegatstree, omegatv0tree, One);
					}
					else	{
						matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
					}

				}
			}
			else	{
				if (gc == 3)	{
					matrixtype = 2;
					kappatree = 0;
					gctree = 0;
					rootgc = 0;
					stattree = 0;
					nucmatrixtree = 0;
					nucmatrix = 0;
					stationary = 0;

					gcindex = index;
					gctree1 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
					index++;
					gctree2 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
					index++;
					gctree3 = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
					index++;
					rootgc1 = new Beta(One,One);
					rootgc2 = new Beta(One,One);
					rootgc3 = new Beta(One,One);
					stattree1 = new GCStatTree(gctree1,rootgc1);
					stattree2 = new GCStatTree(gctree2,rootgc2);
					stattree3 = new GCStatTree(gctree3,rootgc3);
					nucmatrixtree1 = new GTRGCNucMatrixTree(renormrelrate,stattree1,normalise);
					nucmatrixtree2 = new GTRGCNucMatrixTree(renormrelrate,stattree2,normalise);
					nucmatrixtree3 = new GTRGCNucMatrixTree(renormrelrate,stattree3,normalise);
					/*
					nucmatrixtree1 = new GTRGCNucMatrixTree(relrate,stattree1,normalise);
					nucmatrixtree2 = new GTRGCNucMatrixTree(relrate,stattree2,normalise);
					nucmatrixtree3 = new GTRGCNucMatrixTree(relrate,stattree3,normalise);
					*/
					if (omegaratiotree >= 3)	{
						cerr << "multiple omegas with gc3 not implemented yet\n";
						exit(1);
					}
					matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);

				}
				else if (gc == 1)	{
					kappatree = 0;
					gctree1 = 0;
					rootgc1 = 0;
					stattree1 = 0;
					nucmatrixtree1 = 0;

					gctree2 = 0;
					rootgc2 = 0;
					stattree2 = 0;
					nucmatrixtree2 = 0;

					gctree3 = 0;
					rootgc3 = 0;
					stattree3 = 0;
					nucmatrixtree3 = 0;

					nucmatrix = 0;
					stationary = 0;

					gcindex = index;
					gctree = new MeanLogitTreeFromMultiVariate(process,index,MEAN,false);
					index++;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new GTRGCNucMatrixTree(renormrelrate,stattree,normalise);
					// nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
					if (omegaratiotree == 4)	{
						matrixtree = new GCOmega3X3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, synratetvgctree, omegatstree, omegatv0tree, omegatvgctree, One);
						// matrixtree = new GCOmega3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatstree, omegatv0tree, omegatvgctree, One);
					}
					else if (omegaratiotree == 3)	{
						matrixtree = new GCOmega2X2MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, synratetstree, omegatstree, omegatv0tree, One);
					}
					else	{
						matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
					}

				}
				else	{
					kappatree = 0;
					gctree1 = 0;
					rootgc1 = 0;
					stattree1 = 0;
					nucmatrixtree1 = 0;

					gctree2 = 0;
					rootgc2 = 0;
					stattree2 = 0;
					nucmatrixtree2 = 0;

					gctree3 = 0;
					rootgc3 = 0;
					stattree3 = 0;
					nucmatrixtree3 = 0;

					gctree = 0;
					rootgc = 0;
					stattree = 0;
					nucmatrixtree = 0;

					if (gcstat)	{
						rootgc = new Beta(One,One);
						gcstationary = new InstantStat(rootgc);
						stationary = gcstationary;
					}
					else	{
						freestationary = new Dirichlet(Nnuc);
						stationary = freestationary;
					}
					nucmatrix = new GTRRandomSubMatrixWithNormRates(renormrelrate,stationary,normalise);
					// nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
					if (omegaratiotree == 4)	{
						matrixtree = new Omega3X3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, synratetstree, synratetvgctree, omegatstree, omegatv0tree, omegatvgctree, One);
						// matrixtree = new Omega3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatstree, omegatv0tree, omegatvgctree, One);
					}
					else if (omegaratiotree == 3)	{
						matrixtree = new Omega2X2MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, synratetstree, omegatstree, omegatv0tree, One);
					}
					else	{
						matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);
					}

				}
			}

			// make substitution mappings
			pathconjtree = 0;
			phyloprocess = 0;
			if (sample != -1)	{
			if (conjpath == 2)	{
				if (matrixtype == 2)	{
					// pathconjtree = new MG3CodonBranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				}
				else if (matrixtype == 1)	{
					pathconjtree = new MGCodonBranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				}
				else	{
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				}
				phyloprocess = 0;
				if (clampsuffstat)	{
					pathconjtree->ReadFromFile(suffstatfile);
				}
				else	{
					if ((GetNprocs() == 1) || (GetMyid() > 0))	{
						phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
					}
				}
			}
			else if (conjpath == 1)	{
				pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				phyloprocess = 0;
				if (clampsuffstat)	{
					pathconjtree->ReadFromFile(suffstatfile);
				}
				else	{
					if ((GetNprocs() == 1) || (GetMyid() > 0))	{
						phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
					}
				}
			}
			else	{
				if (GetNprocs() > 1)	{
					cerr << "error: mpi only in path conjugate mode\n";
					exit(1);
				}
				pathconjtree = 0;
				if (priorsampling)	{
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
				}
			}
			}
		}

		if (mappingfreq == -1)	{
			if (omegaratiotree && (rawNstate == Nnuc))	{
				mappingfreq = 0.2;
			}
			else	{
				mappingfreq = 1;
			}
		}
		cerr << "mapping freq : " << mappingfreq << '\n';
	}

	int GetMyid()	{
		return myid;
	}

	int GetNprocs()	{
		return nprocs;
	}

	int GetSiteMin()	{
		return sitemin[myid];
	}

	int GetSiteMax()	{
		return sitemax[myid];
	}

	void MakeMPI()	{

		sitemin = new int[nprocs];
		sitemax = new int[nprocs];
		sitemin[0] = sitemax[0] = 0;
		int width = Nsite / (nprocs-1);
		for (int i=1; i<nprocs; i++)	{
			sitemin[i] = (i-1)*width;
			sitemax[i] = i * width;
			if (i == (nprocs-1))	{
				sitemax[i] = Nsite;
			}
			if (!myid)	{
				cerr << i << '\t' << sitemin[i] << '\t' << sitemax[i] << '\n';
			}
		}
	}

	void GlobalUpdateParameters()	{
#ifdef USE_MPI
		ostringstream os;
		os.precision(20);
		ToStream(os);
		string s = os.str();
		int size = s.size();
		char* c = new char[size];
		for (int i=0; i<size; i++)	{
			c[i] = s[i];
		}
		MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(c,s.size(),MPI_CHAR,0,MPI_COMM_WORLD);
		delete[] c;

		// receive total log likelihood
		lnL = 0;
		MPI_Status stat;
		for (int i=1; i<nprocs; i++)	{
			double tmp;
			MPI_Recv(&tmp,1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
			lnL += tmp;
		}
#endif
	}

	void SlaveUpdateParameters()	{
#ifdef USE_MPI
		int size;
		MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
		char* c = new char[size];
		MPI_Bcast(c,size,MPI_CHAR,0,MPI_COMM_WORLD);
		ostringstream os;
		for (int i=0; i<size; i++)	{
			os << c[i];
		}
		delete[] c;
		string s = os.str();
		if (s.size() != size)	{
			cerr << "error in slave update params\n";
			cerr << size << '\t' << s.size() << '\n';
			exit(1);
		}
		istringstream is(s);
		FromStream(is);
		Update();

		// send log likelihood
		SimpleGetLogLikelihood();
		MPI_Send(&lnL,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
#endif
	}

	void GlobalUpdateSuffStat()	{
#ifdef USE_MPI
		if (pathconjtree)	{
		pathconjtree->InactivateSufficientStatistic();
		pathconjtree->ActivateSufficientStatistic();
		pathconjtree->ResetSufficientStatistic();
		MPI_Status stat;
		for (int i=1; i<nprocs; i++)	{
			int size;
			MPI_Recv(&size,1,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
			char* c = new char[size];
			MPI_Recv(c,size,MPI_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
			ostringstream os;
			for (int i=0; i<size; i++)	{
				os << c[i];
			}
			delete[] c;
			string s = os.str();
			if (s.size() != size)	{
				cerr << "error in global update suffstat\n";
				cerr << i << '\t' << size << '\t' << s.size() << '\n';
				exit(1);
			}
			istringstream is(s);
			pathconjtree->AddSuffStatFromStream(is);
		}
		}
		Update();
		// cerr << "master : " << pathconjtree->GetLogProb() << '\n';
#endif
	}

	void SlaveUpdateSuffStat()	{
#ifdef USE_MPI
		if (pathconjtree)	{
		// cerr << "slave : " << pathconjtree->GetLogProb() << '\n';
		ostringstream os;
		os.precision(20);
		pathconjtree->PrintSuffStat(os);
		string s = os.str();
		int size = s.size();
		char* c = new char[size];
		for (int i=0; i<size; i++)	{
			c[i] = s[i];
		}
		MPI_Send(&size,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
		MPI_Send(c,size,MPI_CHAR,0,TAG1,MPI_COMM_WORLD);
		delete[] c;
		}
#endif
	}

	Tree* GetTree() {return tree;}
	Tree* GetFineGrainedTree() {return splittree;}

	void UpdateOmegaTree()	{
		if (omegaratiotree == 5)	{
			((MeanExpTreeFromMultiVariate*) omegatv0tree)->specialUpdate();
		}
		else if (omegaratiotree == 4)	{
			((MeanExpTreeFromMultiVariate*) omegatv0tree)->specialUpdate();
			((MeanExpTreeFromMultiVariate*) omegatvgctree)->specialUpdate();
			((MeanExpTreeFromMultiVariate*) omegatstree)->specialUpdate();
		}
		else if (omegaratiotree == 3)	{
			((MeanExpTreeFromMultiVariate*) omegatv0tree)->specialUpdate();
			((MeanExpTreeFromMultiVariate*) omegatstree)->specialUpdate();
		}
		else if (omegaratiotree == 2)	{
			nonsynratetree->specialUpdate();
			((RatioTree*) omegatree)->specialUpdate();
		}
		else if (omegaratiotree == 1)	{
			((MeanExpTreeFromMultiVariate*) omegatree)->specialUpdate();
		}
		else	{
			cerr << "error : update omega tree called under simple rate model\n";
			exit(1);
		}
	}

	MeanExpTreeFromMultiVariate* GetKappaTree() {return kappatree;}

	MeanExpTreeFromMultiVariate* GetSynRateTree() {
		if (Unconstrained())	{
			cerr << "error : unconstrained model\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(synratetree);
		if (! tmp)	{
			cerr << "error in get synratetree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	MeanExpTreeFromMultiVariate* GetNonSynRateOrOmegaTree() {
		if (omegaratiotree == 2)	{
			MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(nonsynratetree);
			if (! tmp)	{
				cerr << "error in get omegatree: null pointer\n";
				exit(1);
			}
			tmp->specialUpdate();
			return tmp;
		}
		if (! omegaratiotree)	{
			cerr << "error : get omega tree: model is not codon\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatree);
		if (! tmp)	{
			cerr << "error in get omegatree: null pointer\n";
			exit(1);
		}
		tmp->specialUpdate();
		return tmp;
	}

	BranchValPtrTree<Dvar<PosReal> >* GetOmegaTree() {
		return omegatree;
	}


	/*
	void return UpdateOmegaTree() {
		if (lognormalomegatree)	{
			lognormalomegatree->specialUpdate();
		}
		else	{
			MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatree);
			if (! tmp)	{
				cerr << "error in get omegatree: null pointer\n";
				exit(1);
			}
			tmp->specialUpdate();
		}
	}
	*/
	/*
	MeanExpTreeFromMultiVariate* GetOmegaTree() {
		if (! omegaratiotree)	{
			cerr << "error : get omega tree: model is not codon\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatree);
		if (! tmp)	{
			cerr << "error in get omegatree: null pointer\n";
			exit(1);
		}
		return tmp;
	}
	*/

	MeanExpTreeFromMultiVariate* GetOmegaTSTree() {
		if ((omegaratiotree < 3) || (omegaratiotree == 5))	{
			cerr << "error : get omega tree: model is not codon\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatstree);
		if (! tmp)	{
			cerr << "error in get omegatstree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	MeanExpTreeFromMultiVariate* GetOmegaTV0Tree() {
		if (omegaratiotree < 3)	{
			cerr << "error : get omegatv0tree: model is not codon\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatv0tree);
		if (! tmp)	{
			cerr << "error in get omegatree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	MeanExpTreeFromMultiVariate* GetOmegaTVGCTree() {
		if ((omegaratiotree <= 3) || (omegaratiotree == 5))	{
			cerr << "error : get omega tree: model is not codon\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(omegatvgctree);
		if (! tmp)	{
			cerr << "error in get omegatree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	const double* GetRelativeRates() {
		/*
		if (rawNstate == Naa)	{
			return renormrelrate->val().GetArray();
		}
		*/
		if (! renormrelrate)	{
			cerr << "error: relrates not defined\n";
			exit(1);
		}
		return renormrelrate->val().GetArray();
	}

	int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	MeanLogitTreeFromMultiVariate* GetMeanLogitGCTree() {return dynamic_cast<MeanLogitTreeFromMultiVariate*>(gctree);}
	BranchVarTree<UnitReal>* GetGCTree() {return gctree;}
	GCStatTree* GetGCStatTree() {return stattree;}
	MeanLogitTreeFromMultiVariate* GetGC1Tree() {return gctree1;}
	MeanLogitTreeFromMultiVariate* GetGC2Tree() {return gctree2;}
	MeanLogitTreeFromMultiVariate* GetGC3Tree() {return gctree3;}

	bool isGCActivated() {return GetGCTree();}
	bool isGC3Activated() {return GetGC3Tree();}
	int GetGCIndex() {return gcindex;}

	bool SeparateSyn()	{
		return separatesyn;
	}

	bool SeparateOmega()	{
		return separateomega;
	}

	bool HasOmega()	{
		return (omegaratiotree != 0);
	}

	bool Has1Omega()	{
		return ((omegaratiotree == 1) || (omegaratiotree == 2));
	}

	bool Has2Omega()	{
		return (omegaratiotree == 3);
	}

	bool Has3Omega()	{
		return (omegaratiotree >= 4);
	}

	int GetOmegaIndex()	{
		return omegaindex;
	}

	int GetOmegaTsIndex()	{
		return omegatsindex;
	}

	int GetOmegaTv0Index()	{
		return omegatv0index;
	}

	int GetOmegaTvGCIndex()	{
		return omegatvgcindex;
	}

	bool HasTsTv()	{
		return tstvindex != -1;
	}

	bool HasTvGC()	{
		return tvgcindex != -1;
	}

	int GetTsTvIndex()	{
		return tstvindex;
	}

	int GetTvGCIndex()	{
		return tvgcindex;
	}

	TimeLine* GetTimeLine() {return timeline;}
	TimeIntervals* GetTimeIntervals() {return timeintervals;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	Chronogram* GetChronogram() {
		if (Unconstrained())	{
			cerr << "error : chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return chronogram;
	}
	LengthTree* GetLengthTree() {return lengthtree;}
	SplitLengthTree* GetSplitLengthTree()	{
		SplitLengthTree* tmp = dynamic_cast<SplitLengthTree*>(lengthtree);
		if (! tmp)	{
			cerr << "error in GetSplitLengthTree : null pointer\n";
			cerr << lengthtree << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}


	double GetMixAlpha()	{
		if (gammatree)	{
			return MixAlpha->val();
		}
		else	{
			return 1;
		}
	}

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix(int mat = 0) {return sigmaarray->GetVal(mat);}

	CalibratedChronogram* GetCalibratedChronogram()	{
		if (Unconstrained())	{
			cerr << "error : calibrated chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double Update(bool check = false)	{
		double ret = ProbModel::Update();
		/*
		if (phyloprocess)	{
			phyloprocess->Sample();
			ret = ProbModel::Update();
		}
		*/
		return ret;
	}

	double GetLeafStatesLogProb()	{
		double total = 0;
		for (int i=0; i<GetNtaxa(); i++)	{
			for (int j=0; j<Ncont; j++)	{
				if (leafstates[i][j])	{
					total += leafstates[i][j]->GetLogProb();
				}
			}
		}
	}

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

		if (separatesyn)	{
			total += synsigma->GetLogProb();
			total += lognormalsyntree->GetLogProb();
		}
		if (separateomega)	{
			total += omegasigma->GetLogProb();
			total += lognormalomegatree->GetLogProb();
		}

		if (gammamixtree)	{
			total += gammamixtree->GetLogProb();
		}
		if (gammatree)	{
			total += MixAlpha->GetLogProb();
			total += gammatree->GetLogProb();
		}
		if (whitenoise)	{
			total += wnvar->GetLogProb();
			total += wntree->GetLogProb();
		}
		if (jitter == 1)	{
			total += ugam->GetLogProb();
			total += ugamtree->GetLogProb();
		}
		if (jitter == 2)	{
			total += ugam->GetLogProb();
			total += wngamtree->GetLogProb();
		}
		if (contjitter)	{
			total += GetLeafStatesLogProb();
			total += leafvar->GetLogProb();
		}
		if (withtimeline)	{
			total += timelinesigma->GetLogProb();
			total += timeline->GetLogProb();
		}

		// cerr << total << '\n';
		total += DiagArray->GetLogProb();
		total += DiagArray0->GetLogProb();
		total += sigmaarray->GetLogProb();
		total += driftarray->GetLogProb();
		if (withexpdrift)	{
			total += driftphiarray->GetLogProb();
		}
		if (withreldrift)	{
			total += driftphiarray2->GetLogProb();
			total += driftarray2->GetLogProb();
		}
		if (autoregressive)	{
			total += phi->GetLogProb();
			total += mean->GetLogProb();
		}
		// cerr << total << '\n';
		total += process->GetLogProb();
		// cerr << total << '\n';

		if (rawNstate == Naa)	{
			total += exprelrate->GetLogProb();
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
			total += exprelrate->GetLogProb();
			// total += relrate->GetLogProb();
		}
		else if (mutmodel == 2)	{
			total += tsrelrate->GetLogProb();
			total += tvrelrate->GetLogProb();
		}

		if (gc == 3)	{
			total += rootgc1->GetLogProb();
			total += rootgc2->GetLogProb();
			total += rootgc3->GetLogProb();
		}
		if (rootgc)	{
			total += rootgc->GetLogProb();
		}
		if (freestationary)	{
			total += freestationary->GetLogProb();
		}
		// cerr << total << '\n';
		// cerr << "---\n";
		return total;
	}

	double GetLogLikelihood()	{
		if ((GetNprocs() == 1) || (GetMyid() > 0))	{
			return SimpleGetLogLikelihood();
		}
		return lnL;
		// return GlobalGetLogLikelihood();
	}

	double SimpleGetLogLikelihood()	{
		lnL = 0;
		if (priorsampling || (sample == -1))	{
			return 0;
		}
		else if (clampsuffstat)	{
			if (pathconjtree)	{
				lnL = pathconjtree->GetLogProb();
			}
			lnL = phyloprocess->GetLogProb();
		}
		else	{
			if (pathconjtree)	{
				lnL = pathconjtree->GetLogProb();
			}
			lnL = phyloprocess->GetLogProb();
		}
		return lnL;
	}

	double Move(double tuning_modulator = 1)	{

		if (GetNprocs() == 1)	{
			// return ProbModel::Move(tuning_modulator);
			// scheduler.Cycle(1,ncycle,true,true);
			scheduler.Cycle(1,ncycle,false,false);
			return 0;
		}
		double ret = 0;
		if (GetMyid())	{
			scheduler.Cycle(1,1,false,false);
			// ret = ProbModel::Move(tuning_modulator);
			SlaveUpdateSuffStat();
			SlaveUpdateParameters();
		}
		else	{
			GlobalUpdateSuffStat();
			scheduler.Cycle(1,ncycle,false,false);
			// ret = ProbModel::Move(tuning_modulator);
			GlobalUpdateParameters();
		}
		return ret;
	}

	void MasterReceiveSuffStat()	{
		GlobalUpdateSuffStat();
	}

	void MasterMoveParams()	{
		scheduler.Cycle(1,ncycle,false,false);
	}

	void MasterSendParams()	{
		GlobalUpdateParameters();
	}

	void SlaveSendSuffStat()	{
		SlaveUpdateSuffStat();
	}

	void SlaveResampleSub()	{
		scheduler.Cycle(1,1,false,false);
	}

	void SlaveReceiveParams()	{
		SlaveUpdateParameters();
	}


	virtual void MakeScheduler()	{

		if (sample != -1)	{
			if (conjpath)	{
				if (! clampsuffstat)	{
					if ((GetNprocs() == 1) || (GetMyid() > 0))	{
						scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree,mappingfreq),1,"mapping + sufficient stat");
					}
				}
			}
			else	{
				if (phyloprocess)	{
					scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
				}
			}
		}

		if (! GetMyid())	{

		vector <SplitMultiVariateNodeMove*> nodesplitarray;
		vector <SplitMultiVariateBranchMove*> branchsplitarray;
		if (Split())	{
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,10));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,0.1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,0.01));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,10));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,0.1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,0.01));
		}

		for (int i=0; i<nrep; i++)	{
			if (Unconstrained())	{
				scheduler.Register(new SimpleMove(mu,1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(mu,0.1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(syngammatree,3),30,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.3),30,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.03),30,"syngamtree");
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

				// CHRONOCOMP
				/*
				if (separatesyn)	{
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,1),10,"chrono lognormalsynrate comp");
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,0.1),10,"chrono lognormalsynrate comp");
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,0.01),10,"chrono lognormalsynrate comp");
				}
				else	{
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,1),10,"chrono process comp");
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,0.1),10,"chrono process comp");
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,0.01),10,"chrono process comp");
				}
				*/
			}

			if (gammamixtree)	{
				scheduler.Register(new SimpleMove(gammamixtree,10),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,1),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,0.1),10,"gammamixtree");
			}
			if (gammatree)	{
				scheduler.Register(new SimpleMove(MixAlpha,10),10,"mixalpha");
				scheduler.Register(new SimpleMove(MixAlpha,1),10,"mixalpha");
				scheduler.Register(new SimpleMove(MixAlpha,0.1),10,"mixalpha");
				scheduler.Register(new SimpleMove(gammatree,10),10,"gammatree");
				scheduler.Register(new SimpleMove(gammatree,1),10,"gammatree");
				scheduler.Register(new SimpleMove(gammatree,0.1),10,"gammatree");

			}
			if (whitenoise)	{
				scheduler.Register(new SimpleMove(wntree,10),10,"gammatree");
				scheduler.Register(new SimpleMove(wntree,1),10,"gammatree");
				scheduler.Register(new SimpleMove(wntree,0.1),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,10),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,1),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,0.1),10,"gammatree");
			}
			if (isCalibrated())	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");

				/*
				scheduler.Register(new RootMove(GetCalibratedChronogram(),1),10,"root age additive");
				scheduler.Register(new RootMove(GetCalibratedChronogram(),0.1),10,"root age additive");
				scheduler.Register(new RootMove(GetCalibratedChronogram(),0.01),10,"root age additive");
				*/

				/*
				if (separatesyn)	{
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,1),10,"chrono scale lognormalsynrate comp");
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,0.1),10,"chrono scale lognormalsynrate comp");
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,0.01),10,"chrono scale lognormalsynrate comp");
				}
				else	{
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,1),10,"chrono scale process comp");
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");
				}
				*/
			}

			if (separatesyn)	{
				scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");
				scheduler.Register(new SimpleMove(lognormalsyntree,10),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,1),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.1),10,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.01),10,"log normal syn");
			}

			if (separateomega)	{
				scheduler.Register(new SimpleMove(omegasigma,10),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,1),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,0.1),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,0.01),100,"omega sigma");
				scheduler.Register(new SimpleMove(lognormalomegatree,10),10,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,1),10,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,0.1),10,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,0.01),10,"log normal omega");
			}

			if (jitter == 1)	{
				scheduler.Register(new SimpleMove(ugamtree,10),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugamtree,1),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugamtree,0.1),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugam,10),10,"ugam");
				scheduler.Register(new SimpleMove(ugam,1),10,"ugam");
				scheduler.Register(new SimpleMove(ugam,0.1),10,"ugam");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,1),10,"ugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,0.1),10,"ugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,0.01),10,"ugam process comp");
			}

			if (jitter == 2)	{
				scheduler.Register(new SimpleMove(wngamtree,10),10,"wntree");
				scheduler.Register(new SimpleMove(wngamtree,1),10,"wntree");
				scheduler.Register(new SimpleMove(wngamtree,0.1),10,"wntree");
				scheduler.Register(new SimpleMove(ugam,10),10,"wn");
				scheduler.Register(new SimpleMove(ugam,1),10,"wn");
				scheduler.Register(new SimpleMove(ugam,0.1),10,"wn");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,1),10,"wnugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,0.1),10,"wnugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,0.01),10,"wnugam process comp");
			}

			if (contjitter)	{
				scheduler.Register(new SimpleMove(leafvar,1),10,"leafvar");
				scheduler.Register(new SimpleMove(leafvar,0.1),10,"leafvar");
				scheduler.Register(new SimpleMove(leafvar,0.01),10,"leafvar");
			}

			scheduler.Register(new SimpleMove(sigmaarray,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigmaarray,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigmaarray,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigmaarray,0.01),100,"sigma");

			/*
			int n = taxonset->GetNtaxa() * 10;
			scheduler.Register(new MultiVariatePropagateMove(process,1,0.1,0.1),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.5,0.5),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.9,0.9),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.01,0.99,0.99),n,"propmove");
			*/

			if (withtimeline)	{
				scheduler.Register(new SimpleMove(timelinesigma,10),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,1),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,0.1),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,0.01),100,"timelinesigma");

				scheduler.Register(new SimpleMove(timeline,10),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,1),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,0.1),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,0.01),100,"timeline");
			}

			if (Split())	{

				scheduler.Register(nodesplitarray[0],10,"split node process");
				scheduler.Register(nodesplitarray[1],10,"split node process");
				scheduler.Register(nodesplitarray[2],10,"split node process");
				scheduler.Register(nodesplitarray[3],10,"split node process");
				scheduler.Register(branchsplitarray[0],10,"split node process");
				scheduler.Register(branchsplitarray[1],10,"split node process");
				scheduler.Register(branchsplitarray[2],10,"split node process");
				scheduler.Register(branchsplitarray[3],10,"split node process");

			}

			if (withdrift)	{
				scheduler.Register(new SimpleMove(driftarray,10),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,0.1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,0.01),10,"drift");
			}

			if (withexpdrift)	{
				scheduler.Register(new SimpleMove(driftphiarray,10),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,0.1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,0.01),10,"drift exp relax");
			}
			if (withreldrift)	{
				scheduler.Register(new SimpleMove(driftphiarray2,10),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,0.1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,0.01),10,"drift exp relax");

				scheduler.Register(new SimpleMove(driftarray2,10),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,0.1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,0.01),10,"drift");
			}

			// scheduler.Register(new SimpleMove(process,10),10,"multinormal");
			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

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

			if (rawNstate == Naa)	{
				scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.1,2),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.03,4),10,"relrates");
				scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");

				if ((omegaratiotree == 3) || (omegaratiotree == 5))	{
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype,aasimilarityMatrix,1,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype, aasimilarityMatrix,0.1,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype, aasimilarityMatrix,0.01,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
				}
				else	{
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,1,omegaindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.1,omegaindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.01,omegaindex),10,"process relrate comp move");
				}
			}
			else	{

				if ((! mutmodel) || (mutmodel >= 4))	{
				// if (! mutmodel)	{
					/*
					scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
					scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
					*/
					scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.1,1),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.03,2),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.01,4),10,"relrates");
					scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");
				}
				else if (mutmodel == 2)	{
					scheduler.Register(new ProfileMove(tsrelrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(tvrelrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(tvrelrate,0.03,2),10,"relrates");
				}

				if ((mutmodel >= 4) || (omegaratiotree >= 4))	{
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,1,tstvindex,tvgcindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,0.1,tstvindex,tvgcindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,0.1,tstvindex,tvgcindex),10,"process relrate comp move");
				}

			}

			if (gc == 3)	{
				scheduler.Register(new SimpleMove(rootgc1,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc1,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc1,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,0.01),10,"root gc");
			}
			if (rootgc)	{
				scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");
			}
			if (freestationary)	{
				scheduler.Register(new ProfileMove(freestationary,0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(freestationary,0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(freestationary,0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(freestationary,0.001),10,"stat");
			}

			/*
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,1,0),5,"root compensation");
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.1,0),5,"root compensation");
			scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.01,0),5,"root compensation");
			*/

		}
		}
	}

	void CheckCal()	{
		if (iscalib && GetCalibratedChronogram()->CheckBounds())	{
			cerr << "before move:\n";
			cerr << "error : some nodes are outside of calibration intervals\n";
			cerr << GetCalibratedChronogram()->CheckBounds() << '\n';
			exit(1);
		}
	}

	void drawSample()	{
		cerr << "check draw sample\n";
		exit(1);
		cerr << "sample\n";

		if (Unconstrained())	{
			mu->Sample();
			syngammatree->Sample();
		}
		else	{
			if (chronoprior)	{
				Chi->Sample();
				Chi2->Sample();
			}
			chronogram->Sample();
		}

		if (gammamixtree)	{
			gammamixtree->Sample();
		}
		if (gammatree)	{
			MixAlpha->Sample();
			MixAlpha->setval(1.0);
			gammatree->Sample();
		}
		if (whitenoise)	{
			wnvar->Sample();
			wntree->Sample();
		}
		if (jitter == 1)	{
			ugam->Sample();
			ugamtree->Sample();
		}
		if (jitter == 2)	{
			ugam->Sample();
			wngamtree->Sample();
		}
		if (contjitter)	{
			leafvar->Sample();
		}

		if (withtimeline)	{
			timelinesigma->Sample();
			timeline->Sample();
		}

		if (separatesyn)	{
			synsigma->setval(1.0);
			lognormalsyntree->Sample();
		}
		if (separateomega)	{
			omegasigma->setval(1.0);
			lognormalomegatree->Sample();
		}

		DiagArray->Sample();

		sigmaarray->Sample();
		// sigma->SetIdentity();
		driftarray->Sample();
		if (withexpdrift)	{
			driftphiarray->Sample();
		}
		if (withreldrift)	{
			driftphiarray2->Sample();
			driftarray2->Sample();
		}

		if (autoregressive)	{
			phi->Sample();
			mean->Sample();
		}

		process->Sample();

		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		if (! Unconstrained())	{
			GetSynRateTree()->specialUpdate();
		}
		if (omegaratiotree)	{
			UpdateOmegaTree();
		}

		if (rawNstate == Naa)	{
			exprelrate->Sample();
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			exprelrate->Sample();
		}
		else if (mutmodel == 2)	{
			tsrelrate->Sample();
			tvrelrate->Sample();
		}

		if (gc == 3)	{
			rootgc1->Sample();
			rootgc2->Sample();
			rootgc3->Sample();
		}
		if (rootgc)	{
			rootgc->Sample();
		}
		if (freestationary)	{
			freestationary->Sample();
		}

		if (phyloprocess)	{
			phyloprocess->Sample();
		}
	}

	double GetMeanTimeLine(int k)	{
		double mean = 0;
		for (int i=0; i<timeline->GetSize(); i++)	{
			mean += (*(timeline->GetVal(i)))[k];
		}
		mean /= timeline->GetSize();
		return mean;
	}

	double GetTimeLine(int i, int k)	{
		return (*(timeline->GetVal(i)))[k];
	}

	double GetVarTimeLine(int k)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<timeline->GetSize(); i++)	{
			double tmp = (*(timeline->GetVal(i)))[k];
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= timeline->GetSize();
		var /= timeline->GetSize();
		var -= mean * mean;
		return var;
	}

	Var<RealVector>* GetDrift(int i = 0)	{
		return driftarray->GetVal(i);
	}

	Var<PosReal>* GetDriftPhi(int i = 0)	{
		if (!driftphiarray)	{
			cerr << "error in GetDriftPhi: null\n";
			exit(1);
		}
		return driftphiarray->GetVal(i);
	}

	Var<RealVector>* GetDrift2(int i = 0)	{
		return driftarray2->GetVal(i);
	}

	Var<PosReal>* GetDriftPhi2(int i = 0)	{
		if (!driftphiarray2)	{
			cerr << "error in GetDriftPhi: null\n";
			exit(1);
		}
		return driftphiarray2->GetVal(i);
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

	double GetRadConsNormFactor()	{
		double cons = 0;
		double rad = 0;
		int ncons = 0;
		int nrad = 0;
		for (int i=0; i<Naa; i++)	{
			for (int j=0; j<Naa; j++)	{
				if (i != j)	{
					if(aasimilarityMatrix->isRadical(i,j))	{
						// rad += (*exprelrate)[rrindex(i,j,Naa)];
						rad += (*exprelrate)[rrindex(i,j,Naa)] * (*stationary)[j];
						// rad += (*exprelrate)[rrindex(i,j,Naa)] * (*stationary)[i] * (*stationary)[j];
						nrad++;
					}
					else	{
						// cons += (*exprelrate)[rrindex(i,j,Naa)];
						cons += (*exprelrate)[rrindex(i,j,Naa)] * (*stationary)[j];
						ncons++;
					}
				}
			}
		}
		rad /= nrad;
		cons /= ncons;
		return rad / cons;
	}

	double GetMaxdS()	{
		if (! omegaratiotree)	{
			cerr << "error : get max dS called under simple rate model\n";
			exit(1);
		}
		return GetSynRateTree()->GetMax();
	}

	double GetMaxdN()	{
		if (! omegaratiotree)	{
			cerr << "error : get max dN called under simple rate model\n";
			exit(1);
		}
		if (! nonsynratetree)	{
			cerr << "error : cannot call getmaxdN\n";
			exit(1);
		}
		return nonsynratetree->GetMax();
	}

	double GetMaxOmega()	{
		if (! omegaratiotree)	{
			cerr << "error : get max omega called under simple rate model\n";
			exit(1);
		}
		return ((MeanExpTreeFromMultiVariate*) omegatree)->GetMax();
	}

	double GetMeanSynRate()	{
		if (Unconstrained())	{
			return syngammatree->GetMean();
		}
		if (separatesyn)	{
			return lognormalsyntree->GetTotalLength();
		}
		return GetSynRateTree()->GetTotal();
	}

	double GetMeanSynRateTS()	{
		if (mutmodel >=4)	{
			return tstvtree->GetMean();
		}
		return synratetstree->GetMean();
	}

	double GetMeanSynRateTVGC()	{
		if (mutmodel >=4)	{
			return tvgctree->GetMean();
		}
		return synratetvgctree->GetMean();
	}

	double GetTotalTime()	{
		if (Unconstrained())	{
			cerr << "error : get total time under unconstrained\n";
			exit(1);
		}
		return chronogram->GetTotalTime();
	}

	double GetMeanKappa()	{
		if (mutmodel == 3)	{
			double total = 0;
			total += kappatree1->GetMean();
			total += kappatree2->GetMean();
			total += kappatree3->GetMean();
			return total / 3;
		}
		return kappatree->GetMean();
	}

	double GetMeanOmega()	{
		if (! omegaratiotree)	{
			cerr << "error : get mean omega called under simple rate model\n";
			exit(1);
		}
		if (omegaratiotree == 2)	{
			return ((RatioTree*) omegatree)->GetMean();
		}
		else	{
			return ((MeanExpTreeFromMultiVariate*) omegatree)->GetMean();
		}
	}

	double GetMeanOmegaTVGC()	{
		if ((omegaratiotree <= 3) || (omegaratiotree == 5))	{
			cerr << "error : get mean omega called under wrong model\n";
			exit(1);
		}
		return ((MeanExpTreeFromMultiVariate*) omegatvgctree)->GetMean();
	}

	double GetMeanOmegaTS()	{
		if ((omegaratiotree < 3) || (omegaratiotree == 5))	{
			cerr << "error : get mean omega called under wrong model\n";
			exit(1);
		}
		return ((MeanExpTreeFromMultiVariate*) omegatstree)->GetMean();
	}

	double GetMeanGC1()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree1->GetMean();
	}

	double GetRootGC1()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return rootgc1->val();
	}

	double GetMeanGC2()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree2->GetMean();
	}

	double GetRootGC2()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return rootgc2->val();
	}

	double GetMeanGC3()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree3->GetMean();
	}

	double GetRootGC3()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return rootgc3->val();
	}

	double GetMeanGC()	{
		if (gc != 1)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return 0;
		// return gctree->GetMean();
	}

	double GetVarGC()	{
		if (gc != 1)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return 0;
		// return gctree->GetVar();
	}

	double GetRootGC()	{
		if (gc != 1)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return rootgc->val();
	}

	Var<RealVector>* GetMeanVector()	{
		return mean;
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tlength";
		if (omegaratiotree == 5)	{
			os << "\tTs/Tv\tdN/dS_Tv0";
		}
		else if (omegaratiotree == 4)	{
			os << "\tTs/Tv0\tTvGC/Tv0\tdN/dS_Ts\tdN/dS_Tv0\tdN/dS_TvGC";
		}
		else if (omegaratiotree == 3)	{
			os << "\tTs/Tv\tdN/dS_Ts\tdN/dS_Tv";
		}
		else if (omegaratiotree == 2)	{
			os << "\tdN";
		}
		else if (omegaratiotree == 1)	{
			os << "\tdN/dS";
		}
		else if (mutmodel == 4)	{
			os << "\tTs/Tv";
		}
		else if (mutmodel == 5)	{
			os << "\tTs/Tv0\tTvGC/Tv0";
		}
		if (gc == 3)	{
			os << "\tmeangc1\tmeangc2\tmeangc3";
		}
		else if (gc)	{
			os << "\tmeangc";
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		if (isJittered())	{
			os << '\t' << "ugam";
		}
		if (contjitter)	{
			os << '\t' << "leafvar";
		}

		for (int mat=0; mat<Nmat; mat++)	{
			for (int k=0; k<Ncont+L; k++)	{
				for (int l=k+1; l<Ncont+L; l++)	{
					os << '\t' << "sigma_" << k+1 << '_' << l+1;
				}
			}
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << "sigma_" << k+1 << '_' << k+1;
			}
		}
		if (withtimeline)	{
			for (int k=0; k<Ncont+L; k++)	{
				os << "\tmeantimeline" << k+1;
			}
			for (int k=0; k<Ncont+L; k++)	{
				os << "\tvartimeline" << k+1;
			}
		}
		if (whitenoise)	{
			os << '\t' << "wnvar";
		}
		for (int mat=0; mat<Nmat; mat++)	{
			if (withdrift)	{
				os << "\tdim";
				for (int k=0; k<Ncont+L; k++)	{
					os << '\t' << "drift_" << k+1;
				}
			}
			if (withexpdrift)	{
				os << '\t' << "phi";
			}
			if (withreldrift)	{
				for (int k=0; k<Ncont+L; k++)	{
					os << '\t' << "drift2_" << k+1;
				}
				os << '\t' << "phi2";
			}
		}

		if (autoregressive)	{
			os << '\t' << "theta";
			os << '\t' << "dim";
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << "mean_" << k+1;
			}
		}
		else	{
			os << "\tdim";
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << "root_" << k+1;
			}
		}

		if (gc == 3)	{
			os << "\trootgc1\trootgc2\trootgc3";
		}
		else if (gc)	{
			os << "\trootgc";
		}
		else	{
			os << "\tstatent";
		}
		if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			os << "\trrent";
		}
		else if (mutmodel == 2)	{
			os << "\ttsrrent\ttvrrent";
		}
		else	{
			os << "\tmeankappa";
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Ncont+L; k++)	{
				os << "\tdiag" << k;
			}
		}

		if (gammamixtree)	{
			for (int k=0; k<gammamixtree->GetNComponent(); k++)	{
				os << "\tgammamix" << k ;
			}
		}
		if (gammatree)	{
				os << "\tmixalpha";
				os << "\tgammean\tgamvar";
		}

		if (chronoprior)	{
			os << "\tchronologprior";
			os << "\tdelta\tkappa";
			// os << "\tnumerror";
		}

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		CheckCal();
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		/*
		os << '\t' << pathconjtree->GetLogProb();
		os << '\t' << phyloprocess->GetLogProb();
		*/
		os << '\t' << GetMeanSynRate();
		if (omegaratiotree == 5)	{
			os << '\t' << GetMeanSynRateTS();
			os << '\t' << GetMeanOmega();
		}
		else if (omegaratiotree == 4)	{
			os << '\t' << GetMeanSynRateTS();
			os << '\t' << GetMeanSynRateTVGC();
			os << '\t' << GetMeanOmegaTS();
			os << '\t' << GetMeanOmega();
			os << '\t' << GetMeanOmegaTVGC();
		}
		else if (omegaratiotree == 3)	{
			os << '\t' << GetMeanSynRateTS();
			os << '\t' << GetMeanOmegaTS();
			os << '\t' << GetMeanOmega();
		}
		else if (omegaratiotree)	{
			os << '\t' << GetMeanOmega();
		}
		else if (mutmodel == 4)	{
			os << '\t' << GetMeanSynRateTS();
		}
		else if (mutmodel == 5)	{
			os << '\t' << GetMeanSynRateTS();
			os << '\t' << GetMeanSynRateTVGC();
		}
		if (gc == 3)	{
			os << '\t' << GetMeanGC1();
			os << '\t' << GetMeanGC2();
			os << '\t' << GetMeanGC3();
		}
		else if (gc)	{
			os << '\t' << GetMeanGC();
		}
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
			/*
			if (chronoprior)	{
				os << '\t' << Chi->val();
				os << '\t' << Chi2->val();
			}
			*/
		}
		if (isJittered())	{
			os << '\t' << *ugam;
		}
		if (contjitter)	{
			os << '\t' << *leafvar;
		}

		for (int mat=0; mat<Nmat; mat++)	{
			for (int k=0; k<Ncont+L; k++)	{
				for (int l=k+1; l<Ncont+L; l++)	{
					os << '\t' << (*sigmaarray->GetVal(mat))[k][l];
				}
			}
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << (*sigmaarray->GetVal(mat))[k][k];
			}
		}
		if (withtimeline)	{
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << GetMeanTimeLine(k);
			}
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << GetVarTimeLine(k);
			}
		}
		if (whitenoise)	{
			os << '\t' << *wnvar;
		}
		for (int mat=0; mat<Nmat; mat++)	{
			if (withdrift)	{
				os << '\t' << *driftarray->GetVal(mat);
			}
			if (withexpdrift)	{
				os << '\t' << *driftphiarray->GetVal(mat);
			}
			if (withreldrift)	{
				os << '\t' << *driftarray2->GetVal(mat);
				os << '\t' << *driftphiarray2->GetVal(mat);
			}
		}

		os << '\t' << GetMultiVariateProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();

		if (autoregressive)	{
			os << '\t' << phi->val();
			os << '\t' << *mean;
		}

		if (gc == 3)	{
			os << '\t' << GetRootGC1();
			os << '\t' << GetRootGC2();
			os << '\t' << GetRootGC3();
		}
		else if (gc)	{
			os << '\t' << GetRootGC();
		}
		else	{
			os << '\t' << stationary->val().GetEntropy();
		}

		if (rawNstate == Naa)	{
			os << '\t' << renormrelrate->val().GetEntropy();
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			os << '\t' << renormrelrate->val().GetEntropy();
		}
		else if (mutmodel == 2)	{
			os << '\t' << tsrelrate->val().GetEntropy();
			os << '\t' << tvrelrate->val().GetEntropy();
		}
		else	{
			os << '\t' << GetMeanKappa();
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Ncont+L; k++)	{
				os << '\t' << DiagArray->GetVal(k)->val();
			}
		}
		if (gammamixtree)	{
			for (int k=0; k<gammamixtree->GetNComponent(); k++)	{
				os << '\t' << gammamixtree->GetComponent(k)->val();
			}
		}
		if (gammatree)	{
			os << '\t' << MixAlpha->val();
			os << '\t' << gammatree->GetMean() << '\t' << gammatree->GetVar();
		}

		if (chronoprior)	{
			os << '\t' << chronogram->GetLogProb();
			os << '\t' << Chi->val() << '\t' << Chi2->val();
			// os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\n';
		os.flush();
	}

	void PrintEntries(ostream& os, int* array = 0)	{

		int cumul = 0;
		if ((! separatesyn) && (! Unconstrained()))	{
		if ((! array) || (array[cumul] == 1))	{
			os << "dS\n";
		}
		}
		cumul++;
		if (! separateomega)	{
		if (omegaratiotree == 2)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "dN\n";
			}
			cumul++;
		}
		else if (omegaratiotree == 1) {
			if ((! array) || (array[cumul] == 1))	{
				os << "omega\n";
			}
			cumul++;
		}
		else if (omegaratiotree == 3)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega Ts\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega Tv\n";
			}
			cumul++;
		}
		else if (omegaratiotree == 4)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv0\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "TvGC/Tv0\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega Ts\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega Tv0\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega TvGC\n";
			}
			cumul++;
		}
		else if (omegaratiotree == 5)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv0\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "omega Tv0\n";
			}
			cumul++;
		}
		}
		if (mutmodel == 5)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv0\n";
				os << "TvGC/Tv0\n";
			}
			cumul += 2;
		}
		else if (mutmodel == 4)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv\n";
			}
			cumul += 2;
		}
		else if (mutmodel == 3)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv1\n";
				os << "Ts/Tv2\n";
				os << "Ts/Tv3\n";
			}
			cumul += 3;
		}
		else if (mutmodel)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "Ts/Tv\n";
			}
			cumul++;
		}
		if (gc == 3)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "gc1\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "gc2\n";
			}
			cumul++;
			if ((! array) || (array[cumul] == 1))	{
				os << "gc3\n";
			}
			cumul++;
		}
		else if (gc == 1)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "gc\n";
			}
			cumul++;
		}
		for (int k=0; k<Ncont; k++)	{
			if ((! array) || (array[cumul] == 1))	{
				os << GetContinuousData()->GetCharacterName(k) << '\n';
				// os << "character " << k + 1 << '\n';
			}
			cumul++;
		}
	}

	void SampleProcess()	{
		InverseWishartMatrix* mat = dynamic_cast<InverseWishartMatrix*>(sigma);
		if (! mat)	{
			cerr << "error in sample process: null matrix\n";
			exit(1);
		}
		mat->localCorrupt(false);
		mat->localUpdate();
		process->SampleFixRoot();
	}

	int GetNcheck()	{
		int tmp =  4 + (Ncont + L) * (Ncont + L+1) / 2;
		if (autoregressive)	{
			tmp += 3;
		}
		return tmp;
	}

	void GetCheck(double* p)	{
		cerr << "in get check\n";
		exit(1);
		int n = 0;
		GetSynRateTree()->specialUpdate();
		p[n] = GetMeanSynRate();
		n++;
		if (omegaratiotree)	{
			UpdateOmegaTree();
			p[n] = GetMeanOmega();
			n++;
		}
		if (autoregressive)	{
			p[n] = phi->val();
			n++;
			p[n] = mean->GetMean();
			n++;
			p[n] = mean->GetVar();
			n++;
		}

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				p[n] = (*sigma)[k][l];
				n++;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			p[n] =  (*sigma)[k][k];
			n++;
		}
		if (gc == 3)	{
			p[n] =  GetMeanGC1();
			n++;
			p[n] =  GetMeanGC2();
			n++;
			p[n] =  GetMeanGC3();
			n++;
		}
		else if (gc)	{
			p[n] =  GetMeanGC();
			n++;
		}
		else	{
			p[n] =  stationary->val().GetEntropy();
			n++;
		}
		if (rawNstate == Naa)	{
			p[n] = renormrelrate->val().GetEntropy();
			n++;
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		//if (! mutmodel)	{
			p[n] = renormrelrate->val().GetEntropy();
			n++;
		}
		if (n!= GetNcheck())	{
			cerr << "error : non matching Ncheck\n";
			cerr << n << '\t' << GetNcheck() << '\n';
			cerr << Ncont << '\n';
			exit(1);
		}
	}

	void ToStream(ostream& os)	{
		os.precision(15);
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
		if (gammamixtree)	{
			os << *gammamixtree << '\n';
		}
		if (gammatree)	{
			os << *MixAlpha << '\n';
			os << *gammatree << '\n';
		}
		if (withtimeline)	{
			os << *timelinesigma << '\n';
			os << *timeline << '\n';
		}
		if (whitenoise)	{
			os << *wnvar << '\n';
			os << *wntree << '\n';
		}
		if (jitter == 1)	{
			os << *ugam << '\n';
			os << *ugamtree << '\n';
		}
		if (jitter == 2)	{
			os << *ugam << '\n';
			os << *wngamtree << '\n';
		}
		if (contjitter)	{
			os << *leafvar << '\n';
		}

		if (separatesyn)	{
			os << *synsigma << '\n';
			os << *lognormalsyntree << '\n';
		}
		if (separateomega)	{
			os << *omegasigma << '\n';
			os << *lognormalomegatree << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigmaarray << '\n';
		os << *driftarray << '\n';
		if (withexpdrift)	{
			os << *driftphiarray << '\n';
		}
		if (withreldrift)	{
			os << *driftphiarray2 << '\n';
			os << *driftarray2 << '\n';
		}
		os << '\n';
		if (autoregressive)	{
			os << *phi << '\n';
			os << *mean << '\n';
		}
		os << '\n';
		os << *process << '\n';
		if (rawNstate == Naa)	{
			os << *exprelrate << '\n';
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			os << *exprelrate << '\n';
			// os << *relrate << '\n';
		}
		else if (mutmodel == 2)	{
			os << *tsrelrate << '\n';
			os << *tvrelrate << '\n';
		}
		if (gc == 3)	{
			os << *rootgc1 << '\n';
			os << *rootgc2 << '\n';
			os << *rootgc3 << '\n';
		}
		if (rootgc)	{
			os << *rootgc << '\n';
		}
		if (freestationary)	{
			os << *freestationary << '\n';
		}
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
		if (gammamixtree)	{
			is >> *gammamixtree;
		}
		if (gammatree)	{
			is >> *MixAlpha;
			is >> *gammatree;
		}
		if (withtimeline)	{
			is >> *timelinesigma;
			is >> *timeline;
		}
		if (whitenoise)	{
			is >> *wnvar;
			is >> *wntree;
		}
		if (jitter == 1)	{
			is >> *ugam;
			is >> *ugamtree;
		}
		if (jitter == 2)	{
			is >> *ugam;
			is >> *wngamtree;
		}
		if (contjitter)	{
			is >> *leafvar;
		}

		if (separatesyn)	{
			is >> *synsigma;
			is >> *lognormalsyntree;
		}
		if (separateomega)	{
			is >> *omegasigma;
			is >> *lognormalomegatree;
		}

		is >> *DiagArray;
		is >> *sigmaarray;
		is >> *driftarray;

		if (withexpdrift)	{
			is >> *driftphiarray;
		}
		if (withreldrift)	{
			is >> *driftphiarray2;
			is >> *driftarray2;
		}
		if (autoregressive)	{
			is >> *phi;
			is >> *mean;
		}
		is >> *process;
		if (rawNstate == Naa)	{
			is >> *exprelrate;
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			is >> *exprelrate;
			// is >> *relrate;
		}
		else if (mutmodel == 2)	{
			is >> *tsrelrate;
			is >> *tvrelrate;
		}
		if (gc == 3)	{
			is >> *rootgc1;
			is >> *rootgc2;
			is >> *rootgc3;
		}
		if (rootgc)	{
			is >> *rootgc;
		}
		if (freestationary)	{
			is >> *freestationary;
		}
	}

	/*
	int CountParam()	{
		int n = 0;
		// mu
		n++;
		if (Unconstrained())	{
			n += syngammatree->GetNBranchVals();
		}
		else	{
			n += chronogram->GetNnodeVals();
			if (isCalibrated())	{
				// is >> *GetCalibratedChronogram()->GetScale();
				n++;
			}
			if (chronoprior)	{
				// is >> *Chi >> *Chi2;
				n += 2;
			}
		}
		if (gammamixtree)	{
			n += gammamixtree->GetNBranchVals();
		}
		if (gammatree)	{
			// is >> *MixAlpha;
			n++;
			n += gammatree->GetNbranchVals();
		}
		if (withtimeline)	{
			n += timelinesigma->GetNVals();
			n += timeline->GetNvals();
		}
		if (whitenoise)	{
			// is >> *wnvar;
			n++;
			n += wntree->GetNBranchVals();
		}
		if (jitter == 1)	{
			// is >> *ugam;
			n++;
			n += ugamtree->GetNBranchVals();;
		}
		if (jitter == 2)	{
			// is >> *ugam;
			n++;
			n += wngamtree->GetNBranchVals();
		}
		if (contjitter)	{
			// is >> *leafvar;
			n++;
		}

		if (separatesyn)	{
			// is >> *synsigma;
			n++;
			n += lognormalsyntree->GetNNodeVals();
		}
		if (separateomega)	{
			// is >> *omegasigma;
			n++;
			n += lognormalomegatree->GetNNodeVals();
		}

		is >> *DiagArray;
		is >> *sigmaarray;
		is >> *driftarray;
		if (withexpdrift)	{
			is >> *driftphiarray;
		}
		if (withreldrift)	{
			is >> *driftphiarray2;
			is >> *driftarray2;
		}
		if (autoregressive)	{
			is >> *phi;
			is >> *mean;
		}
		is >> *process;
		if (rawNstate == Naa)	{
			is >> *exprelrate;
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			is >> *exprelrate;
			// is >> *relrate;
		}
		else if (mutmodel == 2)	{
			is >> *tsrelrate;
			is >> *tvrelrate;
		}
		if (gc == 3)	{
			is >> *rootgc1;
			is >> *rootgc2;
			is >> *rootgc3;
		}
		if (rootgc)	{
			is >> *rootgc;
		}
		if (freestationary)	{
			is >> *freestationary;
		}
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
		if (gammamixtree)	{
			is >> *gammamixtree;
		}
		if (gammatree)	{
			is >> *MixAlpha;
			is >> *gammatree;
		}
		if (withtimeline)	{
			is >> *timelinesigma;
			is >> *timeline;
		}
		if (whitenoise)	{
			is >> *wnvar;
			is >> *wntree;
		}
		if (jitter == 1)	{
			is >> *ugam;
			is >> *ugamtree;
		}
		if (jitter == 2)	{
			is >> *ugam;
			is >> *wngamtree;
		}
		if (contjitter)	{
			is >> *leafvar;
		}

		if (separatesyn)	{
			is >> *synsigma;
			is >> *lognormalsyntree;
		}
		if (separateomega)	{
			is >> *omegasigma;
			is >> *lognormalomegatree;
		}

		is >> *DiagArray;
		is >> *sigmaarray;
		is >> *driftarray;
		if (withexpdrift)	{
			is >> *driftphiarray;
		}
		if (withreldrift)	{
			is >> *driftphiarray2;
			is >> *driftarray2;
		}
		if (autoregressive)	{
			is >> *phi;
			is >> *mean;
		}
		is >> *process;
		if (rawNstate == Naa)	{
			is >> *exprelrate;
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			is >> *exprelrate;
			// is >> *relrate;
		}
		else if (mutmodel == 2)	{
			is >> *tsrelrate;
			is >> *tvrelrate;
		}
		if (gc == 3)	{
			is >> *rootgc1;
			is >> *rootgc2;
			is >> *rootgc3;
		}
		if (rootgc)	{
			is >> *rootgc;
		}
		if (freestationary)	{
			is >> *freestationary;
		}
	}
	*/

	void GetNucMatrix(ifstream& is)	{
		if (rawNstate == Naa)	{
			cerr << "error in get nuc matrix\n";
			exit(1);
		}
		is >> *renormrelrate;
		// is >> *relrate;
		is >> *stationary;
	}

	double GetCorrelStat(SequenceAlignment* datacopy)	{

		double** empfreq = new double*[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			empfreq[i] = new double[GetNstate()];
		}

		datacopy->GetTaxEmpiricalFreq(empfreq);

		double meangc = 0;
		double meancont = 0;
		double vargc = 0;
		double varcont = 0;
		double cov = 0;
		int weight = 0;
		for (int i=0; i<GetNtaxa(); i++)	{
			if (contdata->GetTaxonSet()->GetTaxonIndex(GetData()->GetTaxonSet()->GetTaxon(i)) != -1)	{
				double tmp = empfreq[i][1] + empfreq[i][2];
				double gc = log(tmp / (1 - tmp));
				double x = log(contdata->GetState(GetData()->GetTaxonSet()->GetTaxon(i),0));
				meangc += gc;
				vargc += gc * gc;
				meancont += x;
				varcont += x * x;
				cov += gc * x;
				weight ++;
			}
		}
		meangc /= weight;
		vargc /= weight;
		vargc -= meangc * meangc;
		meancont /= weight;
		varcont /= weight;
		varcont -= meancont * meancont;
		cov /= weight;
		cov -= meangc * meancont;
		// double correl = cov;
		double correl = cov / sqrt(varcont * vargc);

		for (int i=0; i<GetNtaxa(); i++)	{
			delete[] empfreq[i];
		}
		delete[] empfreq;

		return correl;
	}

	double GetObsCorrelStat()	{
		return GetCorrelStat(GetData());
	}

	double PostPredCompo()	{

		if (GetNstate() != Nnuc)	{
			cerr << "error: only nucleotide data\n";
			exit(1);
		}

		SequenceAlignment* datacopy = new SequenceAlignment(GetData());
		phyloprocess->PostPredSample();
		phyloprocess->GetLeafData(datacopy);

		double correl = GetCorrelStat(datacopy);

		// delete datacopy;

		return correl;
	}

	/*
	void PostPredCompo(string name)	{

		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample();
		phyloprocess->GetLeafData(datacopy);
		ostringstream s;
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
	}
	*/

	void PostPredAli(string name)	{

		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample();
		phyloprocess->GetLeafData(datacopy);
		ostringstream s;
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
	}

	void PostPredSimu(string name)	{
		// set diag elements

		if (clampsuffstat)	{
			cerr << "error in post pred: not available under clamped suff stat\n";
			exit(1);
		}
		cerr << "check post pred simu\n";
		exit(1);

		for (int i=0; i<Ncont; i++)	{
			process->Unclamp(i+L);
		}
		process->GetNodeVal(process->GetRoot()->GetNode())->Clamp();

		process->Sample();
		if (gc == 3)	{
			rootgc1->Sample();
			rootgc2->Sample();
			rootgc3->Sample();
		}
		if (rootgc)	{
			rootgc->Sample();
		}
		if (freestationary)	{
			freestationary->Sample();
		}
		if (rawNstate == Naa)	{
			exprelrate->Sample();
		}
		else if ((! mutmodel) || (mutmodel >= 4))	{
		// if (! mutmodel)	{
			exprelrate->Sample();
			// relrate->Sample();
		}
		else if (mutmodel == 2)	{
			tsrelrate->Sample();
			tvrelrate->Sample();
		}
		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample();
		phyloprocess->GetLeafData(datacopy);
		ostringstream s;
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
		cerr << "COVARIANCE : \n";
		cerr << *(GetCovMatrix()) << '\n';
		cerr << "MAX dS : " << GetMaxdS() << '\n';
		// cerr << "MAX dN : " << GetMaxdN() << '\n';
		cerr << "root : " << *process->GetNodeVal(process->GetRoot()->GetNode()) << '\n';
		if (contdata)	{
			ContinuousData* contdatacopy = new ContinuousData(contdata);
			process->GetLeafData(contdatacopy,2);
			ofstream cos((name + ".cont").c_str());
			contdatacopy->ToStream(cos);
			delete contdatacopy;
		}
		ofstream pros((name + ".param").c_str());
		ToStream(pros);
		delete datacopy;
	}

	void Simulate(string name)	{
		// set diag elements

		if (clampsuffstat)	{
			cerr << "error in simu: not available under clamped suff stat\n";
			exit(1);
		}
		MeanChronogram* meanchrono = new MeanChronogram(GetTree());
		meanchrono->Add(GetChronogram());
		meanchrono->Normalise();
		ofstream tos((name + ".tree").c_str());
		meanchrono->ToStream(tos);

		cerr << "sigma\n";
		sigmaarray->Sample();
		if (autoregressive)	{
			cerr << "??? autoregressive \n";
			exit(1);

			phi->Sample();
			// phi->setval(1);
			mean->Sample();
			/*
			for (int k=0; k<Ncont+L; k++)	{
				(*mean)[k] = 0;
			}
			*/
		}
		cerr << "unclamp leaves\n";
		for (int i=0; i<Ncont; i++)	{
			process->Unclamp(i+L);
		}
		cerr << "clamp root\n";
		process->ClampRoot();
		cerr << "sample process\n";
		process->Sample();

		Update();

		cerr << "new codon data\n";
		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);

		cerr << "ppred sample\n";
		phyloprocess->PostPredSample();

		cerr << "get leaf data\n";

		phyloprocess->GetLeafData(datacopy);
		ofstream os((name + ".ali").c_str());

		cerr << "to stream\n";
		datacopy->ToStream(os);
		if (contdata)	{
			cerr << "cont data\n";
			ContinuousData* contdatacopy = new ContinuousData(contdata);
			process->GetLeafData(contdatacopy,2);
			ofstream cos((name + ".cont").c_str());
			contdatacopy->ToStream(cos);

			cerr << "delete cont data copy\n";
			delete contdatacopy;
		}
		cerr << "param\n";
		ofstream pros((name + ".param").c_str());
		ToStream(pros);

		cerr << "delete data copy\n";
		delete datacopy;
	}

};

#endif
