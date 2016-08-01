
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSequenceAlignment.h"
#include "RYSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
// #include "CoalCalibratedChronogram.h"

#include "BranchProcess.h"
#include "MatrixTree.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"

#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"

#include "Jeffreys.h"

#include "Partition.h"
#include "SplitPartition.h"

#include "SplitLengthTree.h"
#include "SplitChronogram.h"
#include "SplitMultiVariateMove.h"
#include "MultiVariatePropagateMove.h"

#include "PartitionMultiVariateTreeProcess.h"

#include "WhiteNoise.h"

#include "LinRegContSub.h"

#include "WishartArray.h"

#include "ProteinSequenceAlignment.h"
#include "AminoAcidOmegaSubMatrix.h"
#include "SimilarityMatrix.h"
#include "AminoAcidMatrixTree.h"

#include "PosRealVectorMove.h"
#include "RenormalizedPosRealVector.h"
#include "MultiVariateRelRateCompensatoryMove.h"


class LifeHistoryRegressionModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;
	SequenceAlignment* nucdata;
	SequenceAlignment* rydata;
	CodonSequenceAlignment* codondata;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int rawNstate;
	// int Nstate;

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
	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;
	// 2 : bd with cauchy proper lower bounds
	// 3 : rbd with cauchy proper lower bounds

	double softa;

	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	GammaTree* syngammatree;

	LengthTree* lengthtree;
	NodeBranchVarTree<PosReal,PosReal>* finalchrono;

	// if different branches have different scaling factors
	BranchPartition* partition;
	SplitBranchPartition* splitpartition;
	GammaMixTree* gammamixtree;
	GammaTree* gammatree;
	Jeffreys* MixAlpha;

	string mix;

	int Nmat;

	// simplified version of the model: partition-specific drift coefficients
	// but all the rest (contsigma and reg) is constant across partitions
	int Nmat2;

	JeffreysIIDArray* ContDiagArray;
	JeffreysIIDArray* ContDiagArray0;
	JeffreysIIDArray* SubDiagArray;
	JeffreysIIDArray* SubDiagArray0;
	SigmaZero* ContSigmaZero;
	SigmaZero* ContSigmaZero0;
	SigmaZero* SubSigmaZero;
	SigmaZero* SubSigmaZero0;

	IIDArray<CovMatrix>* contsigmaarray;
	IIDArray<CovMatrix>* subsigmaarray;

	BidimIIDNormalArray* regarray;
	Gamma* phi;

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

	MultiVariateTreeProcess* contprocess;
	LinRegNormalProcess* subprocess;

	LengthTree* synratetree;
	// MeanExpTreeFromMultiVariate* synratetree;
	MeanExpTreeFromMultiVariate* nonsynratetree;
	MeanExpTreeFromMultiVariate* kappatree;
	MeanExpTreeFromMultiVariate* kappatree1;
	MeanExpTreeFromMultiVariate* kappatree2;
	MeanExpTreeFromMultiVariate* kappatree3;
	BranchValPtrTree<Dvar<PosReal> >* omegatree;

	// nucleotide mutation matrix is relrate * stationary
	IIDExp* exprelrate;
	RenormalizedPosRealVector* renormrelrate;
	Dirichlet* relrate;
	Dirichlet* tsrelrate;
	Dirichlet* tvrelrate;

	// for homogeneous model
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	int matrixtype;

	// for gc model

	bool whitenoise;
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
	int noreg;

	bool clamproot;
	bool clamptree;

	bool meanexp;

	int omegaindex;
	int omegaratiotree;
	// 0 : rate model
	// 1 : dS omega
	// 2 : dS dN

	bool noregomega;
	bool noregsyn;

	// if true: gc model
	int gc;

	int ry;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int conjpath;
	bool priorsampling;

	bool normalise;

	int nrep;

	int mutmodel;
	// 0 : gtr
	// 1 : hky with variations of Ts/Tv along the tree
	int tstvindex;
	int gcindex;


	int df;
	int Ninterpol;

	double maxdrift;
	double maxphi;
	int uniformprior;

	string bounds;
	bool withdrift;
	bool withexpdrift;
	bool withreldrift;

	int clampsuffstat;
	string suffstatfile;

	bool autoregressive;

	int krkctype;
	SimilarityMatrix * aasimilarityMatrix ;

	bool iscalib;

	public:

	SequenceAlignment* GetData()	{
		if (omegaratiotree && (rawNstate == Nnuc))	{
			cerr << "codon\n";
			return codondata;
		}
		if (rydata)	{
			cerr << "ry\n";
			return rydata;
		}
		return nucdata;
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

	LifeHistoryRegressionModel() {
	}

	BranchPartition* GetPartition()	{
		if (partition && Split())	{
			return splitpartition;
		}
		return partition;
	}

	bool isCalibrated()	{
		return iscalib;
	}

	LifeHistoryRegressionModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double insofta, double priorsigma, int indf, int inmutmodel, int ingc, int inry, bool inclampdiag, int innoreg, bool inautoregressive, int inconjpath, int contdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, string inbounds, string inmix, int inNinterpol, int inwithdrift, int inuniformprior, string insuffstatfile, string rootfile, int inkrkctype, bool sample=true, GeneticCodeType type=Universal)	{

		krkctype = inkrkctype;

		maxphi = 10.0;
		maxdrift = 100.0;
		uniformprior = inuniformprior;
		ry = inry;

		autoregressive = inautoregressive;
		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

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

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		softa = insofta;

		clampdiag = inclampdiag;
		noreg = innoreg;
		if (noreg == 2)	{
			noreg = 0;
			noregsyn = true;
		}
		if (noreg == 3)	{
			noreg = 0;
			noregomega = true;
		}
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		if (omegaratiotree > 2)	{
			noregomega = true;
			omegaratiotree -= 3;
		}
		omegaratiotree = inomegaratiotree;

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

		if (omegaratiotree)	{
			L = 2 + gc;
		}
		else	{
			L = 1 + gc;
		}
		if (mutmodel == 3)	{
			L += 3;
		}
		else if (mutmodel)	{
			L++;
		}

		// get data from file

		cerr << "data\n";
		nucdata = new FileSequenceAlignment(datafile);

		cerr << "data ok\n";

		if (ry)	{
			rydata = new RYSequenceAlignment(nucdata);
			// ofstream os("ry");
			// rydata->ToStream(os);
		}
		else	{
			rydata = 0;
		}

		rawNstate = nucdata->GetNstate();
		cerr << rawNstate << '\n';
		cerr << omegaratiotree  << '\n';
		if (omegaratiotree && (rawNstate == Nnuc))	{
			if (ry)	{
				cerr << "error : ry coding not compatible with codon\n";
				exit(1);
			}
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		cerr << nucdata->GetNsite() << '\n';
		cerr << nucdata << '\t' << GetData() << '\n';
		Nsite = GetData()->GetNsite();
		cerr << "Nsite: " << Nsite << '\n';

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
		matrixtype = 0;
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		cerr << "before : " << tree->GetFullSize(tree->GetRoot()) << '\n';
		cerr << Ninterpol << '\n';
		if (Split())	{
			splittree = new SplitTree(tree,Ninterpol);
		}
		else	{
			splittree = tree;
		}
		cerr << "after  : " << splittree->GetFullSize(splittree->GetRoot()) << '\n';

		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			cerr << "error : life history model requires continuous data\n";
			exit(1);
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
			finalchrono = 0;
		}
		else	{
			PriorMu = new Const<PosReal>(1);
			mu = new Gamma(One,PriorMu);
			mu->ClampAt(1);

			if (calibfile != "None")	{
				iscalib = true;
				double a = rootage * rootage / rootstdev / rootstdev;
				double b = rootage / rootstdev / rootstdev;
				CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

				if (chronoprior == 0)	{
					chronogram = new CalibratedChronogram(tree,mu,a,b,calibset);
				}
				else if (chronoprior == 5)	{
					cerr << "coal deprecated\n";
					exit(1);
					// chronogram = new CoalCalibratedChronogram(tree,mu,coald,coalr,N0,N1,N2,T0,calibset);
				}
				else {
					cerr << "BD\n";
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

			if (Split())	{
				finalchrono = new SplitChronogram(chronogram,GetSplitTree());
				lengthtree = finalchrono;
				// lengthtree = new SplitLengthTree(chronogram,GetSplitTree());
			}
			else	{
				finalchrono = chronogram;
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
			// true means: partitions other than the default one are recursively propagated down to the tips
			partition = new BranchPartition(tree,mix,true);
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

		ContDiagArray = new JeffreysIIDArray(Ncont,mindiag,maxdiag,Zero);
		SubDiagArray = new JeffreysIIDArray(L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			ContDiagArray->setval(1.0);
			SubDiagArray->setval(1.0);
		}
		else	{
			ContDiagArray->ClampAt(priorsigma);
			SubDiagArray->ClampAt(priorsigma);
		}

		ContDiagArray0 = new JeffreysIIDArray(Ncont,mindiag,maxdiag,Zero);
		SubDiagArray0 = new JeffreysIIDArray(L,mindiag,maxdiag,Zero);
		ContDiagArray0->ClampAt(1.0);
		SubDiagArray0->ClampAt(1.0);


		ContSigmaZero = new SigmaZero(ContDiagArray);
		ContSigmaZero0 = new SigmaZero(ContDiagArray0);
		SubSigmaZero = new SigmaZero(SubDiagArray);
		SubSigmaZero0 = new SigmaZero(SubDiagArray0);

		cerr << "sigmaarray\n";
		Nmat2 = 1;
		if (clampdiag)	{
			contsigmaarray = new DiagonalWishartArray(Nmat2, df, ContSigmaZero, ContSigmaZero0);
		}
		else	{
			contsigmaarray = new InverseWishartArray(Nmat2, df, ContSigmaZero, ContSigmaZero0);
		}
		subsigmaarray = new DiagonalWishartArray(Nmat2,df,SubSigmaZero,SubSigmaZero0);

		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(Ncont));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont));
			ifstream is(rootfile.c_str());
			for (int i=0; i<Ncont; i++)	{
				double mean, var;
				is >> mean >> var;
				(*rootmean)[i] = mean;
				(*rootvar)[i] = var;
			}
		}

		cerr << "drift\n";
		// driftarray = new MultiVarArray(Nmat, Zero, ContDiagArray, ContDiagArray0);
		if (uniformprior)	{
			MultiUniArray* tmparray = new MultiUniArray(Nmat, Ncont, Zero, maxdrift);
			driftarray = tmparray;
			if (! withdrift)	{
				tmparray->ClampAtZero();
			}
		}
		else	{
			MultiVarArray* tmparray = new MultiVarArray(Nmat, Zero, ContDiagArray, ContDiagArray0);
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
				driftarray2 = new MultiUniArray(Nmat, Ncont, Zero, maxdrift);
				driftphiarray2 = new PosUniIIDArray(Nmat,One,maxphi);
			}
			else	{
				driftarray2 = new MultiVarArray(Nmat, Zero, ContDiagArray, ContDiagArray0);
				driftphiarray2 = new GammaIIDArray(Nmat,One,One);
			}
			drift2 = driftarray2->GetVal(0);
			driftphi2 = driftphiarray2->GetVal(0);
		}

		cerr << "regarray\n";
		if (autoregressive)	{
			phi = new Gamma(One,One);
			regarray = new BidimIIDNormalArray(Nmat2,L,Ncont+1,Zero,One,One);
		}
		else	{
			phi = 0;
			regarray = new BidimIIDNormalArray(Nmat2,L,Ncont,Zero,One,One);
		}
		regarray->SetAtZero();
		if (noreg)	{
			regarray->ClampAtZero();
		}
		if (noregomega)	{
			for (int i=0; i<Nmat2; i++)	{
				regarray->GetIIDNormal(i,1)->ClampAtZero();
			}
		}

		if (noregsyn)	{
			for (int i=0; i<Nmat2; i++)	{
				regarray->GetIIDNormal(i,0)->ClampAtZero();
			}
		}

		if (gammatree)	{
			contprocess = new PartitionMultiVariateTreeProcess(contsigmaarray,lengthtree,GetPartition(),gammatree,driftarray,driftphiarray,finalchrono, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);

		}
		else	{
			contprocess = new PartitionMultiVariateTreeProcess(contsigmaarray,lengthtree,GetPartition(),gammamixtree,driftarray,driftphiarray,finalchrono, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);
		}

		cerr << "lin reg array\n";
		subprocess = new PartitionLinRegNormalProcess(lengthtree,contprocess,regarray,subsigmaarray,GetPartition(),phi);
		cerr << "ok\n";

		if (whitenoise)	{
			wnvar = new Gamma(One,One);
			wntree = new WhiteNoiseProcess(lengthtree,Zero,wnvar);
		}
		else	{
			wnvar = 0;
			wntree = 0;
		}

		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i,i,contdatatype);
			}
		}

		cerr << "cutoff\n";

		for (int l=0; l<L; l++)	{
			subprocess->CutOff(1,l);
		}

		if (bounds != "None")	{
			FileBoundSet boundset(bounds,contprocess->GetTree());
			contprocess->SetBounds(boundset,0);
		}
		cerr << "ok\n";

		CreateSubstitutionProcess();

		if (phyloprocess)	{
			phyloprocess->Unfold();
		}
		if (sample)	{
			if (Split())	{
				cerr << "before\n";
				SplitLengthTree* tmp = dynamic_cast<SplitLengthTree*> (lengthtree);
				tmp->Check();
			}

			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		// register model
		RootRegister(PriorMu);
		RootRegister(Zero);
		RootRegister(One);
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (exprelrate)	{
			RootRegister(exprelrate);
		}
		if (relrate)	{
			RootRegister(relrate);
		}
		else if (tsrelrate)	{
			RootRegister(tsrelrate);
			RootRegister(tvrelrate);
		}
		if (stationary)	{
			RootRegister(stationary);
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
		if (sample)	{
			Update();
		}

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~LifeHistoryRegressionModel() {}

	void CreateSubstitutionProcess()	{

		matrixtype = 0;

		int offset = 0;
		if (Unconstrained())	{
			offset = 1;
		}
		else	{
			synratetree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,0,INTEGRAL,false,meanexp);
			nonsynratetree = 0;
		}
		omegaindex = 1-offset;

		if (omegaratiotree == 0)	{
			omegatree = 0;
		}
		else if (omegaratiotree == 1)	{
			omegatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,1-offset,MEAN,false,meanexp);
		}
		else	{
			nonsynratetree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,1-offset,INTEGRAL,false,meanexp);
			omegatree = new RatioTree(nonsynratetree,synratetree);
		}

		relrate = 0;
		tsrelrate = 0;
		tvrelrate = 0;
		exprelrate = 0;
		renormrelrate = 0;

		stationary = 0;
		rootgc = 0;

		if (rawNstate == 20)	{
			exprelrate = new IIDExp(Naa*(Naa-1)/2);
			renormrelrate = new RenormalizedPosRealVector(exprelrate);
			LinRegNormal* rootval = subprocess->GetLinRegNormal(subprocess->GetTree()->GetRoot());
			rootval->ClampAt(0,1-offset);
		}
		else if (mutmodel == 0)	{
			if (ry)	{
				relrate = new Dirichlet(1);
			}
			else	{
				relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
			}
		}
		else if (mutmodel == 2)	{
			tsrelrate = new Dirichlet(2);
			tvrelrate = new Dirichlet(4);
		}

		if (! omegaratiotree)	{
			if (mutmodel)	{
				if (gc)	{
					if (whitenoise)	{
						cerr << "whitenoise: still to be implemented\n";
						exit(1);
						// gctree = new MeanJitteredLogitTree(lengthtree,subprocess,1-offset,wntree,false);
					}
					else	{
						// gctree = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,1-offset,MEAN,false);
					}
					gcindex = 1-offset;
					kappatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false,meanexp);
					tstvindex = 2-offset;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new HKYGCNucMatrixTree(kappatree,stattree,tsrelrate,tvrelrate,normalise);
					nucmatrix = 0;
					stationary = 0;
				}
				else	{
					stationary = new Dirichlet(Nnuc);
					kappatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,1-offset,MEAN,false,meanexp);
					tstvindex = 1-offset;
					nucmatrixtree = new HKYNucMatrixTree(kappatree,stationary,tsrelrate,tvrelrate,normalise);
					nucmatrix = 0;
					gctree = 0;
					rootgc = 0;
					stattree = 0;
				}
			}
			else	{
				if (gc)	{
					if (whitenoise)	{
						cerr << "whitenoise: still to be implemented\n";
						exit(1);
						// gctree = new MeanJitteredLogitTree(lengthtree,subprocess,1-offset,wntree,false);
					}
					else	{
						// gctree = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,1-offset,MEAN,false);
					}
					gcindex = 1-offset;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
					nucmatrix = 0;
					stationary = 0;
				}
				else	{
					if (ry)	{
						stationary = new Dirichlet(2);
					}
					else	{
						stationary = new Dirichlet(Nnuc);
					}
					nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
					gctree = 0;
					rootgc = 0;
					stattree = 0;
					nucmatrixtree = 0;
				}
			}

			// make substitution mappings
			if (conjpath)	{
				if (gc || mutmodel)	{
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, nucmatrixtree, GetData());
				}
				else	{
					pathconjtree = new OneMatrixPathConjugateTree(synratetree,nucmatrix,GetData());
				}

				if (clampsuffstat)	{
					cerr << "read suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
				pathconjtree = 0;
				if (priorsampling)	{
					phyloprocess = 0;
				}
				else if (gc || mutmodel)	{
					phyloprocess = new BranchMatrixPhyloProcess(synratetree, nucmatrixtree, GetData());
				}
				else	{
					phyloprocess = new OneMatrixPhyloProcess(synratetree, nucmatrix, GetData());
				}
			}
		}
		else if (rawNstate == Naa)	{

			cerr << "amino acid\n";
			kappatree = 0;
			nucmatrix = 0;

			cerr << "Using the model based on " << krkctype << endl;
			switch(krkctype){
				case 0 : aasimilarityMatrix = new PolarityBasedMatrix();break;
				case 1 : aasimilarityMatrix = new VolumeBasedMatrix(); break;
				case 2 : aasimilarityMatrix = new ChargeBasedMatrix(); break;
				case 3 : aasimilarityMatrix = new PolarityAndVolumeBasedMatrix(); break;
				case 4 : cerr<< " NO Model Specfified" << endl; exit(0);
			}

			aasimilarityMatrix->Affiche();
			cerr << "matrix\n";
			stationary = new Dirichlet(Naa);
			matrixtree = new AminoAcidOmegaMatrixTree(aasimilarityMatrix,renormrelrate, stationary, omegatree, One);

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
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree (synratetree, matrixtree, nucdata);

				if (clampsuffstat)	{
					cerr << "read suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
				pathconjtree = 0;
				if (priorsampling)	{
					phyloprocess = 0;
				}
				phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, nucdata);
			}
			cerr << "aa ok\n";
		}
		else	{

			matrixtype = 1;
			if (mutmodel)	{
				if (gc == 3)	{
					matrixtype = 2;
					if (mutmodel == 3)	{
						gctree1 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false);
						gctree2 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,3-offset,MEAN,false);
						gctree3 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,4-offset,MEAN,false);
						gcindex = 2-offset;
						kappatree1 = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,5-offset,MEAN,false,meanexp);
						kappatree2 = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,6-offset,MEAN,false,meanexp);
						kappatree3 = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,7-offset,MEAN,false,meanexp);
						kappatree = 0;
						tstvindex = 5-offset;
						rootgc1 = new Beta(One,One);
						rootgc2 = new Beta(One,One);
						rootgc3 = new Beta(One,One);
						stattree1 = new GCStatTree(gctree1,rootgc1);
						stattree2 = new GCStatTree(gctree2,rootgc2);
						stattree3 = new GCStatTree(gctree3,rootgc3);
						nucmatrixtree1 = new HKYGCNucMatrixTree(kappatree1,stattree1,tsrelrate,tvrelrate,normalise);
						nucmatrixtree2 = new HKYGCNucMatrixTree(kappatree2,stattree2,tsrelrate,tvrelrate,normalise);
						nucmatrixtree3 = new HKYGCNucMatrixTree(kappatree3,stattree3,tsrelrate,tvrelrate,normalise);
						matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);

						gctree = 0;
						rootgc = 0;
						stattree = 0;
						nucmatrixtree = 0;
						nucmatrix = 0;
						stationary = 0;
					}
					else	{
						gctree1 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false);
						gctree2 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,3-offset,MEAN,false);
						gctree3 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,4-offset,MEAN,false);
						gcindex = 2-offset;
						kappatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,5-offset,MEAN,false,meanexp);
						tstvindex = 5-offset;
						rootgc1 = new Beta(One,One);
						rootgc2 = new Beta(One,One);
						rootgc3 = new Beta(One,One);
						stattree1 = new GCStatTree(gctree1,rootgc1);
						stattree2 = new GCStatTree(gctree2,rootgc2);
						stattree3 = new GCStatTree(gctree3,rootgc3);
						nucmatrixtree1 = new HKYGCNucMatrixTree(kappatree,stattree1,tsrelrate,tvrelrate,normalise);
						nucmatrixtree2 = new HKYGCNucMatrixTree(kappatree,stattree2,tsrelrate,tvrelrate,normalise);
						nucmatrixtree3 = new HKYGCNucMatrixTree(kappatree,stattree3,tsrelrate,tvrelrate,normalise);
						matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);

						gctree = 0;
						rootgc = 0;
						stattree = 0;
						nucmatrixtree = 0;
						nucmatrix = 0;
						stationary = 0;
					}
				}
				else if (gc == 1)	{
					gctree = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false);
					gcindex = 2-offset;
					kappatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,3-offset,MEAN,false,meanexp);
					tstvindex = 3-offset;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new HKYGCNucMatrixTree(kappatree,stattree,tsrelrate,tvrelrate,normalise);
					matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);

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
				}
				else	{
					stationary = new Dirichlet(Nnuc);
					nucmatrix = 0;
					kappatree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false,meanexp);
					tstvindex = 2-offset;
					nucmatrixtree = new HKYNucMatrixTree(kappatree,stationary,tsrelrate,tvrelrate,normalise);
					matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);

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
				}
			}
			else	{
				if (gc == 3)	{
					matrixtype = 2;
					kappatree = 0;

					gctree1 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false);
					gctree2 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,3-offset,MEAN,false);
					gctree3 = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,4-offset,MEAN,false);
					gcindex = 2-offset;
					rootgc1 = new Beta(One,One);
					rootgc2 = new Beta(One,One);
					rootgc3 = new Beta(One,One);
					stattree1 = new GCStatTree(gctree1,rootgc1);
					stattree2 = new GCStatTree(gctree2,rootgc2);
					stattree3 = new GCStatTree(gctree3,rootgc3);
					nucmatrixtree1 = new GTRGCNucMatrixTree(relrate,stattree1,normalise);
					nucmatrixtree2 = new GTRGCNucMatrixTree(relrate,stattree2,normalise);
					nucmatrixtree3 = new GTRGCNucMatrixTree(relrate,stattree3,normalise);
					matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);

					gctree = 0;
					rootgc = 0;
					stattree = 0;
					nucmatrixtree = 0;
					nucmatrix = 0;
					stationary = 0;
				}
				else if (gc == 1)	{
					kappatree = 0;

					gctree = new MeanLogitTreeFromMultiVariate(lengthtree,subprocess,2-offset,MEAN,false);
					gcindex = 2-offset;
					rootgc = new Beta(One,One);
					stattree = new GCStatTree(gctree,rootgc);
					nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
					matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);

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
				}
				else	{
					kappatree = 0;

					stationary = new Dirichlet(Nnuc);
					nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
					matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

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
				}
			}

			// make substitution mappings
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
				if (clampsuffstat)	{
					pathconjtree->ReadFromFile(suffstatfile);
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else if (conjpath == 1)	{
				pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				if (clampsuffstat)	{
					cerr << "read suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
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

	Tree* GetTree() {return tree;}
	Tree* GetFineGrainedTree() {return splittree;}

	int GetTsTvIndex()	{
		if (! mutmodel)	{
			cerr << "error in get tstv index\n";
			exit(1);
		}
		return tstvindex;
	}

	void UpdateOmegaTree()	{
		if (omegaratiotree == 2)	{
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

	BranchVarTree<UnitReal>* GetGCTree() {return gctree;}
	MeanLogitTreeFromMultiVariate* GetGC1Tree() {return gctree1;}
	MeanLogitTreeFromMultiVariate* GetGC2Tree() {return gctree2;}
	MeanLogitTreeFromMultiVariate* GetGC3Tree() {return gctree3;}

	bool isGCActivated() {return GetGCTree();}
	bool isGC3Activated() {return GetGC3Tree();}
	int GetGCIndex() {return gcindex;}

	MultiVariateTreeProcess* GetContProcess() {return contprocess;}
	NodeVarTree<RealVector>* GetSubProcess() {return subprocess;}
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

	CovMatrix* GetContMatrix(int mat = 0) {return contsigmaarray->GetVal(mat);}

	CalibratedChronogram* GetCalibratedChronogram()	{
		if (Unconstrained())	{
			cerr << "error : calibrated chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double GetPhi()	{
		if (phi)	{
			return phi->val();
		}
		return 0;
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
			if ((chronoprior >= 1) && (chronoprior <= 4))	{
				total += Chi->GetLogProb();
				total += Chi2->GetLogProb();
			}
			total += chronogram->GetLogProb();
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

		total += ContDiagArray->GetLogProb();
		total += SubDiagArray->GetLogProb();
		total += contsigmaarray->GetLogProb();
		total += subsigmaarray->GetLogProb();
		total += regarray->GetLogProb();
		if (autoregressive)	{
			total += phi->GetLogProb();
		}
		total += driftarray->GetLogProb();
		if (withexpdrift)	{
			total += driftphiarray->GetLogProb();
		}
		if (withreldrift)	{
			total += driftphiarray2->GetLogProb();
			total += driftarray2->GetLogProb();
		}
		total += contprocess->GetLogProb();
		total += subprocess->GetLogProb();

		if (relrate)	{
			total += relrate->GetLogProb();
		}
		if (exprelrate)	{
			total += exprelrate->GetLogProb();
		}
		if (tsrelrate)	{
			total += tsrelrate->GetLogProb();
			total += tvrelrate->GetLogProb();
		}

		if (gc == 3)	{
			total += rootgc1->GetLogProb();
			total += rootgc2->GetLogProb();
			total += rootgc3->GetLogProb();
		}
		else if (gc)	{
			total += rootgc->GetLogProb();
		}
		else	{
			total += stationary->GetLogProb();
		}
		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (priorsampling)	{
			return 0;
		}
		else if (clampsuffstat)	{
			ret = pathconjtree->GetLogProb();
		}
		else	{
			ret = phyloprocess->GetLogProb();
		}
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			if (! clampsuffstat)	{
				scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
			}
		}
		else	{
			if (phyloprocess)	{
				scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
			}
		}

		vector <SplitMultiVariateNodeMove*> nodesplitarray;
		vector <SplitMultiVariateBranchMove*> branchsplitarray;
		if (Split())	{
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,10));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,0.1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,0.01));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,10));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,0.1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,0.01));
		}

		for (int i=0; i<nrep; i++)	{
			if (Unconstrained())	{
				scheduler.Register(new SimpleMove(mu,1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(mu,0.1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(syngammatree,1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.01),10,"syngamtree");
			}
			else if (! clamptree)	{
				if ((chronoprior >= 1) && (chronoprior <= 4))	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
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
			if (isCalibrated() && (! clamptree))	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
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

			scheduler.Register(new SimpleMove(contprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(subprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(subprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(subprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(subprocess,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(contsigmaarray,10),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigmaarray,1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigmaarray,0.1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigmaarray,0.01),100,"cont sigma");

			scheduler.Register(new SimpleMove(ContDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,0.1),10,"theta");

			scheduler.Register(new SimpleMove(regarray,10),100,"reg");
			scheduler.Register(new SimpleMove(regarray,1),100,"reg");
			scheduler.Register(new SimpleMove(regarray,0.1),100,"reg");
			scheduler.Register(new SimpleMove(regarray,0.01),100,"reg");

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");
			}

			scheduler.Register(new SimpleMove(subsigmaarray,10),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigmaarray,1),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigmaarray,0.1),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigmaarray,0.01),100,"sub sigma");

			scheduler.Register(new SimpleMove(SubDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(SubDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(SubDiagArray,0.1),10,"theta");

			if (rawNstate == Naa)	{
				scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.1,2),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.03,4),10,"relrates");
				scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");

				scheduler.Register(new LinRegRelRateCompensatoryMove(subprocess,exprelrate,aasimilarityMatrix,0.3,omegaindex),20,"process relrate comp move");
				scheduler.Register(new LinRegRelRateCompensatoryMove(subprocess,exprelrate,aasimilarityMatrix,0.03,omegaindex),20,"process relrate comp move");

			}
			if (relrate)	{
				scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
				scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
				scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
			}
			if (tsrelrate)	{
				scheduler.Register(new ProfileMove(tsrelrate,0.1,1),10,"relrates");
				scheduler.Register(new ProfileMove(tvrelrate,0.1,1),10,"relrates");
				scheduler.Register(new ProfileMove(tsrelrate,0.03,2),10,"relrates");
				scheduler.Register(new ProfileMove(tvrelrate,0.03,2),10,"relrates");
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
			else if (gc)	{
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
		}
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";

		if (Unconstrained())	{
			mu->Sample();
			syngammatree->Sample();
		}
		else	{
			if ((chronoprior >= 1) && (chronoprior <= 4))	{
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

		ContDiagArray->Sample();
		SubDiagArray->Sample();

		contsigmaarray->Sample();
		subsigmaarray->Sample();
		regarray->Sample();
		phi->Sample();
		// sigma->SetIdentity();
		driftarray->Sample();
		if (withexpdrift)	{
			driftphiarray->Sample();
		}
		if (withreldrift)	{
			driftphiarray2->Sample();
			driftarray2->Sample();
		}

		contprocess->Sample();
		subprocess->Sample();

		for (int l=0; l<L; l++)	{
			subprocess->CutOff(1,l);
		}
		if (! Unconstrained())	{
			GetSynRateTree()->specialUpdate();
		}
		if (omegaratiotree)	{
			UpdateOmegaTree();
		}

		if (exprelrate)	{
			exprelrate->Sample();
		}
		if (relrate)	{
			relrate->Sample();
		}
		if (tsrelrate)	{
			tsrelrate->Sample();
			tvrelrate->Sample();
		}

		if (gc == 3)	{
			rootgc1->Sample();
			rootgc2->Sample();
			rootgc3->Sample();
		}
		else if (gc)	{
			rootgc->Sample();
		}
		else	{
			stationary->Sample();
		}

		if (phyloprocess)	{
			phyloprocess->Sample();
		}

		cerr << "ok\n";
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
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
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
		return GetSynRateTree()->GetTotal();
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
		if (omegaratiotree)	{
			return ((RatioTree*) omegatree)->GetMean();
		}
		else	{
			return ((MeanExpTreeFromMultiVariate*) omegatree)->GetMean();
		}
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

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		if (omegaratiotree)	{
			os << "\tsynrate\tomega";
		}
		else	{
			os << "\trate";
		}

		// os << "\trootleft\trootright";

		for (int mat=0; mat<Nmat2; mat++)	{
			for (int k=0; k<Ncont; k++)	{
				for (int l=k+1; l<Ncont; l++)	{
					os << '\t' << "cont_" << mat << '_' << k << '_' << l;
				}
			}
			for (int k=0; k<L; k++)	{
				for (int l=0; l<regarray->GetDim(); l++)	{
					os << '\t' << "reg_" << mat << '_' << k << '_' << l;
				}
			}
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << "cont_" << mat << '_' << k << '_' << k;
			}
			for (int k=0; k<L; k++)	{
				os << '\t' << "sub_" << mat << '_' << k << '_' << k;
			}
		}
		if (autoregressive)	{
			os << '\t' << "phi";
		}
		if (whitenoise)	{
			os << '\t' << "wnvar";
		}
		for (int mat=0; mat<Nmat; mat++)	{
			if (withdrift)	{
				os << "\tdim";
				for (int k=0; k<Ncont; k++)	{
					os << '\t' << "drift_" << mat << '_' << k;
				}
			}
			if (withexpdrift)	{
				os << '\t' << "phi_" << mat;
			}
			if (withreldrift)	{
				for (int k=0; k<Ncont; k++)	{
					os << '\t' << "drift2_" << mat << '_' << k;
				}
				os << '\t' << "phi2_" << mat;
			}
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 4))	{
			os << "\tp1\tp2";
		}

		os << "\tdim";
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}

		os << "\tdim";
		for (int k=0; k<L; k++)	{
			os << '\t' << "root_" << k;
		}

		if (gc == 3)	{
			os << "\tmeangc1\tmeangc2\tmeangc3";
			os << "\trootgc1\trootgc2\trootgc3";
		}
		else if (gc)	{
			os << "\tmeangc\trootgc";
		}
		else	{
			os << "\tstatent";
		}
		if (exprelrate || relrate)	{
			os << "\trrent";
		}
		if (tsrelrate)	{
			os << "\ttsrrent\ttvrrent";
		}
		if (kappatree) {
			os << "\tmeankappa";
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

		if ((chronoprior >= 1) && (chronoprior <= 4))	{
			os << "\tnumerror";
		}

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		if (omegaratiotree)	{
			os << '\t' << GetMeanOmega();
		}

		for (int mat=0; mat<Nmat2; mat++)	{
			for (int k=0; k<Ncont; k++)	{
				for (int l=k+1; l<Ncont; l++)	{
					os << '\t' << (*contsigmaarray->GetVal(mat))[k][l];
				}
			}
			for (int k=0; k<L; k++)	{
				for (int l=0; l<regarray->GetDim(); l++)	{
					os << '\t' << (*(*regarray->GetIIDNormalArray(mat))[k])[l];
				}
			}
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << (*contsigmaarray->GetVal(mat))[k][k];
			}
			for (int k=0; k<L; k++)	{
				os << '\t' << (*subsigmaarray->GetVal(mat))[k][k];
			}
		}
		if (autoregressive)	{
			os << '\t' << phi->val();
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
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 4))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}
		os << '\t' << GetContProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();
		os << '\t' << GetSubProcess()->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode())->val();

		if (gc == 3)	{
			os << '\t' << GetMeanGC1();
			os << '\t' << GetMeanGC2();
			os << '\t' << GetMeanGC3();
			os << '\t' << GetRootGC1();
			os << '\t' << GetRootGC2();
			os << '\t' << GetRootGC3();
		}
		else if (gc)	{
			os << '\t' << GetMeanGC();
			os << '\t' << GetRootGC();
		}
		else	{
			os << '\t' << stationary->val().GetEntropy();
		}

		if (exprelrate)	{
			os << '\t' << renormrelrate->val().GetEntropy();
		}
		if (relrate)	{
			os << '\t' << relrate->val().GetEntropy();
		}
		if (tsrelrate)	{
			os << '\t' << tsrelrate->val().GetEntropy();
			os << '\t' << tvrelrate->val().GetEntropy();
		}
		if (kappatree)	{
			os << '\t' << GetMeanKappa();
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

		if ((chronoprior >= 1) && (chronoprior <= 4))	{
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\n';
		os.flush();

	}

	void PrintEntries(ostream& os, int* array = 0)	{

		int cumul = 0;
		if ((! array) || (array[cumul] == 1))	{
			os << "dS\n";
		}
		cumul++;
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
		if (mutmodel == 3)	{
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
		for (int k=0; k<Ncont; k++)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "character " << k + 1 << '\n';
			}
			cumul++;
		}
	}

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
			if ((chronoprior >= 1) && (chronoprior <= 4))	{
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
		if (whitenoise)	{
			os << *wnvar << '\n';
			os << *wntree << '\n';
		}
		os << *ContDiagArray << '\n';
		os << *SubDiagArray << '\n';
		os << *contsigmaarray << '\n';
		os << *subsigmaarray << '\n';
		os << *regarray << '\n';
		if (autoregressive)	{
			os << *phi << '\n';
		}
		os << *driftarray << '\n';
		if (withexpdrift)	{
			os << *driftphiarray << '\n';
		}
		if (withreldrift)	{
			os << *driftphiarray2 << '\n';
			os << *driftarray2 << '\n';
		}
		os << '\n';
		os << *contprocess << '\n';
		os << *subprocess << '\n';
		if (exprelrate)	{
			os << *exprelrate << '\n';
		}
		if (relrate)	{
			os << *relrate << '\n';
		}
		if (tsrelrate)	{
			os << *tsrelrate << '\n';
			os << *tvrelrate << '\n';
		}
		if (gc == 3)	{
			os << *rootgc1 << '\n';
			os << *rootgc2 << '\n';
			os << *rootgc3 << '\n';
		}
		else if (gc)	{
			os << *rootgc << '\n';
		}
		else	{
			os << *stationary << '\n';
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
			if ((chronoprior >= 1) && (chronoprior <= 4))	{
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
		if (whitenoise)	{
			is >> *wnvar;
			is >> *wntree;
		}

		is >> *ContDiagArray;
		is >> *SubDiagArray;
		is >> *contsigmaarray;
		is >> *subsigmaarray;
		is >> *regarray;
		if (autoregressive)	{
			is >>  *phi;
		}
		is >> *driftarray;
		if (withexpdrift)	{
			is >> *driftphiarray;
		}
		if (withreldrift)	{
			is >> *driftphiarray2;
			is >> *driftarray2;
		}
		is >> *contprocess;
		is >> *subprocess;
		if (exprelrate)	{
			is >> *exprelrate;
		}
		if (relrate)	{
			is >> *relrate;
		}
		if (tsrelrate)	{
			is >> *tsrelrate;
			is >> *tvrelrate;
		}
		if (gc == 3)	{
			is >> *rootgc1;
			is >> *rootgc2;
			is >> *rootgc3;
		}
		else if (gc)	{
			is >> *rootgc;
		}
		else	{
			is >> *stationary;
		}
	}

	/*
	void GetNucMatrix(ifstream& is)	{
		is >> *relrate;
		is >> *stationary;
		cerr << *relrate << '\n';
		cerr << *stationary << '\n';
	}
	*/
};

#endif
