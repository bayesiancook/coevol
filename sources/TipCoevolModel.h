
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "IID.h"

#include "PrecisionNormalTreeProcess.h"

#include "BrownianProcess.h"
#include "CommutativeBrownianProcess.h"
#include "BrownianMove.h"

#include "BDCalibratedChronogram.h"
#include "SerialCoalescent.h"
#include "LogisticSerialCoalescent.h"
#include "SerialBirthDeath.h"

#include "F81TransitionMatrix.h"
#include "GeneralConjugatePath.h"

#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "Jeffreys.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "MeanChronogram.h"
#include "AuxCoevol.h"

#include "WhiteNoise.h"
#include "BranchProductProcess.h"
#include "MultiVariateUgamCompMove.h"

class F81MatrixTree : public BranchValPtrTree<RandomTransitionMatrix>	{

	public:

	F81MatrixTree(LengthTree* inlengthtree, Var<Profile>* instationary) {
		SetWithRoot(true);
		lengthtree = inlengthtree;
		stationary = instationary;
		RecursiveCreate(GetRoot());
	}

	~F81MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomTransitionMatrix* CreateBranchVal(const Link* link)	{
		return new RandomF81TransitionMatrix(stationary,lengthtree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	private:

	LengthTree* lengthtree;
	Var<Profile>* stationary;

};


class BranchTransitionMatrixPhyloProcess : public PhyloProcess	{


	protected:

	public:

	BranchTransitionMatrixPhyloProcess(LengthTree* inlengthtree, BranchValPtrTree<RandomTransitionMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(inlengthtree,indata,false)	{
		matrixtree = inmatrixtree;
		lengthtree = inlengthtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		RandomBranchSitePath* tmp = new RandomBranchSitePath(this, matrixtree->GetBranchVal(link->GetBranch()), 0, lengthtree->GetBranchVal(link->GetBranch()));
		return tmp;
	}

	protected:
	BranchValPtrTree<RandomTransitionMatrix>* matrixtree;
	BranchVarTree<PosReal>* lengthtree;
};

class GTRLogNormalModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	SequenceAlignment* morpho2data;
	SequenceAlignment* morpho3data;
	SequenceAlignment* morpho4data;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	int Ncont;
	int dim;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	double rootalpha;
	double rootbeta;

	// chronogram
	CalibratedChronogram* chronogram;

	int scaling;
	// 0 : no scaling
	// 1 : burst followed by exponential decay
	Var<PosReal>* scalet0;
	Const<PosReal>* scalefactorAlpha;
	Const<PosReal>* scalefactorBeta;
	Const<PosReal>* scalerateAlpha;
	Const<PosReal>* scalerateBeta;
	Gamma* scalefactor;
	Gamma* scalerate;

	GlobalScalingFunction* scalefunction;

	Gamma* scalesigma;
	GammaWhiteNoiseProcess* gammascaletree;

	// 2 : gamma white noise scaling

	CalibrationSet* calibset;

	Const<PosReal>* CoalAlpha;
	Const<PosReal>* CoalBeta;
	Gamma* coalrate;

	Const<PosReal>* T0Alpha;
	Const<PosReal>* T0Beta;
	Const<PosReal>* divrateAlpha;
	Const<PosReal>* extrateAlpha;
	Const<PosReal>* K0Alpha;
	Const<PosReal>* massextAlpha;
	Const<PosReal>* K1Alpha;
	Const<PosReal>* divrateBeta;
	Const<PosReal>* extrateBeta;
	Const<PosReal>* K0Beta;
	Const<PosReal>* massextBeta;
	Const<PosReal>* K1Beta;
	Const<PosReal>* T1;

	Gamma* T0;
	Gamma* divrate;
	Gamma* extrate;
	Gamma* K0;
	Beta* massext;
	Gamma* K1;

	Const<PosReal>* bddivrateAlpha;
	Const<PosReal>* bddivrateBeta;
	Const<PosReal>* muAlpha;
	Const<PosReal>* muBeta;
	Const<PosReal>* psiAlpha;
	Const<PosReal>* psiBeta;
	Const<PosReal>* rhoAlpha;
	Const<PosReal>* rhoBeta;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Rvar<PosReal>* Chi;
	Rvar<PosReal>* Chi2;

	Gamma* bddivrate;
	Gamma* mu;
	Gamma* psi;
	Beta* rho;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;
	ConjugateInverseWishart* conjugatesigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* process;
	// ConjugateMultiVariateTreeProcess* process;

	LengthTree* nucratetree;
	LengthTree* morphoratetree;

	BranchProductProcess* productnucratetree;
	BranchProductProcess* productmorphoratetree;

	MeanExpTreeFromMultiVariate* meanexpnucratetree;
	MeanExpTreeFromMultiVariate* meanexpmorphoratetree;

	Const<PosReal>* rateAlpha;
	Const<PosReal>* rateBeta;
	Const<PosReal>* morphorateAlpha;
	Const<PosReal>* morphorateBeta;

	Gamma* wnnucmeanrate;
	Gamma* wnnucsigma;
	Gamma* wnmorphomeanrate;
	Gamma* wnmorphosigma;

	GammaWhiteNoiseProcess* wnnucratetree;
	GammaWhiteNoiseProcess* wnmorphoratetree;

	// 0: simple log normal
	// 1 : white noise
	// 2 : combined multiplicatively
	int clockprior;

	Dirichlet* nucrelrate;
	Dirichlet* nucstationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	Dirichlet* morpho2stationary;
	Dirichlet* morpho3stationary;
	Dirichlet* morpho4stationary;

	F81MatrixTree* morpho2matrixtree;
	F81MatrixTree* morpho3matrixtree;
	F81MatrixTree* morpho4matrixtree;

	// phylo process
	// OneMatrixPhyloProcess* phyloprocess;
	PathConjugateTree* nucpathconjtree;
	PhyloProcess* nucphyloprocess;
	TransitionPathConjugateTree* morpho2pathconjtree;
	TransitionPathConjugateTree* morpho3pathconjtree;
	TransitionPathConjugateTree* morpho4pathconjtree;
	PhyloProcess* morpho2phyloprocess;
	PhyloProcess* morpho3phyloprocess;
	PhyloProcess* morpho4phyloprocess;

	bool conjpath;
	int chronoprior;

	bool clampdiag;
	bool clamptree;
	bool meanexp;

	int L;
	int df;

	int morpho;
	int prior;

	int brownian;
	BrownianProcess* brownianProcess;
	ConjugateBrownianProcess* conjugateBrownianProcess;
	int nSegments;	//Number of segments in each brownian bridge
	Segmentation segm;
	CommutativeBrownianProcess *nuccommutativeProcess;
	CommutativeBrownianProcess *morphocommutativeProcess;
	int conjsigma;
	int fixbl;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRLogNormalModel(string nucdatafile, string morpho2datafile, string morpho3datafile, string morpho4datafile, string contdatafile, int contdatatype, string treefile, string calibfile, double inrootage, double inrootstdev, double divrateval, double extrateval, double massextval, double K0val, double K1val, double divratestdev, double extratestdev, double massextstdev, double K0stdev, double K1stdev, double T1val, double bddivratemean, double mumean, double psimean, double rhomean, double bddivratestdev, double mustdev, double psistdev, double rhostdev, double divcutoff, int Nextant, int inchronoprior, int inclampdiag, int inclamptree, int inmeanexp, int indf, int inbrownian, int innSegments, Segmentation insegm, int inclockprior, double ratemean, double ratestdev, double morphoratemean, double morphoratestdev, int inscaling, double scalet0val, double scalefactorval, double scalefactorstdev, double scalerateval, double scaleratestdev, int inprior, bool sample = true)	{

        fixbl = 0;
        prior = inprior;
		/*
		if (inprior == 0)	{
			prior = 0;
			fixbl = 0;
		}
		else if (inprior == 1)	{
			prior = 1;
			fixbl = 0;
		}
		else if (inprior == 2)	{
			prior = 0;
			fixbl = 1;
		}
		*/

		conjpath = true;
		chronoprior = inchronoprior;

		int bdpriortype = 1;
		// 1: improper prior for upper bounds
		// 2: Cauchy prior for upper bounds
		if (chronoprior == 5)	{
			chronoprior = 4;
			bdpriortype = 2;
		}
		else if (chronoprior == 6)	{
			chronoprior = 4;
			bdpriortype = 3;
		}
		clampdiag = inclampdiag;
		clamptree = inclamptree;
		if (clamptree)	{
			fixbl = 1;
		}
		meanexp = inmeanexp;

		clockprior = inclockprior;

		scaling = inscaling;

		brownian = inbrownian;
		conjsigma = 1;
		// 0: branchwise approx
		// 1: brownian 

		nSegments = innSegments;
		segm = insegm;

		// fetch data from file
		nucdata = new FileSequenceAlignment(nucdatafile);

		morpho = 0;
		if (morpho2datafile != "null")	{
			morpho = 1;
		}

		if (morpho)	{
			morpho2data = new FileSequenceAlignment(morpho2datafile);
			morpho3data = new FileSequenceAlignment(morpho3datafile);
			morpho4data = new FileSequenceAlignment(morpho4datafile);
		}

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		// get continuous data from file
		if (contdatafile != "null")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		if (clockprior == 1)	{
			L = 0;
		}
		else	{
			L = 1;
			if (morpho)	{
				L = 2;
			}
		}

		dim = Ncont + L;
		df = Ncont + L + indf;

		cerr << "tree and data ok\n";

		calibset = 0;
		if (calibfile != "null")	{
			calibset = new FileCalibrationSet(calibfile, tree);
		}
		cerr << "calib ok\n";

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		if (chronoprior == 0)	{// uniform

			double rootage = inrootage;
			double rootstdev = inrootstdev;

			rootalpha = rootage * rootage / rootstdev / rootstdev;
			rootbeta = rootage / rootstdev / rootstdev;

			chronogram = new CalibratedChronogram(tree,One,rootalpha,rootbeta,calibset,false);

		}

		else if (chronoprior == 1)	{ // serial coalescent

			double rootage = inrootage;
			double rootstdev = inrootstdev;
			if (rootstdev)	{
				rootalpha = rootage * rootage / rootstdev / rootstdev;
				rootbeta = rootage / rootstdev / rootstdev;

				CoalAlpha = new Const<PosReal>(rootalpha);
				CoalBeta = new Const<PosReal>(rootbeta);
				coalrate = new Gamma(CoalAlpha,CoalBeta);
			}
			else	{
				CoalAlpha = new Const<PosReal>(1.0);
				CoalBeta = new Const<PosReal>(1.0);
				coalrate = new Gamma(CoalAlpha,CoalBeta);
				coalrate->ClampAt(rootage);
			}

			chronogram = new SerialCoalescent(tree,One,coalrate,calibset,divcutoff);

		}

		else if (chronoprior == 2)	{ // logistic serial coalescent

			T1 = new Const<PosReal>(T1val);

			double rootage = inrootage;
			double rootstdev = inrootstdev;

			rootalpha = rootage * rootage / rootstdev / rootstdev;
			rootbeta = rootage / rootstdev / rootstdev;

			T0Alpha = new Const<PosReal>(rootalpha);
			T0Beta = new Const<PosReal>(rootbeta);
			T0 = new Gamma(T0Alpha,T0Beta);

			if (extratestdev)	{
				double extratealpha = extrateval * extrateval / extratestdev / extratestdev;
				double extratebeta = extrateval / extratestdev / extratestdev;

				extrateAlpha = new Const<PosReal>(extratealpha);
				extrateBeta = new Const<PosReal>(extratebeta);
				extrate = new Gamma(extrateAlpha,extrateBeta);
			}
			else	{
				extrateAlpha = new Const<PosReal>(1.0);
				extrateBeta = new Const<PosReal>(1.0);
				extrate = new Gamma(extrateAlpha,extrateBeta);
				extrate->ClampAt(extrateval);
			}

			if (divratestdev)	{
				double divratealpha = divrateval * divrateval / divratestdev / divratestdev;
				double divratebeta = divrateval / divratestdev / divratestdev;

				divrateAlpha = new Const<PosReal>(divratealpha);
				divrateBeta = new Const<PosReal>(divratebeta);
				divrate = new Gamma(divrateAlpha,divrateBeta);
			}
			else	{
				divrateAlpha = new Const<PosReal>(1.0);
				divrateBeta = new Const<PosReal>(1.0);
				divrate = new Gamma(divrateAlpha,divrateBeta);
				divrate->ClampAt(divrateval);
			}

			if (massextstdev)	{
				double p = massextval;
				double m = massextstdev;
				/*
				double m = p * (1-p) / massextstdev / massextstdev - 1;
				if (m<0)	{
					cerr << "error: beta parameter of mass ext is negative\n";
					exit(1);
				}
				*/
				double massextalpha = p * m;
				double massextbeta = (1-p) * m;

				massextAlpha = new Const<PosReal>(massextalpha);
				massextBeta = new Const<PosReal>(massextbeta);
				massext = new Beta(massextAlpha,massextBeta);
			}
			else	{
				massextAlpha = new Const<PosReal>(1.0);
				massextBeta = new Const<PosReal>(1.0);
				massext = new Beta(massextAlpha,massextBeta);
				massext->ClampAt(massextval);
			}

			if (K0stdev)	{
				double K0alpha = K0val * K0val / K0stdev / K0stdev;
				double K0beta = K0val / K0stdev / K0stdev;

				K0Alpha = new Const<PosReal>(K0alpha);
				K0Beta = new Const<PosReal>(K0beta);
				K0 = new Gamma(K0Alpha,K0Beta);
			}
			else	{
				K0Alpha = new Const<PosReal>(1.0);
				K0Beta = new Const<PosReal>(1.0);
				K0 = new Gamma(K0Alpha,K0Beta);
				K0->ClampAt(K0val);
			}

			if (K1stdev)	{
				double K1alpha = K1val * K1val / K1stdev / K1stdev;
				double K1beta = K1val / K1stdev / K1stdev;

				K1Alpha = new Const<PosReal>(K1alpha);
				K1Beta = new Const<PosReal>(K1beta);
				K1 = new Gamma(K1Alpha,K1Beta);
			}
			else	{
				K1Alpha = new Const<PosReal>(1.0);
				K1Beta = new Const<PosReal>(1.0);
				K1 = new Gamma(K1Alpha,K1Beta);
				K1->ClampAt(K1val);
			}

			chronogram = new LogisticSerialCoalescent(tree,One,divrate,extrate,K0,massext,K1,T0,T1,calibset,divcutoff);

		}

		else if (chronoprior == 3)	{ // serial BD

			double rootage = inrootage;
			double rootstdev = inrootstdev;

			rootalpha = rootage * rootage / rootstdev / rootstdev;
			rootbeta = rootage / rootstdev / rootstdev;

			if (bddivratestdev)	{
				double bddivratealpha = bddivratemean * bddivratemean / bddivratestdev / bddivratestdev;
				double bddivratebeta = bddivratemean / bddivratestdev / bddivratestdev;

				bddivrateAlpha = new Const<PosReal>(bddivratealpha);
				bddivrateBeta = new Const<PosReal>(bddivratebeta);
				bddivrate = new Gamma(bddivrateAlpha,bddivrateBeta);
				bddivrate->setval(0.1);
			}
			else	{
				bddivrateAlpha = new Const<PosReal>(1.0);
				bddivrateBeta = new Const<PosReal>(1.0);
				bddivrate = new Gamma(bddivrateAlpha,bddivrateBeta);
				bddivrate->ClampAt(bddivratemean);
			}

			if (mustdev)	{
				double mualpha = mumean * mumean / mustdev / mustdev;
				double mubeta = mumean / mustdev / mustdev;

				muAlpha = new Const<PosReal>(mualpha);
				muBeta = new Const<PosReal>(mubeta);
				mu = new Gamma(muAlpha,muBeta);
				mu->setval(1.0);
			}
			else	{
				muAlpha = new Const<PosReal>(1.0);
				muBeta = new Const<PosReal>(1.0);
				mu = new Gamma(muAlpha,muBeta);
				mu->ClampAt(mumean);
			}

			if (psistdev)	{
				double psialpha = psimean * psimean / psistdev / psistdev;
				double psibeta = psimean / psistdev / psistdev;

				psiAlpha = new Const<PosReal>(psialpha);
				psiBeta = new Const<PosReal>(psibeta);
				psi = new Gamma(psiAlpha,psiBeta);
				psi->setval(0.1);
			}
			else	{
				psiAlpha = new Const<PosReal>(1.0);
				psiBeta = new Const<PosReal>(1.0);
				psi = new Gamma(psiAlpha,psiBeta);
				psi->ClampAt(psimean);
			}

			if (rhostdev)	{
				double p = rhomean;
				double m = p * (1-p) / rhostdev - 1;
				if (m<0)	{
					cerr << "error: beta parameter of mass ext is negative\n";
					exit(1);
				}
				double rhoalpha = p * m;
				double rhobeta = (1-p) * m;

				rhoAlpha = new Const<PosReal>(rhoalpha);
				rhoBeta = new Const<PosReal>(rhobeta);
				rho = new Beta(rhoAlpha,rhoBeta);
			}
			else	{
				rhoAlpha = new Const<PosReal>(1.0);
				rhoBeta = new Const<PosReal>(1.0);
				rho = new Beta(rhoAlpha,rhoBeta);
				rho->ClampAt(rhomean);
			}

			chronogram = new SerialBirthDeath(tree,One,rootalpha,rootbeta,bddivrate,mu,psi,rho,calibset,divcutoff,Nextant);
		}

		else if (chronoprior == 4)	{ // birth death Rannala and Yang's style

			double rootage = inrootage;
			double rootstdev = inrootstdev;

			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;

			if (bddivratestdev)	{
				double bddivratealpha = bddivratemean * bddivratemean / bddivratestdev / bddivratestdev;
				double bddivratebeta = bddivratemean / bddivratestdev / bddivratestdev;

				bddivrateAlpha = new Const<PosReal>(bddivratealpha);
				bddivrateBeta = new Const<PosReal>(bddivratebeta);
				Chi = new Gamma(bddivrateAlpha,bddivrateBeta);
			}
			else	{
				bddivrateAlpha = new Const<PosReal>(1.0);
				bddivrateBeta = new Const<PosReal>(1.0);
				Chi = new Gamma(bddivrateAlpha,bddivrateBeta);
				Chi->ClampAt(bddivratemean);
			}

			if (mustdev)	{
				double mualpha = mumean * mumean / mustdev / mustdev;
				double mubeta = mumean / mustdev / mustdev;

				muAlpha = new Const<PosReal>(mualpha);
				muBeta = new Const<PosReal>(mubeta);
				Chi2 = new Gamma(muAlpha,muBeta);
			}
			else	{
				muAlpha = new Const<PosReal>(1.0);
				muBeta = new Const<PosReal>(1.0);
				Chi2 = new Gamma(muAlpha,muBeta);
				Chi2->ClampAt(mumean);
			}

			chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,bdpriortype);

		}
		else 	{
			cerr << "error in prior on divergence times\n";
			exit(1);
		}

		if (clamptree)	{
			chronogram->Clamp();
		}

		if (! prior)	{
		if (scaling && brownian)	{
			cerr << "scaling + brownian not yet implemented\n";
			exit(1);
		}

		if (scaling == 1)	{ // burst

			scalet0 = 0;

			if (scalet0val)	{
				scalet0 = new Const<PosReal>(scalet0val);
			}
			else	{
				scalet0 = chronogram->GetScale();
			}

			if (scalefactorstdev)	{
				double scalefactoralpha = scalefactorval * scalefactorval / scalefactorstdev / scalefactorstdev;
				double scalefactorbeta = scalefactorval / scalefactorstdev / scalefactorstdev;

				scalefactorAlpha = new Const<PosReal>(scalefactoralpha);
				scalefactorBeta = new Const<PosReal>(scalefactorbeta);
				scalefactor = new Gamma(scalefactorAlpha,scalefactorBeta);
			}
			else	{
				scalefactorAlpha = new Const<PosReal>(1.0);
				scalefactorBeta = new Const<PosReal>(1.0);
				scalefactor = new Gamma(scalefactorAlpha,scalefactorBeta);
				scalefactor->ClampAt(scalefactorval);
			}

			if (scaleratestdev)	{
				double scaleratealpha = scalerateval * scalerateval / scaleratestdev / scaleratestdev;
				double scaleratebeta = scalerateval / scaleratestdev / scaleratestdev;

				scalerateAlpha = new Const<PosReal>(scaleratealpha);
				scalerateBeta = new Const<PosReal>(scaleratebeta);
				scalerate = new Gamma(scalerateAlpha,scalerateBeta);
			}
			else	{
				scalerateAlpha = new Const<PosReal>(1.0);
				scalerateBeta = new Const<PosReal>(1.0);
				scalerate = new Gamma(scalerateAlpha,scalerateBeta);
				scalerate->ClampAt(scalerateval);
			}

			scalefunction = new BurstScalingFunction(scalet0,scalefactor,scalerate);
		}
		else if (scaling == 2)	{
			scalesigma = new Gamma(One,One);
			gammascaletree = new GammaWhiteNoiseProcess(chronogram,scalesigma,0,MEAN,true);
		}

		if (dim)	{

			double mindiag = 0.001;
			double maxdiag = 1000;
			DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
			DiagArray->setval(1.0);
			// create a diagonal matrix, with the kappa_i along the diagonal
			sigmaZero = new SigmaZero(DiagArray);

			// create covariance matrix
			// from an inverse wishart of parameter sigma0
			cerr << "sigma\n";
			if (clampdiag)	{
				sigma = new DiagonalCovMatrix(sigmaZero, df);
			}
			else	{
				conjugatesigma = new ConjugateInverseWishart(sigmaZero, df);
				sigma = conjugatesigma;
			}

		}

		cerr << "process\n";

		if (brownian)	{

			if (scaling)	{
				cerr << "error: scaling with brownian process not yet implemented\n";
				exit(1);
			}

			cerr << "fine-grained brownian process\n";
			cerr << "nseg : " << nSegments << '\n';
			BrownianBridge::setNTreeSegments(nSegments);
			BrownianBridge::setSegmentation(segm);

			if (! clampdiag)	{
				conjugateBrownianProcess = new ConjugateBrownianProcess(chronogram,conjugatesigma);
				// conjugateBrownianProcess = new ConjugateBrownianProcess(chronogram,conjugatesigma);
				brownianProcess = conjugateBrownianProcess;
			}
			else	{
				brownianProcess = new BrownianProcess(chronogram,sigma);
				// brownianProcess = new BrownianProcess(chronogram,sigma);
			}

			process = brownianProcess->GetInstantProcess();

			for (int i=0; i<Ncont; i++)	{
				brownianProcess->GetInstantProcess()->SetAndClamp(contdata,L+i,i);
			}

			// just for numerical stability of the starting point
			for (int l=0; l<L; l++)	{
					brownianProcess->GetInstantProcess()->CutOff(1,l);
			}

			nuccommutativeProcess = new CommutativeBrownianProcess(brownianProcess,0);
			if (morpho)	{
				morphocommutativeProcess = new CommutativeBrownianProcess(brownianProcess,1);
			}

			if (clockprior == 0)	{

				nucratetree = nuccommutativeProcess;
				morphoratetree = morphocommutativeProcess;
			}
			else if (clockprior == 1)	{

				cerr << "error: brownian process only with Brownian clock\n";
				exit(1);
			}
			else if (clockprior == 2)	{

				wnnucsigma = new Gamma(One,One);
				wnnucratetree = new GammaWhiteNoiseProcess(chronogram,wnnucsigma,0,MEAN,true);
				productnucratetree = new BranchProductProcess(nuccommutativeProcess,wnnucratetree);
				nucratetree = productnucratetree;

				if (morpho)	{
					wnmorphosigma = new Gamma(One,One);
					wnmorphoratetree = new GammaWhiteNoiseProcess(chronogram,wnmorphosigma,0,MEAN,true);
					productmorphoratetree = new BranchProductProcess(morphocommutativeProcess,wnmorphoratetree);
					morphoratetree = productmorphoratetree;
				}

			}
			else	{
				cerr << "error in clock prior\n";
				exit(1);
			}
		}

		else	{

			if (dim)	{
				if (clampdiag)	{
					if (scaling == 1)	{
						process = new MultiVariateTreeProcess(sigma,chronogram,chronogram->GetScale(),scalefunction);
					}
					else if (scaling == 2)	{
						process = new MultiVariateTreeProcess(sigma,chronogram,gammascaletree);
					}
					else	{
						process = new MultiVariateTreeProcess(sigma,chronogram);
					}
				}
				else	{
					if (scaling == 1)	{
						process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram,chronogram->GetScale(),scalefunction);
					}
					else if (scaling == 2)	{
						process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram,gammascaletree);
					}
					else	{
						process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram);
					}
				}

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
			}

			// create the branch lengths resulting from combining
			// the times given by the chronogram with the rate (first entry of the multivariate process)

			cerr << "branch rates\n";
			if (clockprior == 0)	{

				meanexpnucratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
				nucratetree = meanexpnucratetree;

				if (morpho)	{
					meanexpmorphoratetree = new MeanExpTreeFromMultiVariate(process,1,INTEGRAL,false,meanexp);
					morphoratetree = meanexpmorphoratetree;
				}
			}
			else if (clockprior == 1)	{

				wnnucsigma = new Gamma(One,One);
				double ratealpha = ratemean * ratemean / ratestdev / ratestdev;
				double ratebeta = ratemean / ratestdev / ratestdev;
				rateAlpha = new Const<PosReal>(ratealpha);
				rateBeta = new Const<PosReal>(ratebeta);
				wnnucmeanrate = new Gamma(rateAlpha,rateBeta);
				wnnucratetree = new GammaWhiteNoiseProcess(chronogram,wnnucsigma,wnnucmeanrate,INTEGRAL);
				nucratetree = wnnucratetree;

				if (morpho)	{
					wnmorphosigma = new Gamma(One,One);
					double morphoratealpha = morphoratemean * morphoratemean / morphoratestdev / morphoratestdev;
					double morphoratebeta = morphoratemean / morphoratestdev / morphoratestdev;
					morphorateAlpha = new Const<PosReal>(morphoratealpha);
					morphorateBeta = new Const<PosReal>(morphoratebeta);
					wnmorphomeanrate = new Gamma(morphorateAlpha,morphorateBeta);
					wnmorphoratetree = new GammaWhiteNoiseProcess(chronogram,wnmorphosigma,wnmorphomeanrate,INTEGRAL);
					morphoratetree = wnmorphoratetree;
				}
			}
			else if (clockprior == 2)	{

				wnnucsigma = new Gamma(One,One);
				meanexpnucratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
				wnnucratetree = new GammaWhiteNoiseProcess(chronogram,wnnucsigma,0,MEAN,true);
				productnucratetree = new BranchProductProcess(meanexpnucratetree,wnnucratetree);
				nucratetree = productnucratetree;

				if (morpho)	{
					wnmorphosigma = new Gamma(One,One);
					meanexpmorphoratetree = new MeanExpTreeFromMultiVariate(process,1,INTEGRAL,false,meanexp);
					wnmorphoratetree = new GammaWhiteNoiseProcess(chronogram,wnmorphosigma,0,MEAN,true);
					productmorphoratetree = new BranchProductProcess(meanexpmorphoratetree,wnmorphoratetree);
					morphoratetree = productmorphoratetree;
				}

			}
			else	{
				cerr << "error in clock prior\n";
				exit(1);
			}
		}

		// substitution matrix
		int Nstate = nucdata->GetNstate();
		nucrelrate = new Dirichlet(Nstate*(Nstate-1)/2);
		nucstationary = new Dirichlet(Nstate);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(nucrelrate,nucstationary,false);

		if (morpho)	{
			morpho2stationary = new Dirichlet(2);
			morpho3stationary = new Dirichlet(3);
			morpho4stationary = new Dirichlet(4);
			morpho2stationary->setuniform();
			morpho3stationary->setuniform();
			morpho4stationary->setuniform();
			morpho2matrixtree = new F81MatrixTree(morphoratetree,morpho2stationary);
			morpho3matrixtree = new F81MatrixTree(morphoratetree,morpho3stationary);
			morpho4matrixtree = new F81MatrixTree(morphoratetree,morpho4stationary);
		}

		// a phylogenetic process
		cerr << "phyloprocess\n";
		if (conjpath)	{
			// a phylogenetic process
			nucpathconjtree = new OneMatrixPathConjugateTree(nucratetree,nucmatrix,nucdata);
			nucphyloprocess = new PathConjugatePhyloProcess(nucpathconjtree);
			if (morpho)	{
				morpho2pathconjtree = new BranchMatrixTransitionPathConjugateTree(morphoratetree,morpho2matrixtree,morpho2data);
				morpho3pathconjtree = new BranchMatrixTransitionPathConjugateTree(morphoratetree,morpho3matrixtree,morpho3data);
				morpho4pathconjtree = new BranchMatrixTransitionPathConjugateTree(morphoratetree,morpho4matrixtree,morpho4data);
				morpho2phyloprocess = new TransitionPathConjugatePhyloProcess(morpho2pathconjtree);
				morpho3phyloprocess = new TransitionPathConjugatePhyloProcess(morpho3pathconjtree);
				morpho4phyloprocess = new TransitionPathConjugatePhyloProcess(morpho4pathconjtree);
			}
		}
		else	{
			nucphyloprocess = new OneMatrixPhyloProcess(nucratetree,nucmatrix,nucdata);
			if (morpho)	{
				morpho2phyloprocess = new BranchTransitionMatrixPhyloProcess(morphoratetree,morpho2matrixtree,morpho2data);
				morpho3phyloprocess = new BranchTransitionMatrixPhyloProcess(morphoratetree,morpho3matrixtree,morpho3data);
				morpho4phyloprocess = new BranchTransitionMatrixPhyloProcess(morphoratetree,morpho4matrixtree,morpho4data);
			}
		}

		cerr << "unfold\n";
		nucphyloprocess->Unfold();
		if (morpho)	{
			morpho2phyloprocess->Unfold();
			morpho3phyloprocess->Unfold();
			morpho4phyloprocess->Unfold();
		}
		if (sample)	{
			nucphyloprocess->Sample();
			if (morpho)	{
				morpho2phyloprocess->Sample();
				morpho3phyloprocess->Sample();
				morpho4phyloprocess->Sample();
			}
		}

		}

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		if (chronoprior == 1)	{
			RootRegister(CoalAlpha);
			RootRegister(CoalBeta);
		}
		else if (chronoprior == 2)	{
			RootRegister(divrateAlpha);
			RootRegister(extrateAlpha);
			RootRegister(massextAlpha);
			RootRegister(K0Alpha);
			RootRegister(K1Alpha);
			RootRegister(T0Alpha);
			RootRegister(divrateBeta);
			RootRegister(extrateBeta);
			RootRegister(massextBeta);
			RootRegister(K0Beta);
			RootRegister(K1Beta);
			RootRegister(T0Beta);

			RootRegister(T1);
		}
		else if (chronoprior == 3)	{
			RootRegister(bddivrateAlpha);
			RootRegister(bddivrateBeta);
			RootRegister(muAlpha);
			RootRegister(muBeta);
			RootRegister(psiAlpha);
			RootRegister(psiBeta);
			RootRegister(rhoAlpha);
			RootRegister(rhoBeta);
		}
		else if (chronoprior == 4)	{
			RootRegister(bddivrateAlpha);
			RootRegister(bddivrateBeta);
			RootRegister(muAlpha);
			RootRegister(muBeta);
		}
		if (! prior)	{
		if (scaling == 1)	{
			RootRegister(scalet0);
			RootRegister(scalefactorAlpha);
			RootRegister(scalefactorBeta);
			RootRegister(scalerateAlpha);
			RootRegister(scalerateBeta);
		}
		RootRegister(nucrelrate);
		RootRegister(nucstationary);
		if (morpho)	{
			RootRegister(morpho2stationary);
			RootRegister(morpho3stationary);
			RootRegister(morpho4stationary);
		}
		if (clockprior == 1)	{
			RootRegister(rateAlpha);
			RootRegister(rateBeta);
			if (morpho)	{
				RootRegister(morphorateAlpha);
				RootRegister(morphorateBeta);
			}
		}
		}

		Register();

		cerr << "scheduler\n";

		MakeScheduler();

		cerr << "sample\n";
		if (sample)	{
			Update();
			TraceHeader(cerr);
			Trace(cerr);
		}

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

	void FastUpdate()	{

		chronogram->specialUpdate();
		
		if (! prior)	{
		if (brownian)	{

			nuccommutativeProcess->specialUpdate();
			if (morpho)	{
				morphocommutativeProcess->specialUpdate();
			}

			if (clockprior == 2)	{
				productnucratetree->specialUpdate();
				if (morpho)	{
					productmorphoratetree->specialUpdate();
				}
			}
		}
		else	{

			if (clockprior != 1)	{
				meanexpnucratetree->specialUpdate();
				if (morpho)	{
					meanexpmorphoratetree->specialUpdate();
				}
			}
			if (clockprior == 2)	{
				productnucratetree->specialUpdate();
				if (morpho)	{
					productmorphoratetree->specialUpdate();
				}
			}
		}
		}
	}

	// MeanExpTreeFromMultiVariate* GetNucRateTree() {return nucratetree;}
	// MeanExpTreeFromMultiVariate* GetMorphoRateTree() {return morphoratetree;}
	LengthTree* GetNucRateTree() {return nucratetree;}
	LengthTree* GetMorphoRateTree() {return morphoratetree;}

	CalibratedChronogram* GetChronogram() {return chronogram;}

	CovMatrix* GetCovMatrix() {

		if (! dim)	{
			cerr << "error: no cov matrix\n";
			exit(1);
		}
		return sigma;
	}

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

	MultiVariateTreeProcess* GetMultiVariateProcess() {

		if (! dim)	{
			cerr << "error: no multivariate process\n";
			exit(1);
		}
		if (brownian)	{
			return brownianProcess->GetInstantProcess();
		}
		return process;
	}

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	int GetDim() {return dim;}

	int GetNcont() {return Ncont;}

	int Morpho() {return morpho;}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;
		if (chronoprior == 1)	{
			total += coalrate->GetLogProb();
		}
		else if (chronoprior == 2)	{
			total += T0->GetLogProb();
			total += K0->GetLogProb();
			total += K1->GetLogProb();
			total += extrate->GetLogProb();
			total += divrate->GetLogProb();
			total += massext->GetLogProb();
		}
		else if (chronoprior == 3)	{
			total += bddivrate->GetLogProb();
			total += mu->GetLogProb();
			total += psi->GetLogProb();
			if ((rho->val() > 0) && (rho->val() < 1))	{
				total += rho->GetLogProb();
			}
		}
		else if (chronoprior == 4)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();
		if (! prior)	{
		if (scaling == 1)	{
			total += scalefactor->GetLogProb();
			total += scalerate->GetLogProb();
		}
		else if (scaling == 2)	{
			total += scalesigma->GetLogProb();
			total += gammascaletree->GetLogProb();
		}
		if (dim)	{
			total += DiagArray->GetLogProb();
			total += sigma->GetLogProb();
			if (brownian)	{
				total += brownianProcess->GetLogProb();
			}
			else	{
				total += process->GetLogProb();
			}
		}
		if (clockprior == 1)	{
			total += wnnucsigma->GetLogProb();
			total += wnnucmeanrate->GetLogProb();
			total += wnnucratetree->GetLogProb();
			if (morpho)	{
				total += wnmorphosigma->GetLogProb();
				total += wnmorphomeanrate->GetLogProb();
				total += wnmorphoratetree->GetLogProb();
			}
		}
		else if (clockprior == 2)	{
			total += wnnucsigma->GetLogProb();
			total += wnnucratetree->GetLogProb();
			if (morpho)	{
				total += wnmorphosigma->GetLogProb();
				total += wnmorphoratetree->GetLogProb();
			}
		}
		total += nucrelrate->GetLogProb();
		total += nucstationary->GetLogProb();
		}
		return total;
	}

	double GetLogLikelihood()	{
		if (prior)	{
			return 0;
		}
		double ret = 0;
		ret += nucphyloprocess->GetLogProb();
		if (morpho)	{
			ret += morpho2phyloprocess->GetLogProb();
			ret += morpho3phyloprocess->GetLogProb();
			ret += morpho4phyloprocess->GetLogProb();
		}
		return ret;
	}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	Tree* GetTree()	{
		return tree;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		if (! prior)	{
		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(nucphyloprocess,nucpathconjtree),1,"mapping + sufficient stat");
			if (morpho)	{
				scheduler.Register(new DSemiConjugateMappingMove(morpho2phyloprocess,morpho2pathconjtree),1,"mapping + sufficient stat");
				scheduler.Register(new DSemiConjugateMappingMove(morpho3phyloprocess,morpho3pathconjtree),1,"mapping + sufficient stat");
				scheduler.Register(new DSemiConjugateMappingMove(morpho4phyloprocess,morpho4pathconjtree),1,"mapping + sufficient stat");
			}
		}
		else	{
			scheduler.Register(new SimpleMove(nucphyloprocess,1),1,"mapping");
			if (morpho)	{
				scheduler.Register(new SimpleMove(morpho2phyloprocess,1),1,"mapping");
				scheduler.Register(new SimpleMove(morpho3phyloprocess,1),1,"mapping");
				scheduler.Register(new SimpleMove(morpho4phyloprocess,1),1,"mapping");
			}
		}
		}

		int nrep = conjpath ? 30 : 1;
		for (int i=0; i<nrep; i++)	{

			if (! prior)	{
			if (brownian)	{
				if (clampdiag)	{

					scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.1),10,"PureBrownian");
					scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.01),10,"PureBrownian");
					scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),1),10,"InstantValues");
					scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),0.1),10,"InstantValues");
				 
					scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
					scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
					scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
					scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
				}
				else	{
 					if (conjsigma == 2)	{
						scheduler.Register(new ConjugateBrownianExternalMove(conjugatesigma, conjugateBrownianProcess, 0.1,10) ,1,"External Conjugate Brownian sigma");
						scheduler.Register(new ConjugateBrownianExternalMove(conjugatesigma, conjugateBrownianProcess, 0.01,10) ,1,"External Conjugate Brownian sigma");
						scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.1,10) ,1,"External Conjugate Brownian all branch sigma");
						scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.01,10) ,1,"External Conjugate Brownian all branch sigma");
						scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.001,10) ,1,"External Conjugate Brownian all branch sigma");
					}
					else	{

						scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.1),10,"PureBrownian");
						scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.01),10,"PureBrownian");
						scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),1),10,"InstantValues");
						scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),0.1),10,"InstantValues");
					 
						scheduler.Register(new BrownianSigmaMove(conjugatesigma, brownianProcess, 0.1), 10, "linear sigma brownian");
						scheduler.Register(new BrownianSigmaMove(conjugatesigma, brownianProcess, 0.01), 10, "linear sigma brownian");
						scheduler.Register(new BrownianSigmaMove(conjugatesigma, brownianProcess, 0.001), 10, "linear sigma brownian");

						scheduler.Register(new ConjugateBrownianSigmaMove(conjugatesigma, brownianProcess), 1, "sigma(Gibbs)");
					}
				}
			}
			else	{

				if (dim)	{
					scheduler.Register(new SimpleMove(process,1),10,"multinormal");
					scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
					scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");

					if (clampdiag)	{
						scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
						scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
						scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
						scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");

					}
					else	{
						scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
						scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
						scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
						scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");
					}
				}
			}

			if (dim)	{
				scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
				scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
				scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");
			}

			}

			if (chronoprior == 1)	{
				scheduler.Register(new SimpleMove(coalrate,1),10,"coal rate");
				scheduler.Register(new SimpleMove(coalrate,0.1),10,"coal rate");
			}
			else if (chronoprior == 2)	{
				scheduler.Register(new SimpleMove(T0,1),10,"T0");
				scheduler.Register(new SimpleMove(T0,0.1),10,"T0");
				scheduler.Register(new SimpleMove(K0,1),10,"K0");
				scheduler.Register(new SimpleMove(K0,0.1),10,"K0");
				scheduler.Register(new SimpleMove(K1,1),10,"K1");
				scheduler.Register(new SimpleMove(K1,0.1),10,"K1");
				scheduler.Register(new SimpleMove(divrate,1),10,"divrate");
				scheduler.Register(new SimpleMove(divrate,0.1),10,"divrate");
				scheduler.Register(new SimpleMove(extrate,1),10,"extrate");
				scheduler.Register(new SimpleMove(extrate,0.1),10,"extrate");
				scheduler.Register(new SimpleMove(massext,1),10,"massext");
				scheduler.Register(new SimpleMove(massext,0.1),10,"massext");
			}
			else if (chronoprior == 3)	{
				scheduler.Register(new SimpleMove(bddivrate,1),10,"bddivrate");
				scheduler.Register(new SimpleMove(bddivrate,0.1),10,"bddivrate");
				scheduler.Register(new SimpleMove(mu,1),10,"mu");
				scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
				scheduler.Register(new SimpleMove(psi,1),10,"psi");
				scheduler.Register(new SimpleMove(psi,0.1),10,"psi");
				scheduler.Register(new SimpleMove(rho,1),10,"rho");
				scheduler.Register(new SimpleMove(rho,0.1),10,"rho");
			}
			else if (chronoprior == 4)	{
				scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
			}

			if (!prior)	{

			if (clockprior)	{

				scheduler.Register(new SimpleMove(wnnucratetree,1),10,"white noise");
				scheduler.Register(new SimpleMove(wnnucratetree,0.1),10,"white noise");
				scheduler.Register(new SimpleMove(wnnucratetree,0.01),10,"white noise");

				scheduler.Register(new SimpleMove(wnnucsigma,1),10,"wn sigma");
				scheduler.Register(new SimpleMove(wnnucsigma,0.1),10,"wn sigma");
				scheduler.Register(new SimpleMove(wnnucsigma,0.01),10,"wn sigma");

				if (morpho)	{
					scheduler.Register(new SimpleMove(wnmorphoratetree,1),10,"morpho white noise");
					scheduler.Register(new SimpleMove(wnmorphoratetree,0.1),10,"morpho white noise");
					scheduler.Register(new SimpleMove(wnmorphoratetree,0.01),10,"morpho white noise");

					scheduler.Register(new SimpleMove(wnmorphosigma,1),10,"wn sigma");
					scheduler.Register(new SimpleMove(wnmorphosigma,0.1),10,"wn sigma");
					scheduler.Register(new SimpleMove(wnmorphosigma,0.01),10,"wn sigma");
				}
			}

			if (clockprior == 2)	{

				if (brownian)	{

					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,brownianProcess->GetInstantProcess(),0,1),10,"wnprocess comp");
					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,brownianProcess->GetInstantProcess(),0,0.1),10,"wnprocess comp");
					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,brownianProcess->GetInstantProcess(),0,0.01),10,"wnprocess comp");
					if (morpho)	{
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,brownianProcess->GetInstantProcess(),1,1),10,"wnprocess comp");
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,brownianProcess->GetInstantProcess(),1,0.1),10,"wnprocess comp");
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,brownianProcess->GetInstantProcess(),1,0.01),10,"wnprocess comp");
					}

				}
				else	{

					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,process,0,1),10,"wnprocess comp");
					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,process,0,0.1),10,"wnprocess comp");
					scheduler.Register(new MultiVariateUgamCompMove(wnnucratetree,process,0,0.01),10,"wnprocess comp");

					if (morpho)	{
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,process,1,1),10,"wnprocess comp");
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,process,1,0.1),10,"wnprocess comp");
						scheduler.Register(new MultiVariateUgamCompMove(wnmorphoratetree,process,1,0.01),10,"wnprocess comp");
					}
				}
			}

			if (clockprior == 1)	{

				scheduler.Register(new SimpleMove(wnnucmeanrate,1),10,"mean rate");
				scheduler.Register(new SimpleMove(wnnucmeanrate,0.1),10,"mean rate");
				scheduler.Register(new SimpleMove(wnnucmeanrate,0.01),10,"mean rate");

				if (morpho)	{
					scheduler.Register(new SimpleMove(wnmorphomeanrate,1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnmorphomeanrate,0.1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnmorphomeanrate,0.01),10,"mean rate");
				}

			}
			
			if (scaling == 1)	{

				scheduler.Register(new SimpleMove(scalefactor,10),10,"scalefactor");
				scheduler.Register(new SimpleMove(scalefactor,1),10,"scalefactor");
				scheduler.Register(new SimpleMove(scalefactor,0.1),10,"scalefactor");

				scheduler.Register(new SimpleMove(scalerate,10),10,"scalerate");
				scheduler.Register(new SimpleMove(scalerate,1),10,"scalerate");
				scheduler.Register(new SimpleMove(scalerate,0.1),10,"scalerate");
			}
			else if (scaling == 2)	{

				scheduler.Register(new SimpleMove(gammascaletree,1),10,"scale gamma tree");
				scheduler.Register(new SimpleMove(gammascaletree,0.1),10,"scale gamma tree");
				scheduler.Register(new SimpleMove(gammascaletree,0.01),10,"scale gamma tree");

				scheduler.Register(new SimpleMove(scalesigma,10),10,"scalesigma");
				scheduler.Register(new SimpleMove(scalesigma,1),10,"scalesigma");
				scheduler.Register(new SimpleMove(scalesigma,0.1),10,"scalesigma");

			}

			}

			if (! fixbl)	{
			if(brownian && (segm == SEGM_ABSOLUTE)) {

				cerr << "absolute segmentation not yet implemented\n";
				exit(1);

				/*
				scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 1), 10, "horizontal move");
				scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 0.1), 10, "horizontal move");
				scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 0.01), 10, "horizontal move");

				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				*/
			}
			else	{

				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.3),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.03),10,"chrono");

				/*
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				*/

				/*
				scheduler.Register(new AllInternalNodesMove(chronogram,1),10,"all int chrono");
				scheduler.Register(new AllInternalNodesMove(chronogram,0.1),10,"all int chrono");
				scheduler.Register(new AllInternalNodesMove(chronogram,0.01),10,"all int chrono");
				*/

				if (chronoprior == 4)	{

					scheduler.Register(new SimpleMove(chronogram->GetScale(),1),10,"root age");
					scheduler.Register(new SimpleMove(chronogram->GetScale(),0.1),10,"root age");
					scheduler.Register(new SimpleMove(chronogram->GetScale(),0.01),10,"root age");
					scheduler.Register(new SimpleMove(chronogram->GetScale(),0.001),10,"root age");
				}
				else	{

					scheduler.Register(new ChronoRootOnlyMove(chronogram,1),10,"chrono root only move");
					scheduler.Register(new ChronoRootOnlyMove(chronogram,0.1),10,"chrono root only move");
					scheduler.Register(new ChronoRootOnlyMove(chronogram,0.03),10,"chrono root only move");
				}
			}
			}

			if (! prior)	{
			scheduler.Register(new ProfileMove(nucrelrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(nucrelrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(nucrelrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(nucstationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(nucstationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(nucstationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(nucstationary,0.001),10,"stat");
			}
		}
	}

	void drawSample()	{
		cerr << "sample\n";
		exit(1);
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetRootAge()	{
		return GetCalibratedChronogram()->GetScale()->val();
	}

	double GetCoalRate()	{
		return coalrate->val();
	}

	void TraceChrono(ostream& os)	{
		chronogram->Newick(os);
	}

	double GetT0()	{
		return T0->val() + chronogram->GetScale()->val();
	}

	double GetNucTotalLength()	{
		/*
		if (brownian)	{
			return brownianProcess->GetIntegralRate(0);
		}
		*/
		return GetNucRateTree()->GetTotalLength();
	}

	double GetMorphoTotalLength()	{

		return GetMorphoRateTree()->GetTotalLength();
	}

	double GetTotalMorphoRateVariance()	{

		if (! morpho)	{
			return 0;
		}
		return GetVariance(morphoratetree,true);
	}

	double GetLnMorphoRateVariance()	{

		if ((! morpho) || (! dim))	{
			return 0;
		}
		if (brownian)	{
			return GetVariance(morphocommutativeProcess,true);
		}
		return GetVariance(meanexpmorphoratetree,true);
	}

	double GetWnMorphoRateVariance()	{

		if ((! morpho) || (! clockprior))	{
			return 0;
		}
		bool integral = (clockprior == 1);
		return GetVariance(wnmorphoratetree,integral);
	}

	double GetTotalRateVariance()	{

		return GetVariance(nucratetree,true);
	}

	double GetLnRateVariance()	{

		if (! dim)	{
			return 0;
		}
		if (brownian)	{
			return GetVariance(nuccommutativeProcess,true);
		}
		return GetVariance(meanexpnucratetree,true);
	}

	double GetWnRateVariance()	{

		if (! clockprior)	{
			return 0;
		}
		bool integral = (clockprior == 1);
		return GetVariance(wnnucratetree,integral);
	}

	double GetVariance(LengthTree* lengthtree, bool integral)	{

		double mean = 0;
		double var = 0;
		double weight = 0;
		RecursiveGetMeanAndVar(lengthtree,GetTree()->GetRoot(),mean,var,weight,integral);
		mean /= weight;
		var /= weight;
		var -= mean*mean;
		return var;
	}
		
	void RecursiveGetMeanAndVar(LengthTree* lengthtree, const Link* from, double& mean, double& var, double& weight, bool integral)	{

		if (! from->isRoot())	{

			double rate = 1;
			double time = chronogram->GetBranchVal(from->GetBranch())->val();
			if (time > 1e-8)	{
				if (integral)	{
					double length = lengthtree->GetBranchVal(from->GetBranch())->val();
					rate = length / time;
				}
				else	{
					rate = lengthtree->GetBranchVal(from->GetBranch())->val();
				}
				double lograte = (rate < 1e-8) ? 0 : log(rate);

				weight += time;
				mean += time * lograte;
				var += time * lograte * lograte;
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetMeanAndVar(lengthtree,link->Out(),mean,var,weight,integral);
		}
	}

	void GetInternalTipVar(double& pint, double& ptip, double& meantot, double& meanint, double& meantip, double& vartot, double& varint, double& vartip)	{

		if (!clockprior)	{
			return;
		}

		double mint0 = 0;
		double mint1 = 0;
		double mint2 = 0;
		double mtip0 = 0;
		double mtip1 = 0;
		double mtip2 = 0;

		RecursiveGetInternalTipMoments(GetTree()->GetRoot(),mint0,mint1,mint2,mtip0,mtip1,mtip2);

		double mtot0 = mint0 + mtip0;
		double mtot1 = mint1 + mtip1;
		double mtot2 = mint2 + mtip2;

		meantot = mtot1 / mtot0;
		vartot = mtot2 / mtot0 - meantot*meantot;

		meanint = mint1 / mint0;
		varint = mint2 / mint0 - meanint*meanint;

		meantip = mtip1 / mtip0;
		vartip = mtip2 / mtip0 - meantip*meantip;

		pint = mint0 / mtot0;
		ptip = mtip0 / mtot0;
	}

	void RecursiveGetInternalTipMoments(const Link* from, double& mint0, double& mint1, double& mint2, double& mtip0, double& mtip1, double& mtip2)	{

		if (! from->isRoot())	{

			double length = wnnucratetree->GetBranchVal(from->GetBranch())->val();
			double time = chronogram->GetBranchVal(from->GetBranch())->val();
			if (time < 1e-8)	{
				time = 1e-6;
				/*
				cerr << "error in get rate: " << time << '\n';
				exit(1);
				*/
			}
			double rate = length / time;
			if (from->isLeaf())	{
				mtip0 += time;
				mtip1 += time * rate;
				mtip2 += time * rate * rate;
			}
			else	{
				mint0 += time;
				mint1 += time * rate;
				mint2 += time * rate * rate;
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetInternalTipMoments(link->Out(),mint0,mint1,mint2,mtip0,mtip1,mtip2);
		}
	}

	double GetGlobalMeanInternalTipDiff()	{
		return GetMeanInternalTipDiff(nucratetree);
	}

	double GetWnMeanInternalTipDiff()	{
		if (! clockprior)	{
			return 0;
		}
		return GetMeanInternalTipDiff(wnnucratetree);
	}

	double GetLnMeanInternalTipDiff()	{
		if (! dim)	{
			return 0;
		}
		if (brownian)	{
			return GetMeanInternalTipDiff(nuccommutativeProcess);
		}
		return GetMeanInternalTipDiff(meanexpnucratetree);
	}

	double GetMeanInternalTipDiff(LengthTree* lengthtree)	{
		double meanint = 0;
		double nint = 0;
		double meantip = 0;
		double ntip = 0;
		RecursiveGetMeanInternalTipDiff(lengthtree,GetTree()->GetRoot(),meanint,nint,meantip,ntip);
		double mean = meanint + meantip;
		double n = nint + ntip;
		meanint /= nint;
		meantip /= ntip;
		mean /= n;
		return (meanint - meantip) / mean;
	}

	void RecursiveGetMeanInternalTipDiff(LengthTree* lengthtree, const Link* from, double& meanint, double& nint, double& meantip, double& ntip)	{

		if (! from->isRoot())	{
			double length = lengthtree->GetBranchVal(from->GetBranch())->val();
			// double length = wnnucratetree->GetBranchVal(from->GetBranch())->val();
			double time = chronogram->GetBranchVal(from->GetBranch())->val();
			if (time < 1e-8)	{
				time = 1e-6;
				/*
				cerr << "error in get rate: " << time << '\n';
				exit(1);
				*/
			}
			double rate = length / time;
			if (from->isLeaf())	{
				meantip += time * rate;
				ntip += time;
			}
			else	{
				meanint += time * rate;
				nint += time;
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetMeanInternalTipDiff(lengthtree,link->Out(),meanint,nint,meantip,ntip);
		}
	}

	BranchVarTree<PosReal>* GetWnNucRateTree()	{
		if (! clockprior)	{
			return 0;
		}
		return wnnucratetree;
	}

	double GetWnNucSigma()	{
		return wnnucsigma->val();
	}

	double* GetRelRate()	{
		return nucrelrate->GetArray();
	}

	double* GetStationary()	{
		return nucstationary->GetArray();
	}

	double GetNstate()	{
		return nucdata->GetNstate();
	}

	double GetLength(const Link* link)	{
		return nucratetree->GetBranchVal(link->GetBranch())->val();
	}

	double GetWnRate(const Link* link)	{
		double length = wnnucratetree->GetBranchVal(link->GetBranch())->val();
		double time = chronogram->GetBranchVal(link->GetBranch())->val();
		if (time < 1e-8)	{
			time = 1e-6;
			/*
			cerr << "error in get rate: " << time << '\n';
			exit(1);
			*/
		}
		return length / time;
	}

	double GetRate(const Link* link)	{
		double length = nucratetree->GetBranchVal(link->GetBranch())->val();
		double time = chronogram->GetBranchVal(link->GetBranch())->val();
		if (time < 1e-8)	{
			time = 1e-6;
			/*
			cerr << "error in get rate: " << time << '\n';
			exit(1);
			*/
		}
		return length / time;
	}

	double GetMaxRate()	{
		return RecursiveGetMaxRate(GetTree()->GetRoot());
	}

	double RecursiveGetMaxRate(const Link* from)	{

		double max = 0;
		if (! from->isRoot())	{
			max = GetRate(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMaxRate(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		return max;
	}

	double GetMinRate()	{
		return RecursiveGetMinRate(GetTree()->GetRoot());
	}

	double RecursiveGetMinRate(const Link* from)	{

		double min = -1;
		if (! from->isRoot())	{
			min = GetRate(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMinRate(link->Out());
			if ((min == -1) || (min > tmp))	{
				min = tmp;
			}
		}
		return min;
	}

	double GetMaxLength()	{
		return RecursiveGetMaxLength(GetTree()->GetRoot());
	}

	double RecursiveGetMaxLength(const Link* from)	{

		double max = 0;
		if (! from->isRoot())	{
			max = GetLength(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMaxLength(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		return max;
	}

	double GetMinLength()	{
		return RecursiveGetMinLength(GetTree()->GetRoot());
	}

	double RecursiveGetMinLength(const Link* from)	{

		double min = -1;
		if (! from->isRoot())	{
			min = GetLength(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMinLength(link->Out());
			if ((min == -1) || (min > tmp))	{
				min = tmp;
			}
		}
		return min;
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

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\trootage";
		if (chronoprior == 1)	{
			os << "\tcoalrate";
		}
		else if (chronoprior == 2)	{
			os << "\tDT0\tK0\tK1\tmassext\tdivrate\textrate";
		}
		else if (chronoprior == 3)	{
			os << "\tbddivrate\tmu\tpsi\trho";
		}
		else if (chronoprior == 4)	{
			os << "\tchi\tchi2";
		}
		if (! prior)	{
		if (clockprior == 1)	{
			os << "\twnmean\twnvar";
			if (morpho)	{
				os << "\twnmorphomean\twnmorphovvar";
			}
		}
		else if (clockprior == 2)	{
			os << "\twnvar";
			if (morpho)	{
				os << "\twnmorphovar";
			}
		}
		if (scaling == 1)	{
			os << "\tscalefactor\tscalerate";
		}
		else if (scaling == 2)	{
			os << "\tscalesigma";
		}
		os << "\tnuclength";
		if (morpho)	{
			os << "\tmorpholength";
		}
		if (dim)	{
		if (!clampdiag)	{
			for (int k=0; k<Ncont+L; k++)	{
				for (int l=k+1; l<Ncont+L; l++)	{
					os << '\t' << "sigma_" << k << '_' << l;
				}
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		}
		}
		os << '\n';
		os.flush();
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		if (chronogram->CheckBounds())	{
			exit(1);
		}
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetRootAge();
		if (chronoprior == 1)	{
			os << '\t' << GetCoalRate();
		}
		else if (chronoprior == 2)	{
			os << '\t' << T0->val();
			os << '\t' << K0->val();
			os << '\t' << K1->val();
			os << '\t' << massext->val();
			os << '\t' << divrate->val();
			os << '\t' << extrate->val();
		}
		else if (chronoprior == 3)	{
			os << '\t' << bddivrate->val() << '\t' << mu->val() << '\t' << psi->val() << '\t' << rho->val();
		}
		else if (chronoprior == 4)	{
			os << '\t' << Chi->val();
			os << '\t' << Chi2->val();
		}
		if (! prior)	{
		if (clockprior == 1)	{
			os << '\t' << wnnucmeanrate->val() << '\t' << wnnucsigma->val();
			if (morpho)	{
				os << '\t' << wnmorphomeanrate->val() << '\t' << wnmorphosigma->val();
			}
		}
		else if (clockprior == 2)	{
			os << '\t' << wnnucsigma->val();
			if (morpho)	{
				os << '\t' << wnmorphosigma->val();
			}
		}
		if (scaling == 1)	{
			os << '\t' << scalefactor->val() << '\t' << scalerate->val();
		}
		else if (scaling == 2)	{
			os << '\t' << scalesigma->val();
		}
		os << '\t' << GetNucTotalLength();
		if (morpho)	{
			os << '\t' << GetMorphoTotalLength();
		}
		if (dim)	{
		if (!clampdiag)	{
			for (int k=0; k<Ncont+L; k++)	{
				for (int l=k+1; l<Ncont+L; l++)	{
					os << '\t' << (*sigma)[k][l];
				}
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		}
		}
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		if (chronoprior == 1)	{
			os << *coalrate << '\n';
		}
		else if (chronoprior == 2)	{
			os << *T0 << '\t' << *K0 << '\t' << *K1 << '\t' << *massext << '\t' << *divrate << '\t' << *extrate << '\n';
		}
		else if (chronoprior == 3)	{
			os << *bddivrate << '\t' << *mu << '\t' << *psi << '\t' << *rho << '\n';
		}
		else if (chronoprior == 4)	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		os << *chronogram << '\n';
		os << *chronogram->GetScale() << '\n';
		if (! prior)	{
		if (dim)	{
			os << *DiagArray << '\n';
			os << *sigma << '\n';
			if (brownian)	{
				os << *brownianProcess->GetPureBrownianProcess() << '\n';
				os << *brownianProcess->GetInstantProcess() << '\n';
			}
			else	{
				os << *process << '\n';
			}
		}
		if (clockprior == 1)	{
			os << *wnnucmeanrate << '\n';
			os << *wnnucsigma << '\n';
			os << *wnnucratetree << '\n';
			if (morpho)	{
				os << *wnmorphomeanrate << '\n';
				os << *wnmorphosigma << '\n';
				os << *wnmorphoratetree << '\n';
			}
		}
		else if (clockprior == 2)	{
			os << *wnnucsigma << '\n';
			os << *wnnucratetree << '\n';
			if (morpho)	{
				os << *wnmorphosigma << '\n';
				os << *wnmorphoratetree << '\n';
			}
		}
		if (scaling == 1)	{
			os << *scalefactor << '\t' << *scalerate << '\n';
		}
		else if (scaling == 2)	{
			os << *scalesigma << '\n';
			os << *gammascaletree << '\n';
		}
		os << *nucrelrate << '\n';
		os << *nucstationary << '\n';
		}
	}

	void FromStream(istream& is)	{
		if (chronoprior == 1)	{
			is >> *coalrate;
		}
		else if (chronoprior == 2)	{
			is >> *T0 >> *K0 >> *K1 >> *massext >> *divrate >> *extrate;
		}
		else if (chronoprior == 3)	{
			is >> *bddivrate >> *mu >> *psi >> *rho;
		}
		else if (chronoprior == 4)	{
			is >> *Chi >> *Chi2;
		}
		is >> *chronogram;
		is >> *chronogram->GetScale();
		if (! prior)	{
		if (dim)	{
			is >> *DiagArray;
			is >> *sigma;
			if (brownian)	{
				is >> *brownianProcess->GetPureBrownianProcess();
				is >> *brownianProcess->GetInstantProcess();
			}
			else	{
				is >> *process;
			}
		}
		if (clockprior == 1)	{
			is >> *wnnucmeanrate;
			is >> *wnnucsigma;
			is >> *wnnucratetree;
			if (morpho)	{
				is >> *wnmorphomeanrate;
				is >> *wnmorphosigma;
				is >> *wnmorphoratetree;
			}
		}
		else if (clockprior == 2)	{
			is >> *wnnucsigma;
			is >> *wnnucratetree;
			if (morpho)	{
				is >> *wnmorphosigma;
				is >>*wnmorphoratetree;
			}
		}
		if (scaling == 1)	{
			is >> *scalefactor >> *scalerate;
		}
		else if (scaling == 2)	{
			is >> *scalesigma;
			is >> *gammascaletree;
		}
		is >> *nucrelrate;
		is >> *nucstationary;
		}
	}

	void SimulateFromFixedParam(string name, double minlength, double maxlength, double mintotallength, double maxtotallength, string treefile, string paramfile, int i)	{

		if (morpho)	{
			cerr << "in simu: morpho not yet implemented\n";
			exit(1);
		}
		// set diag elements

		ostringstream s;
		s << treefile << i << ".treeparam";
		ifstream tis(s.str().c_str());
		tis >> *chronogram;
		tis >> *chronogram->GetScale();
		chronogram->specialUpdate();
		cerr << "tree depth : " << chronogram->GetScale()->val() << '\n';

		ifstream is(paramfile.c_str());

		if (dim)	{

			string tmp;
			is >> tmp;
			if (tmp != "covariances")	{
				cerr << "error: param file has not the correct format\n";
				exit(1);
			}

			is >> *sigma;
		}

		if (clockprior)	{
			is >> *wnnucsigma;
		}
		is >> *nucrelrate;
		is >> *nucstationary;

		Update();

		do	{
		if (dim)	{
			for (int i=0; i<Ncont; i++)	{
				process->Unclamp(i+L);
			}
			process->ClampRoot();

			if (brownian)	{
				brownianProcess->Sample();
			}
			else	{
				process->Sample();
			}
		}

		if (clockprior)	{

			wnnucratetree->Sample();
		}

		FastUpdate();
		cerr << GetNucTotalLength() << '\t' << GetMinRate() << '\t' << GetMaxRate() << '\t' << GetMinLength() << '\t' << GetMaxLength() << '\n';

		}
		while ((GetNucTotalLength() > maxtotallength) || (GetNucTotalLength() < mintotallength) || (GetMinLength() < minlength) || (GetMaxLength() > maxlength));

		Update();
		cerr << "after update: " << GetNucTotalLength() << '\n';

		SequenceAlignment* datacopy = new SequenceAlignment(nucdata);

		nucphyloprocess->PostPredSample();

		nucphyloprocess->GetLeafData(datacopy);
		ofstream os((name + ".ali").c_str());

		datacopy->ToStream(os);
		if (contdata)	{
			ContinuousData* contdatacopy = new ContinuousData(contdata);
			process->GetLeafData(contdatacopy,L);
			ofstream cos((name + ".cont").c_str());
			contdatacopy->ToStream(cos);

			delete contdatacopy;
		}
		ofstream pros((name + ".fullparam").c_str());
		ToStream(pros);

		delete datacopy;
	}

	// void Simulate(string name, double minlength, double maxlength)	{
	void Simulate(string name, double minlength, double maxlength, double mintotallength, double maxtotallength)	{

		if (morpho)	{
			cerr << "in simu: morpho not yet implemented\n";
			exit(1);
		}
		// set diag elements

		MeanChronogram* meanchrono = new MeanChronogram(GetTree(),false,true,false);
		meanchrono->Add(GetChronogram());
		meanchrono->Normalise();
		ofstream tos((name + ".tree").c_str());
		meanchrono->ToStream(tos);

		ofstream dos((name + ".dates.tab").c_str());
		meanchrono->Tabulate(dos);

		ofstream vos((name + ".param").c_str());

		if (dim)	{
			vos << *sigma << '\n';
			vos << '\n';
			sigma->PrintCorrelationCoefficients(vos);
			vos << '\n';
		}

		if (clockprior)	{
			vos << *wnnucsigma << '\n';
		}

		do	{
		if (dim)	{
			for (int i=0; i<Ncont; i++)	{
				process->Unclamp(i+L);
			}
			process->ClampRoot();

			if (brownian)	{
				brownianProcess->Sample();
			}
			else	{
				process->Sample();
			}
		}

		if (clockprior)	{

			wnnucratetree->Sample();
		}

		FastUpdate();
		cerr << GetNucTotalLength() << '\t' << GetMinRate() << '\t' << GetMaxRate() << '\t' << GetMinLength() << '\t' << GetMaxLength() << '\n';
		// cerr << GetNucTotalLength() << '\n';

		}
		// while ((GetNucTotalLength() > maxlength) || (GetNucTotalLength() < minlength));
		while ((GetNucTotalLength() > maxtotallength) || (GetNucTotalLength() < mintotallength) || (GetMinLength() < minlength) || (GetMaxLength() > maxlength));

		Update();
		cerr << "after update: " << GetNucTotalLength() << '\n';

		SequenceAlignment* datacopy = new SequenceAlignment(nucdata);

		nucphyloprocess->PostPredSample();

		nucphyloprocess->GetLeafData(datacopy);
		ofstream os((name + ".ali").c_str());

		datacopy->ToStream(os);
		if (contdata)	{
			ContinuousData* contdatacopy = new ContinuousData(contdata);
			process->GetLeafData(contdatacopy,L);
			ofstream cos((name + ".cont").c_str());
			contdatacopy->ToStream(cos);

			delete contdatacopy;
		}
		ofstream pros((name + ".fullparam").c_str());
		ToStream(pros);

		delete datacopy;
	}

};

