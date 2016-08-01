
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"

#include "BDCalibratedChronogram.h"
#include "SerialCoalescent.h"
#include "LogisticSerialCoalescent.h"
#include "SerialBirthDeath.h"
#include "LogisticSerialBirthDeath.h"
#include "WhiteNoise.h"

#include "AlphaStableProcess.h"

#include "BranchProductProcess.h"
#include "MultiVariateUgamCompMove.h"

#include "F81TransitionMatrix.h"
#include "GeneralConjugatePath.h"

#include "MeanChronogram.h"

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
	TaxonSet* taxonset;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;

	double rootalpha;
	double rootbeta;

	// chronogram
	CalibratedChronogram* chronogram;

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

	Gamma* lnnucsigma;
	LogNormalTreeProcess* lnnucratetree;

	Gamma* lnmorphosigma;
	LogNormalTreeProcess* lnmorphoratetree;

	LengthTree* nucratetree;
	LengthTree* morphoratetree;

	BranchProductProcess* productnucratetree;
	BranchProductProcess* productmorphoratetree;

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

	PosUniform* asrho;
	// Gamma* asrho;
	Gamma* sigma;
	Gamma* beta;
	PosUniform* alpha;
	PosUniform* bralpha;
	AlphaStableProcess* alphastableprocess;
	AlphaStableProcess* brownianalphastableprocess;
	AlphaStableExpIntegralTree* asnucratetree;
	DoubleAlphaStableExpIntegralTree* doubleasnucratetree;

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
	int clockprior;
	double fixalpha;

	int nSegments;

	int morpho;

	int prior;
	int clamptree;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRLogNormalModel(string nucdatafile, string morpho2datafile, string morpho3datafile, string morpho4datafile, string treefile, string calibfile, double inrootage, double inrootstdev, double divrateval, double extrateval, double massextval, double K0val, double K1val, double divratestdev, double extratestdev, double massextstdev, double K0stdev, double K1stdev, double T1val, double bddivratemean, double mumean, double psimean, double rhomean, double bddivratestdev, double mustdev, double psistdev, double rhostdev, double divcutoff, int Nextant, int inchronoprior, int inclockprior, double ratemean, double ratestdev, double morphoratemean, double morphoratestdev, int inclamptree, int inprior, int innSegments, double infixalpha, bool sample = true)	{

		prior = inprior;
		clamptree = inclamptree;

		conjpath = true;
		chronoprior = inchronoprior;
		clockprior = inclockprior;
		nSegments = innSegments;
		fixalpha = infixalpha;

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

		cerr << "tree and data ok\n";

		calibset = 0;
		if (calibfile != "null")	{
			calibset = new FileCalibrationSet(calibfile, tree);
		}
		else	{
			calibset = new CalibrationSet(tree);
		}
		cerr << "calib ok\n";

		// ----------
		// construction of the graph
		// ----------

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
			if (clockprior == 0)	{
				// a log normal process on that tree
				lnnucsigma = new Gamma(One,One);
				lnnucratetree = new LogNormalTreeProcess(chronogram,lnnucsigma,INTEGRAL);
				nucratetree = lnnucratetree;

				if (morpho)	{
					lnmorphosigma = new Gamma(One,One);
					lnmorphoratetree = new LogNormalTreeProcess(chronogram,lnmorphosigma,INTEGRAL);
					morphoratetree = lnmorphoratetree;
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

				lnnucsigma = new Gamma(One,One);
				lnnucratetree = new LogNormalTreeProcess(chronogram,lnnucsigma,INTEGRAL);

				wnnucsigma = new Gamma(One,One);
				wnnucratetree = new GammaWhiteNoiseProcess(chronogram,wnnucsigma,0,MEAN,true);
				productnucratetree = new BranchProductProcess(lnnucratetree,wnnucratetree);
				nucratetree = productnucratetree;

				if (morpho)	{
					lnmorphosigma = new Gamma(One,One);
					lnmorphoratetree = new LogNormalTreeProcess(chronogram,lnmorphosigma,INTEGRAL);

					wnmorphosigma = new Gamma(One,One);
					wnmorphoratetree = new GammaWhiteNoiseProcess(chronogram,wnmorphosigma,0,MEAN,true);
					productmorphoratetree = new BranchProductProcess(lnmorphoratetree,wnmorphoratetree);
					morphoratetree = productmorphoratetree;
				}
			}
			else if (clockprior == 3)	{

				alpha = new PosUniform(One,2);
				if (fixalpha)	{
					alpha->ClampAt(fixalpha);
				}
				else	{
					alpha->setval(1.0);
				}
				beta = new Gamma(One,One);
				alphastableprocess = new AlphaStableProcess(chronogram,alpha,beta,nSegments);
				asnucratetree = new AlphaStableExpIntegralTree(chronogram,alphastableprocess);
				nucratetree = asnucratetree;

				if (morpho)	{
					cerr << "alpha stable morpho not yet implemented\n";
					exit(1);
				}
			}
			else if (clockprior == 4)	{

				alpha = new PosUniform(One,2);
				if (fixalpha)	{
					alpha->ClampAt(fixalpha);
				}
				else	{
					alpha->setval(1.0);
				}
				bralpha = new PosUniform(One,2);
				bralpha->ClampAt(2.0);
				beta = new Gamma(One,One);
				lnnucsigma = new Gamma(One,One);
				alphastableprocess = new AlphaStableProcess(chronogram,alpha,beta,nSegments);
				brownianalphastableprocess = new AlphaStableProcess(chronogram,bralpha,lnnucsigma,nSegments);
				doubleasnucratetree = new DoubleAlphaStableExpIntegralTree(chronogram,alphastableprocess,brownianalphastableprocess);
				nucratetree = doubleasnucratetree;
			}
			else if (clockprior == 5)	{

				alpha = new PosUniform(One,2);
				if (fixalpha)	{
					alpha->ClampAt(fixalpha);
				}
				else	{
					alpha->setval(1.0);
				}
				beta = new Gamma(One,One);
				alphastableprocess = new AlphaStableProcess(chronogram,alpha,beta,nSegments);
				asnucratetree = new AlphaStableExpIntegralTree(chronogram,alphastableprocess);

				if (morpho)	{
					cerr << "alpha stable morpho not yet implemented\n";
					exit(1);
				}

				wnnucsigma = new Gamma(One,One);
				wnnucratetree = new GammaWhiteNoiseProcess(chronogram,wnnucsigma,0,MEAN,true);
				productnucratetree = new BranchProductProcess(asnucratetree,wnnucratetree);
				nucratetree = productnucratetree;
			}
			else if (clockprior == 6)	{

				alpha = new PosUniform(One,2);
				alpha->ClampAt(2.0);
				asrho = new PosUniform(One,10);
				sigma = new Gamma(One,One);
				beta = new Gamma(One,One);
				alphastableprocess = new AlphaStableProcess(chronogram,alpha,beta,nSegments,asrho,sigma);
				asnucratetree = new AlphaStableExpIntegralTree(chronogram,alphastableprocess);
				nucratetree = asnucratetree;

				if (morpho)	{
					cerr << "alpha stable morpho not yet implemented\n";
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
				nucpathconjtree = 0;
				if (morpho)	{
					morpho2pathconjtree = 0;
					morpho3pathconjtree = 0;
					morpho4pathconjtree = 0;
				}
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
			if (clockprior == 0)	{
				RootRegister(lnnucratetree->GetRootRate());
				if (morpho)	{
					RootRegister(lnmorphoratetree->GetRootRate());
				}
			}
			else if (clockprior == 1)	{
				RootRegister(rateAlpha);
				RootRegister(rateBeta);
				if (morpho)	{
					RootRegister(morphorateAlpha);
					RootRegister(morphorateBeta);
				}
			}
			else if (clockprior == 2)	{
				RootRegister(lnnucratetree->GetRootRate());
				if (morpho)	{
					RootRegister(lnmorphoratetree->GetRootRate());
				}
			}
			RootRegister(nucrelrate);
			RootRegister(nucstationary);
			if (morpho)	{
				RootRegister(morpho2stationary);
				RootRegister(morpho3stationary);
				RootRegister(morpho4stationary);
			}
		}
		Register();

		MakeScheduler();

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
		if (clockprior >= 3)	{
			asnucratetree->specialUpdate();
			productnucratetree->specialUpdate();
		}
	}

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
			if (clockprior == 0)	{
				total += lnnucsigma->GetLogProb();
				total += lnnucratetree->GetLogProb();
				if (morpho)	{
					total += lnmorphosigma->GetLogProb();
					total += lnmorphoratetree->GetLogProb();
				}
			}

			else if (clockprior == 1)	{
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
				total += lnnucsigma->GetLogProb();
				total += lnnucratetree->GetLogProb();
				total += wnnucsigma->GetLogProb();
				total += wnnucratetree->GetLogProb();
				if (morpho)	{
					total += lnmorphosigma->GetLogProb();
					total += lnmorphoratetree->GetLogProb();
					total += wnmorphosigma->GetLogProb();
					total += wnmorphoratetree->GetLogProb();
				}
			}

			else if (clockprior == 3)	{
				total += beta->GetLogProb();
				total += alphastableprocess->GetLogProb();
			}

			else if (clockprior == 4)	{
				total += beta->GetLogProb();
				total += lnnucsigma->GetLogProb();
				total += alphastableprocess->GetLogProb();
				total += brownianalphastableprocess->GetLogProb();
			}

			else if (clockprior == 5)	{
				total += beta->GetLogProb();
				total += alphastableprocess->GetLogProb();
				total += wnnucsigma->GetLogProb();
				total += wnnucratetree->GetLogProb();
			}

			else if (clockprior == 3)	{
				total += beta->GetLogProb();
				total += asrho->GetLogProb();
				total += sigma->GetLogProb();
				total += alphastableprocess->GetLogProb();
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

	LogNormalTreeProcess* GetLnNucRateTree()	{
		if ((clockprior != 0) && (clockprior != 2))	{
			cerr << "error: no log normal tree process\n";
			exit(1);
		}
		return lnnucratetree;
	}

	GammaWhiteNoiseProcess* GetWnNucRateTree()	{
		if ((clockprior != 1) && (clockprior != 2) && (clockprior != 5))	{
			cerr << "error: no white noise process\n";
			exit(1);
		}
		return wnnucratetree;
	}

	LengthTree* GetNucRateTree()	{return nucratetree;}
	LengthTree* GetMorphoRateTree()	{return morphoratetree;}

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
			scheduler.Register(new AllInternalNodesMove(chronogram,0.001),10,"all int chrono");
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

			if (! prior)	{
				if (clockprior == 0)	{

					scheduler.Register(new SimpleMove(lnnucratetree,1),10,"lognormal");
					scheduler.Register(new SimpleMove(lnnucratetree,0.1),10,"lognormal");
					scheduler.Register(new SimpleMove(lnnucratetree,0.01),10,"lognormal");

					scheduler.Register(new SimpleMove(lnnucsigma,10),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,1),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,0.1),10,"sigma");

					if (morpho)	{
						scheduler.Register(new SimpleMove(lnmorphoratetree,1),10,"morpho lognormal");
						scheduler.Register(new SimpleMove(lnmorphoratetree,0.1),10,"morpho lognormal");
						scheduler.Register(new SimpleMove(lnmorphoratetree,0.01),10,"morpho lognormal");

						scheduler.Register(new SimpleMove(lnmorphosigma,10),10,"morphosigma");
						scheduler.Register(new SimpleMove(lnmorphosigma,1),10,"morpho sigma");
						scheduler.Register(new SimpleMove(lnmorphosigma,0.1),10,"morpho sigma");
					}

				}
				else if (clockprior == 1)	{

					scheduler.Register(new SimpleMove(wnnucratetree,1),10,"white noise");
					scheduler.Register(new SimpleMove(wnnucratetree,0.1),10,"white noise");
					scheduler.Register(new SimpleMove(wnnucratetree,0.01),10,"white noise");

					scheduler.Register(new SimpleMove(wnnucmeanrate,1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucmeanrate,0.1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucmeanrate,0.01),10,"mean rate");

					scheduler.Register(new SimpleMove(wnnucsigma,1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucsigma,0.1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucsigma,0.01),10,"mean rate");

					if (morpho)	{
						scheduler.Register(new SimpleMove(wnmorphoratetree,1),10,"morpho white noise");
						scheduler.Register(new SimpleMove(wnmorphoratetree,0.1),10,"morpho white noise");
						scheduler.Register(new SimpleMove(wnmorphoratetree,0.01),10,"morpho white noise");

						scheduler.Register(new SimpleMove(wnmorphomeanrate,1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphomeanrate,0.1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphomeanrate,0.01),10,"morpho mean rate");

						scheduler.Register(new SimpleMove(wnmorphosigma,1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphosigma,0.1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphosigma,0.01),10,"morpho mean rate");
					}
				}
				else if (clockprior == 2)	{

					scheduler.Register(new SimpleMove(lnnucratetree,1),10,"lognormal");
					scheduler.Register(new SimpleMove(lnnucratetree,0.1),10,"lognormal");
					scheduler.Register(new SimpleMove(lnnucratetree,0.01),10,"lognormal");

					scheduler.Register(new SimpleMove(lnnucsigma,10),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,1),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,0.1),10,"sigma");

					scheduler.Register(new SimpleMove(wnnucratetree,1),10,"white noise");
					scheduler.Register(new SimpleMove(wnnucratetree,0.1),10,"white noise");
					scheduler.Register(new SimpleMove(wnnucratetree,0.01),10,"white noise");

					scheduler.Register(new SimpleMove(wnnucsigma,1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucsigma,0.1),10,"mean rate");
					scheduler.Register(new SimpleMove(wnnucsigma,0.01),10,"mean rate");

					if (morpho)	{
						scheduler.Register(new SimpleMove(lnmorphoratetree,1),10,"morpho lognormal");
						scheduler.Register(new SimpleMove(lnmorphoratetree,0.1),10,"morpho lognormal");
						scheduler.Register(new SimpleMove(lnmorphoratetree,0.01),10,"morpho lognormal");

						scheduler.Register(new SimpleMove(lnmorphosigma,10),10,"morphosigma");
						scheduler.Register(new SimpleMove(lnmorphosigma,1),10,"morpho sigma");
						scheduler.Register(new SimpleMove(lnmorphosigma,0.1),10,"morpho sigma");

						scheduler.Register(new SimpleMove(wnmorphoratetree,1),10,"morpho white noise");
						scheduler.Register(new SimpleMove(wnmorphoratetree,0.1),10,"morpho white noise");
						scheduler.Register(new SimpleMove(wnmorphoratetree,0.01),10,"morpho white noise");

						scheduler.Register(new SimpleMove(wnmorphosigma,1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphosigma,0.1),10,"morpho mean rate");
						scheduler.Register(new SimpleMove(wnmorphosigma,0.01),10,"morpho mean rate");
					}

					scheduler.Register(new UniVariateUgamCompMove(wnnucratetree,lnnucratetree,1),10,"wnprocess comp");
					scheduler.Register(new UniVariateUgamCompMove(wnnucratetree,lnnucratetree,0.1),10,"wnprocess comp");
					scheduler.Register(new UniVariateUgamCompMove(wnnucratetree,lnnucratetree,0.01),10,"wnprocess comp");

					if (morpho)	{
						scheduler.Register(new UniVariateUgamCompMove(wnmorphoratetree,lnmorphoratetree,1),10,"wnprocess comp");
						scheduler.Register(new UniVariateUgamCompMove(wnmorphoratetree,lnmorphoratetree,0.1),10,"wnprocess comp");
						scheduler.Register(new UniVariateUgamCompMove(wnmorphoratetree,lnmorphoratetree,0.01),10,"wnprocess comp");
					}

				}

				else if ((clockprior == 3) || (clockprior == 5) || (clockprior == 6))	{

					scheduler.Register(new SimpleMove(alphastableprocess,1),10,"as branch");
					scheduler.Register(new SimpleMove(alphastableprocess,0.1),10,"as branch");
					scheduler.Register(new SimpleMove(alphastableprocess,0.01),10,"as branch");

					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,1),10,"as node");
					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,0.1),10,"as node");
					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,0.01),10,"as node");

					scheduler.Register(new SimpleMove(beta,10),10,"beta");
					scheduler.Register(new SimpleMove(beta,1),10,"beta");
					scheduler.Register(new SimpleMove(beta,0.1),10,"beta");

					/*
					scheduler.Register(new AlphaStableAlphaRescaleMove(alphastableprocess,alpha,0.1),10,"as rescale");
					scheduler.Register(new AlphaStableAlphaRescaleMove(alphastableprocess,alpha,0.01),10,"as rescale");
					scheduler.Register(new AlphaStableAlphaRescaleMove(alphastableprocess,alpha,0.001),10,"as rescale");

					scheduler.Register(new AlphaStableAlphaRescaleBridgeMove(alphastableprocess,alpha,0.1),10,"as rescale bridge");
					scheduler.Register(new AlphaStableAlphaRescaleBridgeMove(alphastableprocess,alpha,0.01),10,"as rescale bridge");
					scheduler.Register(new AlphaStableAlphaRescaleBridgeMove(alphastableprocess,alpha,0.001),10,"as rescale bridge");
					*/

					scheduler.Register(new AlphaStableRescaleMove(alphastableprocess,beta,0.1),10,"as rescale");
					scheduler.Register(new AlphaStableRescaleMove(alphastableprocess,beta,0.01),10,"as rescale");
					scheduler.Register(new AlphaStableRescaleMove(alphastableprocess,beta,0.001),10,"as rescale");

					scheduler.Register(new AlphaStableRescaleBridgeMove(alphastableprocess,beta,0.1),10,"as rescale bridge");
					scheduler.Register(new AlphaStableRescaleBridgeMove(alphastableprocess,beta,0.01),10,"as rescale bridge");
					scheduler.Register(new AlphaStableRescaleBridgeMove(alphastableprocess,beta,0.001),10,"as rescale bridge");

					if (clockprior == 5)	{

						scheduler.Register(new SimpleMove(wnnucratetree,1),10,"white noise");
						scheduler.Register(new SimpleMove(wnnucratetree,0.1),10,"white noise");
						scheduler.Register(new SimpleMove(wnnucratetree,0.01),10,"white noise");

						scheduler.Register(new SimpleMove(wnnucsigma,1),10,"mean rate");
						scheduler.Register(new SimpleMove(wnnucsigma,0.1),10,"mean rate");
						scheduler.Register(new SimpleMove(wnnucsigma,0.01),10,"mean rate");

						scheduler.Register(new ASUgamCompMove(wnnucratetree,alphastableprocess,1),10,"wnprocess comp");
						scheduler.Register(new ASUgamCompMove(wnnucratetree,alphastableprocess,0.1),10,"wnprocess comp");
						scheduler.Register(new ASUgamCompMove(wnnucratetree,alphastableprocess,0.01),10,"wnprocess comp");

					}

					if (clockprior == 6)	{

						scheduler.Register(new SimpleMove(sigma,1),10,"sigma");
						scheduler.Register(new SimpleMove(sigma,0.1),10,"sigma");
						scheduler.Register(new SimpleMove(sigma,0.01),10,"sigma");

						scheduler.Register(new SimpleMove(asrho,1),10,"rho");
						scheduler.Register(new SimpleMove(asrho,0.1),10,"rho");
						scheduler.Register(new SimpleMove(asrho,0.01),10,"rho");

					}

					scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
					scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
					scheduler.Register(new SimpleMove(alpha,0.01),10,"alpha");

				}

				else if (clockprior == 4)	{

					scheduler.Register(new SimpleMove(alphastableprocess,1),10,"as branch");
					scheduler.Register(new SimpleMove(alphastableprocess,0.1),10,"as branch");
					scheduler.Register(new SimpleMove(alphastableprocess,0.01),10,"as branch");

					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,1),10,"as node");
					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,0.1),10,"as node");
					scheduler.Register(new AlphaStableNodeMove(alphastableprocess,0.01),10,"as node");

					scheduler.Register(new SimpleMove(brownianalphastableprocess,1),10,"brownian branch");
					scheduler.Register(new SimpleMove(brownianalphastableprocess,0.1),10,"brownian branch");
					scheduler.Register(new SimpleMove(brownianalphastableprocess,0.01),10,"brownian branch");

					scheduler.Register(new AlphaStableNodeMove(brownianalphastableprocess,1),10,"brownian node");
					scheduler.Register(new AlphaStableNodeMove(brownianalphastableprocess,0.1),10,"brownian node");
					scheduler.Register(new AlphaStableNodeMove(brownianalphastableprocess,0.01),10,"brownian node");

					scheduler.Register(new SimpleMove(beta,10),10,"beta");
					scheduler.Register(new SimpleMove(beta,1),10,"beta");
					scheduler.Register(new SimpleMove(beta,0.1),10,"beta");

					scheduler.Register(new SimpleMove(lnnucsigma,10),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,1),10,"sigma");
					scheduler.Register(new SimpleMove(lnnucsigma,0.1),10,"sigma");

					scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
					scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
					scheduler.Register(new SimpleMove(alpha,0.01),10,"alpha");

				}

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

	double GetTotalRateVariance()	{

		return GetVariance(nucratetree,true);
	}

	double GetLnRateVariance()	{

		if (clockprior >= 3)	{
			return GetVariance(asnucratetree,true);
		}
		return GetVariance(lnnucratetree,true);
	}

	double GetWnRateVariance()	{

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

	double GetNucTotalLength()	{
		return GetNucRateTree()->GetTotalLength();
	}

	double GetMorphoTotalLength()	{
		return GetMorphoRateTree()->GetTotalLength();
	}

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

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\trootage";
		if (chronoprior == 1)	{
			os << "\tcoalrate";
		}
		else if (chronoprior == 2)	{
			os << "\tT0\tK0\tK1\tmassext\tdivrate\textrate";
		}
		else if (chronoprior == 3)	{
			os << "\tbddivrate\tmu\tpsi\trho";
		}
		else if (chronoprior == 4)	{
			os << "\tchi\tchi2";
		}
		if (! prior)	{
			if (clockprior == 0)	{
				os << "\tlnsigma";
				if (morpho)	{
					os << "\tlnmorphosigma";
				}
			}
			else if (clockprior == 1)	{
				os << "\twnsigma\tmeanrate";
				if (morpho)	{
					os << "\twnmorphosigma\tmorphomeanrate";
				}
			}
			else if (clockprior == 2)	{
				os << "\tlnsigma\twnsigma";
				if (morpho)	{
					os << "\tlnmorphosigma\twnmorphosigma";
				}
			}
			else if (clockprior == 3)	{
				os << "\talpha\tbeta";
			}
			else if (clockprior == 4)	{
				os << "\talpha\tbeta\tsigma";
			}
			else if (clockprior == 5)	{
				os << "\talpha\tbeta\twnsigma";
			}
			else if (clockprior == 6)	{
				os << "\talpha\tbeta\trho\tsigma";
			}
				
		}
		os << '\n';
		os.flush();
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		chronogram->CheckBounds();
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetRootAge();
		if (chronoprior == 1)	{
			os << '\t' << GetCoalRate();
		}
		else if (chronoprior == 2)	{
			os << '\t' << GetT0();
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
			if (clockprior == 0)	{
				os << '\t' << lnnucsigma->val();
				if (morpho)	{
					os << '\t' << lnmorphosigma->val();
				}
			}
			else if (clockprior == 1)	{
				os << '\t' << wnnucsigma->val();
				os << '\t' << wnnucmeanrate->val();
				if (morpho)	{
					os << '\t' << wnmorphosigma->val();
					os << '\t' << wnmorphomeanrate->val();
				}
			}
			else if (clockprior == 2)	{
				os << '\t' << lnnucsigma->val();
				os << '\t' << wnnucsigma->val();
				if (morpho)	{
					os << '\t' << lnmorphosigma->val();
					os << '\t' << wnmorphosigma->val();
				}
			}
			else if (clockprior == 3)	{
				os << '\t' << alpha->val();
				os << '\t' << beta->val();
			}
			else if (clockprior == 4)	{
				os << '\t' << alpha->val();
				os << '\t' << beta->val();
				os << '\t' << lnnucsigma->val();
			}
			else if (clockprior == 5)	{
				os << '\t' << alpha->val();
				os << '\t' << beta->val();
				os << '\t' << wnnucsigma->val();
			}
			else if (clockprior == 6)	{
				os << '\t' << alpha->val();
				os << '\t' << beta->val();
				os << '\t' << asrho->val();
				os << '\t' << sigma->val();
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
			if (clockprior == 0)	{
				os << *lnnucsigma << '\n';
				os << *lnnucratetree << '\n';
			}
			else if (clockprior == 1)	{
				os << *wnnucsigma << '\n';
				os << *wnnucmeanrate << '\n';
				os << *wnnucratetree << '\n';
			}
			else if (clockprior == 2)	{
				os << *lnnucsigma << '\n';
				os << *lnnucratetree << '\n';
				os << *wnnucsigma << '\n';
				os << *wnnucratetree << '\n';
			}
			else if (clockprior == 3)	{
				os << *alpha << '\n';
				os << *beta << '\n';
				os << *alphastableprocess << '\n';
			}
			else if (clockprior == 4)	{
				os << *alpha << '\n';
				os << *beta << '\n';
				os << *lnnucsigma << '\n';
				os << *alphastableprocess << '\n';
				os << *brownianalphastableprocess << '\n';
			}
			else if (clockprior == 5)	{
				os << *alpha << '\n';
				os << *beta << '\n';
				os << *alphastableprocess << '\n';
				os << *wnnucsigma << '\n';
				os << *wnnucratetree << '\n';
			}
			else if (clockprior == 6)	{
				os << *alpha << '\n';
				os << *beta << '\n';
				os << *asrho << '\n';
				os << *sigma << '\n';
				os << *alphastableprocess << '\n';
			}
			if (morpho)	{
				if (clockprior == 0)	{
					os << *lnmorphosigma << '\n';
					os << *lnmorphoratetree << '\n';
				}
				else if (clockprior == 1)	{
					os << *wnmorphosigma << '\n';
					os << *wnmorphomeanrate << '\n';
					os << *wnmorphoratetree << '\n';
				}
				else if (clockprior == 2)	{
					os << *lnmorphosigma << '\n';
					os << *lnmorphoratetree << '\n';
					os << *wnmorphosigma << '\n';
					os << *wnmorphoratetree << '\n';
				}
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
			if (clockprior == 0)	{
				is >> *lnnucsigma;
				is >> *lnnucratetree;
			}
			else if (clockprior == 1)	{
				is >> *wnnucsigma;
				is >> *wnnucmeanrate;
				is >> *wnnucratetree;
			}
			else if (clockprior == 2)	{
				is >> *lnnucsigma;
				is >> *lnnucratetree;
				is >> *wnnucsigma;
				is >> *wnnucratetree;
			}
			else if (clockprior == 3)	{
				is >> *alpha;
				is >> *beta;
				is >> *alphastableprocess;
			}
			else if (clockprior == 4)	{
				is >> *alpha;
				is >> *beta;
				is >> *lnnucsigma;
				is >> *alphastableprocess;
				is >> *brownianalphastableprocess;
			}
			else if (clockprior == 5)	{
				is >> *alpha;
				is >> *beta;
				is >> *alphastableprocess;
				is >> *wnnucsigma;
				is >> *wnnucratetree;
			}
			else if (clockprior == 6)	{
				is >> *alpha;
				is >> *beta;
				is >> *asrho;
				is >> *sigma;
				is >> *alphastableprocess;
			}
			if (morpho)	{
				if (clockprior == 0)	{
					is >> *lnmorphosigma;
					is >> *lnmorphoratetree;
				}
				else if (clockprior == 1)	{
					is >> *wnmorphosigma;
					is >> *wnmorphomeanrate;
					is >> *wnmorphoratetree;
				}
				else if (clockprior == 2)	{
					is >> *lnmorphosigma;
					is >> *lnmorphoratetree;
					is >> *wnmorphosigma;
					is >> *wnmorphoratetree;
				}
			}
			is >> *nucrelrate;
			is >> *nucstationary;
		}
	}

	void Simulate(string name, double minlength, double maxlength)	{

		if (morpho)	{
			cerr << "in simu: morpho not yet implemented\n";
			exit(1);
		}
		// set diag elements

		MeanChronogram* meanchrono = new MeanChronogram(GetTree());
		meanchrono->Add(chronogram);
		meanchrono->Normalise();
		ofstream tos((name + ".tree").c_str());
		meanchrono->ToStream(tos);

		ofstream dos((name + ".dates.tab").c_str());
		meanchrono->Tabulate(dos);

		ofstream vos((name + ".param").c_str());
		vos << *wnnucsigma << '\n';
		vos << *wnnucmeanrate << '\n';

		do	{
			if (clockprior == 0)	{
				lnnucratetree->Sample();
			}
			else if (clockprior == 1)	{
				wnnucratetree->Sample();
			}
			// FastUpdate();
			cerr << GetNucTotalLength() << '\n';
		}
		while ((GetNucTotalLength() > maxlength) || (GetNucTotalLength() < minlength));

		Update();
		cerr << "after update: " << GetNucTotalLength() << '\n';

		SequenceAlignment* datacopy = new SequenceAlignment(nucdata);

		nucphyloprocess->PostPredSample();

		nucphyloprocess->GetLeafData(datacopy);
		ofstream os((name + ".ali").c_str());

		datacopy->ToStream(os);
		ofstream pros((name + ".param").c_str());
		ToStream(pros);

		delete datacopy;
	}
};

