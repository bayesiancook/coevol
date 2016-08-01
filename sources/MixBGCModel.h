
#ifndef BGC
#define BGC

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "AutoRegressiveNormalTreeProcess.h"
#include "CodonSequenceAlignment.h"
#include "RYSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
// #include "CoalCalibratedChronogram.h"

#include "BranchProcess.h"
#include "NodeProcess.h"
#include "MatrixTree.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"

#include "WhiteNoise.h"

#include "GeneralConjugatePath.h"

#include "Jeffreys.h"

#include "Partition.h"
#include "SplitPartition.h"

#include "SplitLengthTree.h"
#include "SplitChronogram.h"
#include "SplitMultiVariateMove.h"
#include "MultiVariatePropagateMove.h"

#include "PartitionMultiVariateTreeProcess.h"

#include "LinRegContSub.h"

#include "BGCSubMatrix.h"

#include "BGC.h"
#include "BGCCpGSubMatrix.h"


#include "TripletSequenceAlignment.h"

class BGCModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;
	FileSequenceAlignment** nucdata;
	SequenceAlignment** data;
	TripletSequenceAlignment** tripletdata;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	bool triplet;

	// number of columns
	int Ngene;
	int* Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	// int Nstate;

	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;
	Const<PosReal>* Ten;
	Const<PosReal>* Tenth;

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

	bool iscalib;

	// chronogram
	Chronogram* chronogram;

	LengthTree* lengthtree;
	NodeBranchVarTree<PosReal,PosReal>* finalchrono;

	Jeffreys* synsigma;
	LogNormalTreeProcess* synratetree;

	JeffreysIIDArray* ContDiagArray;
	SigmaZero* ContSigmaZero;
	Rvar<CovMatrix>* contsigma;

	JeffreysIIDArray* BGCDiagArray;
	SigmaZero* BGCSigmaZero;
	DiagonalCovMatrix* bgcsigma;

	IIDNormalArray* regarray;
	//Gamma* phi;
	Jeffreys* phi;

	Jeffreys* alpha;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* contprocess;
	LinRegNormalProcess* logbgcprocess;
	Normal* logbgcoffset;
	BGCProcess* bgcprocess;


	Jeffreys* alpha1;
	Jeffreys* beta1;
	Jeffreys* alpha2;
	// Jeffreys* beta2;

	Gamma** rootbgc;
	Gamma** rootbgcoffset;

	Jeffreys* lambda;
	Jeffreys* at2cg;
	Jeffreys* at2gc;
	Jeffreys* at2ta;
	Jeffreys* gc2cg;

	Jeffreys* cg2at;
	Jeffreys* cg2gc;
	Jeffreys* cg2ta;

	Jeffreys* ctlambda;
	Jeffreys* ctat2cg;
	Jeffreys* ctat2gc;
	Jeffreys* ctat2ta;
	Jeffreys* ctgc2cg;

	Jeffreys* ctcg2at;
	Jeffreys* ctcg2gc;
	Jeffreys* ctcg2ta;

	Jeffreys* cpgrate;

	MeanExpTreeFromMultiVariate* lambdaprocess;

	int discn;
	int discgam;

	// Jeffreys** rootbgc;

	// Gamma* rho;
	// Gamma* genealpha;
	Jeffreys* rho;
	Jeffreys* rhosigma;
	LogNormalTreeProcess* rhoprocess;
	Jeffreys* genealpha;
	Jeffreys* genealpha2;
	RandomMutSubMatrix* mutmatrix;
	RandomNonRevMutSubMatrix* nonrevmutmatrix;
	RandomMutSubMatrix* ctmutmatrix;
	RandomNonRevMutSubMatrix* ctnonrevmutmatrix;

	NucMatrixTree* mutmatrixtree;

	// 0 : ar
	// 1 : gamma
	// 2 : wn
	int geneprocesstype;
	LengthTree** geneprocess;
	// LengthTree** geneprocess;
	AutoRegressiveLogNormalTreeProcess** argeneprocess;
	LogNormalTreeProcess** lngeneprocess;
	GammaTree** gammageneprocess;
	GammaWhiteNoiseProcess** wngeneprocess;
	BranchProductProcess** mixedgeneprocess;
	BranchProductProcess** branchproductprocess;

	bool clampbgc;

	NucMatrixTree** nucmatrixtree;

	// phylo process
	BranchMatrixPathConjugateTree** pathconjtree;
	PhyloProcess** phyloprocess;

	bool meanexp;

	bool conjpath;
	bool priorsampling;

	bool normalise;

	int nrep;

	int df;
	int Ninterpol;

	int clampsuffstat;
	string suffstatfile;

	double fixalpha;
	bool withalphaprocess;
	bool withhotspot;
	bool clamphotspot;

	int timelambda;

	int Nreg;

	MeanExpTreeFromMultiVariate* alphaprocess;

	bool clamptree;
	bool clampdiag;
	bool autoregressive;

	bool clampreg;

	bool clampbgcoffset;

	int flexrho;

	public:

	SequenceAlignment* GetData(int gene)	{
		return data[gene];
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

	bool Reversible()	{
		return (discn == 0);
	}

	bool isCalibrated()	{
		return iscalib;
	}

	BGCModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, int indf, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inNinterpol, string insuffstatfile, string rootfile, double infixalpha, int inclampbgcoffset, int inflexrho, double inlambda, double inat2cg, double inat2gc, double ingc2cg, double incg2at, double incg2gc, double incg2ta, int indiscgam, int indiscn, bool intriplet, bool inclampreg, bool sample=true)	{

		clampreg = inclampreg;
		triplet = intriplet;

		discn = indiscn;
		discgam = indiscgam;

		if (inclampbgcoffset == 0)	{
			clampbgcoffset = false;
			clampbgc = false;
		}
		else if (inclampbgcoffset == 1)	{
			clampbgcoffset = true;
			clampbgc = false;
		}
		else if (inclampbgcoffset == 2)	{
			clampbgcoffset = false;
			clampbgc = true;
		}
		else if (inclampbgcoffset == 3)	{
			clampbgcoffset = true;
			clampbgc = true;
		}

		if (inflexrho == -1)	{
			geneprocesstype = 1;
		}
		else if (inflexrho == -2)	{
			geneprocesstype = 2;
		}
		else if (inflexrho == -3)	{
			geneprocesstype = 3;
		}
		else if (inflexrho == -4)	{
			geneprocesstype = 4;
		}
		else	{
			flexrho = inflexrho;
		}

		clamptree = inclamptree;
		fixalpha = infixalpha;
		clamphotspot = false;
		if ((fixalpha == -3) || (fixalpha == -4))	{
			if (fixalpha == -4)	{
				clamphotspot = true;
			}
			fixalpha = -1;
			withalphaprocess = true;
			withhotspot = true;
			Nreg = 2;
		}
		else if (fixalpha == -2)	{
			fixalpha = -1;
			withalphaprocess = true;
			withhotspot = false;
			Nreg = 2;
		}
		else	{
			withalphaprocess = false;
			Nreg = 1;
		}

		timelambda = 0;
		if (inlambda == -2)	{
			inlambda = -1;
			timelambda = 1;
			Nreg ++;
		}

		autoregressive = inautoregressive;
		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		Ninterpol = inNinterpol;
		df = indf;

		chronoprior = inchronoprior;
		iscalib = false;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		meanexp = inmeanexp;

		// get data from file

		ifstream is(datafile.c_str());
		is >> Ngene;
		Nsite = new int[Ngene];
		nucdata = new FileSequenceAlignment*[Ngene];
		data = new SequenceAlignment*[Ngene];
		tripletdata = 0;
		if (triplet)	{
			tripletdata = new TripletSequenceAlignment*[Ngene];
		}
		for (int gene=0; gene<Ngene; gene++)	{
			string filename;
			is >> filename;
			nucdata[gene] = new FileSequenceAlignment(filename);
			if (triplet)	{
				tripletdata[gene] = new TripletSequenceAlignment(nucdata[gene]);
				data[gene] = tripletdata[gene];
			}
			else	{
				data[gene] = nucdata[gene];
			}
			Nsite[gene] = GetData(gene)->GetNsite();
			// Nstate = GetData(gene)->GetNstate();
			cerr << filename << '\t' << Nsite[gene] << '\n';
		}

		priorsampling = false;

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else if (inconjpath == 2)	{
			conjpath = false;
			priorsampling = true;
		}
		else	{
			conjpath = inconjpath;
		}
		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 10 : 1;
		}
		cerr << "nrep : " << nrep << '\n';
		cerr << conjpath << '\n';
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		taxonset = nucdata[0]->GetTaxonSet();

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
		Ten = new Const<PosReal>(10);
		Tenth = new Const<PosReal>(0.1);

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		if (calibfile != "None")	{
			iscalib = true;

			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			if (rootage == -1)	{
				a = b = -1;
			}
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
			}
			else {
				cerr << "BD\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,chronoprior);
			}
		}
		else	{
			chronogram = new Chronogram(tree,One);
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

		double mindiag = 0.001;
		double maxdiag = 1000;

		synsigma = new Jeffreys(mindiag,maxdiag,Zero);
		synsigma->setval(1);
		synratetree = new LogNormalTreeProcess(lengthtree,synsigma,INTEGRAL);
		synratetree->Reset();

		ContDiagArray = new JeffreysIIDArray(Ncont,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			ContDiagArray->setval(1.0);
		}
		else	{
			ContDiagArray->ClampAt(priorsigma);
		}

		ContSigmaZero = new SigmaZero(ContDiagArray);

		if (clampdiag)	{
			contsigma = new DiagonalCovMatrix(ContSigmaZero,Ncont+df);
		}
		else	{
			contsigma = new InverseWishartMatrix(ContSigmaZero,Ncont+df);
		}

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

		if (clampdiag)	{
			contprocess = new MultiVariateTreeProcess(contsigma,lengthtree,0,0,rootmean,rootvar);
		}
		else	{
			contprocess = new MultiVariateTreeProcess(contsigma,lengthtree,0,0,rootmean, rootvar);
		}

		BGCDiagArray = new JeffreysIIDArray(Nreg,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			BGCDiagArray->setval(1.0);
		}
		else	{
			BGCDiagArray->ClampAt(priorsigma);
		}

		BGCSigmaZero = new SigmaZero(BGCDiagArray);
		bgcsigma = new DiagonalCovMatrix(BGCSigmaZero,Nreg+df);

		double min = 0.001;
		double max = 1000;
		cerr << "regarray\n";
		regarray = new IIDNormalArray(Nreg,Ncont,Zero,One);
		if (autoregressive)	{
			phi = new Jeffreys(min,max,Zero);
			// phi = new Gamma(One,One);
		}
		else	{
			phi = 0;
		}
		cerr << Ncont << '\n';
		regarray->SetAtZero();
		if (clampreg)	{
			cerr << "clamp\n";
			regarray->ClampAtZero();
		}
		/*
		regarray->GetIIDNormal(0)->ClampAt(-1,1);
		regarray->GetIIDNormal(0)->ClampAt(1,2);
		*/
		// regarray->GetIIDNormal(0)->ClampAt(0,1);
		// regarray->GetIIDNormal(0)->ClampAt(0,2);

		if (withhotspot && clamphotspot)	{
			regarray->GetIIDNormal(0)->ClampAt(0,1);
			regarray->GetIIDNormal(0)->ClampAt(0,2);
			regarray->GetIIDNormal(1)->ClampAt(0,0);
		}
		cerr << "bgc process\n";

		logbgcprocess = new LinRegNormalProcess(lengthtree,contprocess,regarray,bgcsigma,phi);
		for (int i=0; i<Nreg; i++)	{
			logbgcprocess->CutOff(0.1,i);
			(*logbgcprocess->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode()))[i] = 0;
		}
		logbgcprocess->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode())->Clamp();

		alpha = 0;
		if (fixalpha != 0)	{
			double min = 0.0001;
			double max = 10;
			alpha = new Jeffreys(min,max,Zero);
			// alpha = new Gamma(One,One);
			cerr << *alpha << '\n';
			if (fixalpha != -1)	{
				alpha->ClampAt(fixalpha);
			}
		}
		alphaprocess = 0;
		if (withalphaprocess)	{
			alphaprocess = new MeanExpTreeFromMultiVariate(lengthtree,logbgcprocess,1,MEAN,false,meanexp);
		}

		lambdaprocess = 0;
		if (timelambda)	{
			lambdaprocess = new MeanExpTreeFromMultiVariate(lengthtree,logbgcprocess,1,MEAN,false,meanexp);
		}

		/*
		double min = 0.1;
		double max = 10;
		*/
		logbgcoffset = new Normal(Zero,Ten);
		logbgcoffset->setval(0);

		if (Reversible())	{
			lambda = new Jeffreys(min,max,Zero);
			if (inlambda == -1)	{
				lambda->setval(1.0);
			}
			else	{
				lambda->ClampAt(inlambda);
			}

			at2cg = new Jeffreys(min,max,Zero);
			at2gc = new Jeffreys(min,max,Zero);
			at2ta = new Jeffreys(min,max,Zero);
			gc2cg = new Jeffreys(min,max,Zero);

			at2ta->ClampAt(1.0);
			if (inat2cg == -1)	{
				at2cg->setval(1.0);
			}
			else	{
				at2cg->ClampAt(inat2cg);
			}
			if (inat2gc == -1)	{
				at2gc->setval(1.0);
			}
			else	{
				at2gc->ClampAt(inat2gc);
			}
			if (ingc2cg == -1)	{
				gc2cg->setval(1.0);
			}
			else	{
				gc2cg->ClampAt(ingc2cg);
			}
		}
		else	{
			at2cg = new Jeffreys(min,max,Zero);
			at2gc = new Jeffreys(min,max,Zero);
			at2ta = new Jeffreys(min,max,Zero);
			cg2at = new Jeffreys(min,max,Zero);
			cg2gc = new Jeffreys(min,max,Zero);
			cg2ta = new Jeffreys(min,max,Zero);

			at2ta->ClampAt(1.0);
			if (inat2cg == -1)	{
				at2cg->setval(1.0);
			}
			else	{
				at2cg->ClampAt(inat2cg);
			}
			if (inat2gc == -1)	{
				at2gc->setval(1.0);
			}
			else	{
				at2gc->ClampAt(inat2gc);
			}
			if (incg2at == -1)	{
				cg2at->setval(1.0);
			}
			else	{
				cg2at->ClampAt(incg2at);
			}
			if (incg2gc == -1)	{
				cg2gc->setval(1.0);
			}
			else	{
				cg2gc->ClampAt(incg2gc);
			}
			if (incg2ta == -1)	{
				cg2ta->setval(1.0);
			}
			else	{
				cg2ta->ClampAt(incg2ta);
			}
		}

		if (triplet)	{
			if (Reversible())	{
				ctlambda = new Jeffreys(min,max,Zero);
				ctlambda->setval(1.0);

				ctat2cg = new Jeffreys(min,max,Zero);
				ctat2gc = new Jeffreys(min,max,Zero);
				ctat2ta = new Jeffreys(min,max,Zero);
				ctgc2cg = new Jeffreys(min,max,Zero);

				// ctat2ta->ClampAt(1.0);
				ctat2ta->setval(1.0);
				ctat2cg->setval(1.0);
				ctat2gc->setval(1.0);
				ctgc2cg->setval(1.0);
			}
			else	{
				ctat2cg = new Jeffreys(min,max,Zero);
				ctat2gc = new Jeffreys(min,max,Zero);
				ctat2ta = new Jeffreys(min,max,Zero);
				ctcg2at = new Jeffreys(min,max,Zero);
				ctcg2gc = new Jeffreys(min,max,Zero);
				ctcg2ta = new Jeffreys(min,max,Zero);

				// ctat2ta->ClampAt(1.0);
				ctat2ta->setval(1.0);
				ctat2cg->setval(1.0);
				ctat2gc->setval(1.0);
				ctcg2at->setval(1.0);
				ctcg2gc->setval(1.0);
				ctcg2ta->setval(1.0);
			}
			cpgrate = new Jeffreys(min,max,Zero);
			cpgrate->setval(2.0);
		}

		alpha1 = new Jeffreys(min,max,Zero);
		beta1 = new Jeffreys(min,max,Zero);
		alpha2 = new Jeffreys(min,max,Zero);
		// beta2 = new Jeffreys(min,max,Zero);
		alpha1->setval(1);
		beta1->setval(1);
		alpha2->setval(1);
		// beta2->setval(1);
		rootbgc = new Gamma*[Ngene];
		rootbgcoffset = new Gamma*[Ngene];
		// rootbgc = new Jeffreys*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			// rootbgc[gene] = new Gamma(One,Tenth);
			// rootbgc[gene] = new Jeffreys(min,max,Zero);
			rootbgc[gene] = new Gamma(alpha1,beta1);
			rootbgcoffset[gene] = new Gamma(alpha2,alpha2);
			rootbgc[gene]->setval(1.0);
			if (clampbgcoffset)	{
				rootbgcoffset[gene]->ClampAt(1.0);
			}
		}

		mutmatrix = 0;
		nonrevmutmatrix = 0;
		ctmutmatrix = 0;
		ctnonrevmutmatrix= 0;
		mutmatrixtree = 0;
		if (timelambda)	{
			if (! Reversible())	{
				cerr << "error : time-dependent mutation process not yet implemented in the non reversible case\n";
				exit(1);
			}
			else	{
				mutmatrixtree = new MutMatrixTree(lambdaprocess,lambda,at2cg,at2gc,at2ta,gc2cg,true);
			}
		}
		else	{
			cerr << "mutmatrix\n";
			if (Reversible())	{
				mutmatrix = new RandomMutSubMatrix(lambda,at2cg,at2gc,at2ta,gc2cg,true);
			}
			else	{
				nonrevmutmatrix = new RandomNonRevMutSubMatrix(at2cg,at2gc,at2ta,cg2at,cg2gc,cg2ta,false,discn);
			}
			if (triplet)	{
				if (Reversible())	{
					ctmutmatrix = new RandomMutSubMatrix(ctlambda,ctat2cg,ctat2gc,ctat2ta,ctgc2cg,true);
				}
				else	{
					ctnonrevmutmatrix = new RandomNonRevMutSubMatrix(ctat2cg,ctat2gc,ctat2ta,ctcg2at,ctcg2gc,ctcg2ta,false,discn);
				}
			}
		}

		cerr << "bggcprocess\n";
		bgcprocess = new BGCProcess(logbgcprocess,logbgcoffset);

		// genealpha = new Gamma(One,Tenth);
		// rho = new Gamma(One,Tenth);
		genealpha = new Jeffreys(min,max,Zero);
		genealpha2 = new Jeffreys(min,max,Zero);
		rho = 0;
		rhosigma = 0;
		rhoprocess = 0;
		rho = new Jeffreys(min,max,Zero);
		genealpha->setval(1);
		genealpha2->setval(1);
		rho->setval(0.1);
		if (flexrho)	{
			rhosigma = new Jeffreys(mindiag,maxdiag,Zero);
			rhosigma->setval(1);
			rhoprocess = new LogNormalTreeProcess(lengthtree,rhosigma,INTEGRAL);
			rhoprocess->ClampRootAt(0);
			// rho->ClampAt(1);
		}
		argeneprocess = 0;
		lngeneprocess = 0;
		gammageneprocess = 0;
		wngeneprocess = 0;
		mixedgeneprocess = 0;
		geneprocess = new LengthTree*[Ngene];
		if (geneprocesstype == 1)	{
			genealpha->setval(10);
			gammageneprocess = new GammaTree*[Ngene];
		}
		else if (geneprocesstype == 2)	{
			genealpha->setval(10);
			wngeneprocess = new GammaWhiteNoiseProcess*[Ngene];
		}
		else if (geneprocesstype == 3)	{
			genealpha->setval(1);
			lngeneprocess = new LogNormalTreeProcess*[Ngene];
		}
		else if (geneprocesstype == 4)	{
			wngeneprocess = new GammaWhiteNoiseProcess*[Ngene];
			argeneprocess = new AutoRegressiveLogNormalTreeProcess*[Ngene];
			mixedgeneprocess = new BranchProductProcess*[Ngene];
		}
		else	{
			argeneprocess = new AutoRegressiveLogNormalTreeProcess*[Ngene];
		}
		branchproductprocess = new BranchProductProcess*[Ngene];

		nucmatrixtree = new NucMatrixTree*[Ngene];

		for (int gene=0; gene<Ngene; gene++)	{
			cerr << gene << '\n';
			if (geneprocesstype == 1)	{
				gammageneprocess[gene] = new GammaTree(GetFineGrainedTree(),genealpha,genealpha,false);
				// gammageneprocess[gene]->Reset();
				geneprocess[gene] = gammageneprocess[gene];
			}
			else if (geneprocesstype == 2)	{
				genealpha->setval(100);
				wngeneprocess[gene] = new GammaWhiteNoiseProcess(lengthtree,genealpha);
				wngeneprocess[gene]->Reset();
				geneprocess[gene] = wngeneprocess[gene];
			}
			else if (geneprocesstype == 3)	{
				lngeneprocess[gene] = new LogNormalTreeProcess(lengthtree,genealpha,MEAN);
				lngeneprocess[gene]->ClampRootAt(0);
				lngeneprocess[gene]->Reset();
				geneprocess[gene] = lngeneprocess[gene];
			}
			else if (geneprocesstype == 4)	{
				argeneprocess[gene] = new AutoRegressiveLogNormalTreeProcess(lengthtree,genealpha,rho,MEAN,false);
				argeneprocess[gene]->Reset();
				genealpha2->setval(100);
				wngeneprocess[gene] = new GammaWhiteNoiseProcess(lengthtree,genealpha2);
				wngeneprocess[gene]->Reset();
				mixedgeneprocess[gene] = new BranchProductProcess(argeneprocess[gene],wngeneprocess[gene],0,clampbgc);
				geneprocess[gene] = mixedgeneprocess[gene];
			}
			else if (flexrho)	{
				argeneprocess[gene] = new AutoRegressiveLogNormalTreeProcess(rhoprocess,genealpha,rho,MEAN,false);
				argeneprocess[gene]->Reset();
				geneprocess[gene] = argeneprocess[gene];
				// geneprocess[gene] = new FlexRhoAutoRegressiveLogNormalTreeProcess(lengthtree,genealpha,rhoprocess,rho,MEAN,false);
			}
			else	{
				argeneprocess[gene] = new AutoRegressiveLogNormalTreeProcess(lengthtree,genealpha,rho,MEAN,false);
				argeneprocess[gene]->Reset();
				geneprocess[gene] = argeneprocess[gene];
			}
			cerr << "branchproduct\n";
			branchproductprocess[gene] = new BranchProductProcess(bgcprocess,geneprocess[gene],rootbgcoffset[gene],clampbgc);
			cerr << "nucmatrix\n";
			if (Reversible())	{
				if (triplet)	{
					nucmatrixtree[gene] = new BGCCpGMatrixTree(branchproductprocess[gene],mutmatrix,ctmutmatrix,cpgrate,alpha,rootbgc[gene],normalise,discgam);
				}
				else	{
					if (alphaprocess)	{
						if (withhotspot)	{
							nucmatrixtree[gene] = new HotSpotBGCMatrixTree(branchproductprocess[gene],mutmatrix,alphaprocess,alpha,rootbgc[gene],normalise);
						}
						else	{
							nucmatrixtree[gene] = new AlphaBGCMatrixTree(branchproductprocess[gene],mutmatrix,alphaprocess,alpha,rootbgc[gene],normalise);
						}
					}
					else	{
						nucmatrixtree[gene] = new BGCMatrixTree(branchproductprocess[gene],mutmatrix,alpha,rootbgc[gene],normalise,discgam,mutmatrixtree);
					}
				}
			}
			else	{
				if (triplet)	{
					nucmatrixtree[gene] = new NonRevBGCCpGMatrixTree(synratetree,branchproductprocess[gene],nonrevmutmatrix,ctnonrevmutmatrix,cpgrate,alpha,rootbgc[gene],normalise,discgam,discn);
				}
				else	{
					nucmatrixtree[gene] = new NonRevBGCMatrixTree(synratetree,branchproductprocess[gene],nonrevmutmatrix,alpha,rootbgc[gene],normalise,discgam,discn);
				}
			}
			cerr << "ok\n";
		}

		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i,i,contdatatype);
			}
		}

		// make substitution mappings
		if (conjpath)	{
			pathconjtree = new BranchMatrixPathConjugateTree*[Ngene];
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				pathconjtree[gene] = new BranchMatrixPathConjugateTree(synratetree, nucmatrixtree[gene], GetData(gene));
				phyloprocess[gene] = new PathConjugatePhyloProcess(pathconjtree[gene]);
			}
		}
		else	{
			pathconjtree = 0;
			if (priorsampling)	{
				phyloprocess = 0;
			}
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene] = new BranchMatrixPhyloProcess(synratetree, nucmatrixtree[gene], GetData(gene));
			}
		}

		if (phyloprocess)	{
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene]->Unfold();
			}
		}
		if (sample)	{
			if (Split())	{
				cerr << "before\n";
				SplitLengthTree* tmp = dynamic_cast<SplitLengthTree*> (lengthtree);
				tmp->Check();
			}

			if (phyloprocess)	{
				for (int gene=0; gene<Ngene; gene++)	{
					phyloprocess[gene]->Sample();
				}
			}
		}

		// register model
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(Ten);
		RootRegister(Tenth);
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (rootmean)	{
			RootRegister(rootmean);
			RootRegister(rootvar);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
			if (isCalibrated())	{
				cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
				cerr << "scale progeny : " << GetCalibratedChronogram()->GetScale()->down.size() << '\n';
			}
		}

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BGCModel() {}

	Tree* GetTree() {return tree;}
	Tree* GetFineGrainedTree() {return splittree;}

	// BranchVarTree<UnitReal>* GetGCTree() {return gcprocess;}

	MultiVariateTreeProcess* GetContProcess() {return contprocess;}
	NodeVarTree<RealVector>* GetLogBGCProcess() {return logbgcprocess;}
	Rvar<Real>* GetLogBGCOffset() {return logbgcoffset;}

	int GetNcont()	{
		return Ncont;
	}

	int GetNreg()	{
		return Nreg;
	}

	BGCProcess* GetBGCProcess()	{
		return bgcprocess;
	}

	BGCMatrixTree* GetBGCMatrixTree(int gene)	{
		if (triplet)	{
			cerr << "error in get bgc matrix tree\n";
			exit(1);
		}
		BGCMatrixTree* tmp = dynamic_cast<BGCMatrixTree*>(nucmatrixtree[gene]);
		if (! tmp)	{
			cerr << "error in get bgc matrix tree\n";
			exit(1);
		}
		return tmp;
	}

	BGCCpGMatrixTree* GetBGCCpGMatrixTree(int gene)	{
		if (!triplet)	{
			cerr << "error in get bgc matrix tree\n";
			exit(1);
		}
		BGCCpGMatrixTree* tmp = dynamic_cast<BGCCpGMatrixTree*>(nucmatrixtree[gene]);
		if (! tmp)	{
			cerr << "error in get bgc matrix tree\n";
			exit(1);
		}
		return tmp;
	}

	Chronogram* GetChronogram() {
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

	ContinuousData* GetContinuousData() {return contdata;}

	CovMatrix* GetContMatrix() {return contsigma;}

	CalibratedChronogram* GetCalibratedChronogram()	{
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

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();

		total += synsigma->GetLogProb();
		total += synratetree->GetLogProb();

		total += ContDiagArray->GetLogProb();
		total += BGCDiagArray->GetLogProb();
		total += contsigma->GetLogProb();
		total += bgcsigma->GetLogProb();
		total += regarray->GetLogProb();
		if (autoregressive)	{
			total += phi->GetLogProb();
		}
		if (alpha)	{
			total += alpha->GetLogProb();
		}
		total += contprocess->GetLogProb();
		total += logbgcprocess->GetLogProb();
		if (Reversible())	{
			total += lambda->GetLogProb();
			total += at2cg->GetLogProb();
			total += at2gc->GetLogProb();
			total += at2ta->GetLogProb();
			total += gc2cg->GetLogProb();
			if (triplet)	{
				total += ctlambda->GetLogProb();
				total += ctat2cg->GetLogProb();
				total += ctat2gc->GetLogProb();
				total += ctat2ta->GetLogProb();
				total += ctgc2cg->GetLogProb();
			}
		}
		else	{
			total += at2cg->GetLogProb();
			total += at2gc->GetLogProb();
			total += at2ta->GetLogProb();
			total += cg2at->GetLogProb();
			total += cg2gc->GetLogProb();
			total += cg2ta->GetLogProb();
			if (triplet)	{
				total += ctat2cg->GetLogProb();
				total += ctat2gc->GetLogProb();
				total += ctat2ta->GetLogProb();
				total += ctcg2at->GetLogProb();
				total += ctcg2gc->GetLogProb();
				total += ctcg2ta->GetLogProb();
			}
		}
		if (triplet)	{
			total += cpgrate->GetLogProb();
		}

		total += logbgcoffset->GetLogProb();

		total += rho->GetLogProb();
		if (flexrho)	{
			rhosigma->GetLogProb();
			rhoprocess->GetLogProb();
		}
		total += genealpha->GetLogProb();
		total += alpha1->GetLogProb();
		total += beta1->GetLogProb();
		total += alpha2->GetLogProb();
		// total += beta2->GetLogProb();
		for (int gene=0; gene<Ngene; gene++)	{
			total += rootbgcoffset[gene]->GetLogProb();
			total += rootbgc[gene]->GetLogProb();
			if (geneprocesstype == 1)	{
				total += gammageneprocess[gene]->GetLogProb();
			}
			else if (geneprocesstype == 2)	{
				total += wngeneprocess[gene]->GetLogProb();
			}
			else if (geneprocesstype == 3)	{
				total += lngeneprocess[gene]->GetLogProb();
			}
			else if (geneprocesstype == 4)	{
				total += genealpha2->GetLogProb();
				total += wngeneprocess[gene]->GetLogProb();
				total += argeneprocess[gene]->GetLogProb();
			}
			else	{
				total += argeneprocess[gene]->GetLogProb();
			}
		}
		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (priorsampling)	{
			return 0;
		}
		else if (clampsuffstat)	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += pathconjtree[gene]->GetLogProb();
			}
		}
		else	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += phyloprocess[gene]->GetLogProb();
			}
		}
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			if (! clampsuffstat)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new DSemiConjugateMappingMove(phyloprocess[gene],pathconjtree[gene]),1,"mapping + sufficient stat");
				}
			}
		}
		else	{
			if (phyloprocess)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new SimpleMove(phyloprocess[gene],1),1,"mapping");
				}
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
			if (! clamptree)	{
				if ((chronoprior >= 1) && (chronoprior <= 3))	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
			}

			if (isCalibrated() && (! clamptree))	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");
			}

			scheduler.Register(new SimpleMove(synratetree,10),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.01),10,"synratetree");

			scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");

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

			scheduler.Register(new SimpleMove(contprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.01),10,"multinormal");

			// scheduler.Register(new SimpleMove(logbgcprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(logbgcprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(logbgcprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(logbgcprocess,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(contsigma,10),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.01),100,"cont sigma");

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

			scheduler.Register(new SimpleMove(bgcsigma,10),100,"sub sigma");
			scheduler.Register(new SimpleMove(bgcsigma,1),100,"sub sigma");
			scheduler.Register(new SimpleMove(bgcsigma,0.1),100,"sub sigma");
			scheduler.Register(new SimpleMove(bgcsigma,0.01),100,"sub sigma");

			scheduler.Register(new SimpleMove(BGCDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(BGCDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(BGCDiagArray,0.1),10,"theta");

			for (int gene=0; gene<Ngene; gene++)	{
				scheduler.Register(new SimpleMove(rootbgcoffset[gene],1),10,"root bgc");
				scheduler.Register(new SimpleMove(rootbgcoffset[gene],0.1),10,"root bgc");
				scheduler.Register(new SimpleMove(rootbgcoffset[gene],0.01),10,"root bgc");
				scheduler.Register(new SimpleMove(rootbgc[gene],1),10,"root bgc");
				scheduler.Register(new SimpleMove(rootbgc[gene],0.1),10,"root bgc");
				scheduler.Register(new SimpleMove(rootbgc[gene],0.01),10,"root bgc");

				if (geneprocesstype == 1)	{
					scheduler.Register(new SimpleMove(gammageneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(gammageneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(gammageneprocess[gene],0.01),10,"geneprocess");
					if (clampbgcoffset)	{
						scheduler.Register(new MultiplicativeCompensatoryMove(gammageneprocess[gene],0,1),10,"geneprocess mult");
						scheduler.Register(new MultiplicativeCompensatoryMove(gammageneprocess[gene],0,0.1),10,"geneprocess mult");
						scheduler.Register(new MultiplicativeCompensatoryMove(gammageneprocess[gene],0,0.01),10,"geneprocess mult");
					}
				}
				else if (geneprocesstype == 2)	{
					scheduler.Register(new SimpleMove(wngeneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(wngeneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(wngeneprocess[gene],0.01),10,"geneprocess");
				}
				else if (geneprocesstype == 3)	{
					scheduler.Register(new SimpleMove(lngeneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(lngeneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(lngeneprocess[gene],0.01),10,"geneprocess");
				}
				else if (geneprocesstype == 4)	{
					scheduler.Register(new SimpleMove(wngeneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(wngeneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(wngeneprocess[gene],0.01),10,"geneprocess");
					scheduler.Register(new SimpleMove(argeneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(argeneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(argeneprocess[gene],0.01),10,"geneprocess");
				}
				else	{
					scheduler.Register(new SimpleMove(argeneprocess[gene],1),10,"geneprocess");
					scheduler.Register(new SimpleMove(argeneprocess[gene],0.1),10,"geneprocess");
					scheduler.Register(new SimpleMove(argeneprocess[gene],0.01),10,"geneprocess");
				}
			}

			scheduler.Register(new SimpleMove(alpha1,10),10,"alpha1");
			scheduler.Register(new SimpleMove(alpha1,1),10,"alpha1");
			scheduler.Register(new SimpleMove(alpha1,0.1),10,"alpha1");
			scheduler.Register(new SimpleMove(alpha1,0.01),10,"alpha1");

			scheduler.Register(new SimpleMove(beta1,10),10,"beta1");
			scheduler.Register(new SimpleMove(beta1,1),10,"beta1");
			scheduler.Register(new SimpleMove(beta1,0.1),10,"beta1");
			scheduler.Register(new SimpleMove(beta1,0.01),10,"beta1");

			scheduler.Register(new SimpleMove(alpha2,10),10,"alpha2");
			scheduler.Register(new SimpleMove(alpha2,1),10,"alpha2");
			scheduler.Register(new SimpleMove(alpha2,0.1),10,"alpha2");
			scheduler.Register(new SimpleMove(alpha2,0.01),10,"alpha2");

			/*
			scheduler.Register(new SimpleMove(beta2,10),10,"beta2");
			scheduler.Register(new SimpleMove(beta2,1),10,"beta2");
			scheduler.Register(new SimpleMove(beta2,0.1),10,"beta2");
			scheduler.Register(new SimpleMove(beta2,0.01),10,"beta2");
			*/

			if (geneprocesstype == 4)	{
				scheduler.Register(new SimpleMove(genealpha2,1),10,"genealpha2");
				scheduler.Register(new SimpleMove(genealpha2,0.1),10,"genealpha2");
				scheduler.Register(new SimpleMove(genealpha2,0.01),10,"genealpha2");
			}
			scheduler.Register(new SimpleMove(genealpha,1),10,"genealpha");
			scheduler.Register(new SimpleMove(genealpha,0.1),10,"genealpha");
			scheduler.Register(new SimpleMove(genealpha,0.01),10,"genealpha");

			scheduler.Register(new SimpleMove(rho,1),10,"rho");
			scheduler.Register(new SimpleMove(rho,0.1),10,"rho");
			scheduler.Register(new SimpleMove(rho,0.01),10,"rho");

			if (flexrho)	{
				scheduler.Register(new SimpleMove(rhosigma,1),10,"rho sigma");
				scheduler.Register(new SimpleMove(rhosigma,0.1),10,"rho sigma");
				scheduler.Register(new SimpleMove(rhosigma,0.01),10,"rho sigma");

				scheduler.Register(new SimpleMove(rhoprocess,1),10,"rho process");
				scheduler.Register(new SimpleMove(rhoprocess,0.1),10,"rho process");
				scheduler.Register(new SimpleMove(rhoprocess,0.01),10,"rho process");
			}

			if (alpha)	{
				scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
				scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
				scheduler.Register(new SimpleMove(alpha,0.01),10,"alpha");
			}

			if (Reversible())	{
				scheduler.Register(new SimpleMove(at2cg,1),10,"at2cg");
				scheduler.Register(new SimpleMove(at2cg,0.1),10,"at2cg");
				scheduler.Register(new SimpleMove(at2cg,0.01),10,"at2cg");

				scheduler.Register(new SimpleMove(at2gc,1),10,"at2gc");
				scheduler.Register(new SimpleMove(at2gc,0.1),10,"at2gc");
				scheduler.Register(new SimpleMove(at2gc,0.01),10,"at2gc");

				scheduler.Register(new SimpleMove(at2ta,1),10,"at2ta");
				scheduler.Register(new SimpleMove(at2ta,0.1),10,"at2ta");
				scheduler.Register(new SimpleMove(at2ta,0.01),10,"at2ta");

				scheduler.Register(new SimpleMove(gc2cg,1),10,"gc2cg");
				scheduler.Register(new SimpleMove(gc2cg,0.1),10,"gc2cg");
				scheduler.Register(new SimpleMove(gc2cg,0.01),10,"gc2cg");

				scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
				scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
				scheduler.Register(new SimpleMove(lambda,0.01),10,"lambda");

				if (triplet)	{
					scheduler.Register(new SimpleMove(ctat2cg,1),10,"at2cg");
					scheduler.Register(new SimpleMove(ctat2cg,0.1),10,"at2cg");
					scheduler.Register(new SimpleMove(ctat2cg,0.01),10,"at2cg");

					scheduler.Register(new SimpleMove(ctat2gc,1),10,"at2gc");
					scheduler.Register(new SimpleMove(ctat2gc,0.1),10,"at2gc");
					scheduler.Register(new SimpleMove(ctat2gc,0.01),10,"at2gc");

					scheduler.Register(new SimpleMove(ctat2ta,1),10,"at2ta");
					scheduler.Register(new SimpleMove(ctat2ta,0.1),10,"at2ta");
					scheduler.Register(new SimpleMove(ctat2ta,0.01),10,"at2ta");

					scheduler.Register(new SimpleMove(ctgc2cg,1),10,"gc2cg");
					scheduler.Register(new SimpleMove(ctgc2cg,0.1),10,"gc2cg");
					scheduler.Register(new SimpleMove(ctgc2cg,0.01),10,"gc2cg");

					scheduler.Register(new SimpleMove(ctlambda,1),10,"lambda");
					scheduler.Register(new SimpleMove(ctlambda,0.1),10,"lambda");
					scheduler.Register(new SimpleMove(ctlambda,0.01),10,"lambda");
				}

			}
			else	{
				scheduler.Register(new SimpleMove(at2cg,1),10,"at2cg");
				scheduler.Register(new SimpleMove(at2cg,0.1),10,"at2cg");
				scheduler.Register(new SimpleMove(at2cg,0.01),10,"at2cg");

				scheduler.Register(new SimpleMove(at2gc,1),10,"at2gc");
				scheduler.Register(new SimpleMove(at2gc,0.1),10,"at2gc");
				scheduler.Register(new SimpleMove(at2gc,0.01),10,"at2gc");

				scheduler.Register(new SimpleMove(at2ta,1),10,"at2ta");
				scheduler.Register(new SimpleMove(at2ta,0.1),10,"at2ta");
				scheduler.Register(new SimpleMove(at2ta,0.01),10,"at2ta");

				scheduler.Register(new SimpleMove(cg2at,1),10,"cg2at");
				scheduler.Register(new SimpleMove(cg2at,0.1),10,"cg2cat");
				scheduler.Register(new SimpleMove(cg2at,0.01),10,"cg2at");

				scheduler.Register(new SimpleMove(cg2gc,1),10,"cg2gc");
				scheduler.Register(new SimpleMove(cg2gc,0.1),10,"cg2gc");
				scheduler.Register(new SimpleMove(cg2gc,0.01),10,"cg2gc");

				scheduler.Register(new SimpleMove(cg2ta,1),10,"cg2ta");
				scheduler.Register(new SimpleMove(cg2ta,0.1),10,"cg2cta");
				scheduler.Register(new SimpleMove(cg2ta,0.01),10,"cg2ta");

				if (triplet)	{
					scheduler.Register(new SimpleMove(ctat2cg,1),10,"at2cg");
					scheduler.Register(new SimpleMove(ctat2cg,0.1),10,"at2cg");
					scheduler.Register(new SimpleMove(ctat2cg,0.01),10,"at2cg");

					scheduler.Register(new SimpleMove(ctat2gc,1),10,"at2gc");
					scheduler.Register(new SimpleMove(ctat2gc,0.1),10,"at2gc");
					scheduler.Register(new SimpleMove(ctat2gc,0.01),10,"at2gc");

					scheduler.Register(new SimpleMove(ctat2ta,1),10,"at2ta");
					scheduler.Register(new SimpleMove(ctat2ta,0.1),10,"at2ta");
					scheduler.Register(new SimpleMove(ctat2ta,0.01),10,"at2ta");

					scheduler.Register(new SimpleMove(ctcg2at,1),10,"cg2at");
					scheduler.Register(new SimpleMove(ctcg2at,0.1),10,"cg2cat");
					scheduler.Register(new SimpleMove(ctcg2at,0.01),10,"cg2at");

					scheduler.Register(new SimpleMove(ctcg2gc,1),10,"cg2gc");
					scheduler.Register(new SimpleMove(ctcg2gc,0.1),10,"cg2gc");
					scheduler.Register(new SimpleMove(ctcg2gc,0.01),10,"cg2gc");

					scheduler.Register(new SimpleMove(ctcg2ta,1),10,"cg2ta");
					scheduler.Register(new SimpleMove(ctcg2ta,0.1),10,"cg2cta");
					scheduler.Register(new SimpleMove(ctcg2ta,0.01),10,"cg2ta");
				}
			}
			if (triplet)	{
				scheduler.Register(new SimpleMove(cpgrate,1),10,"cpgrate");
				scheduler.Register(new SimpleMove(cpgrate,0.1),10,"cpgrate");
				scheduler.Register(new SimpleMove(cpgrate,0.01),10,"cpgrate");
			}

			scheduler.Register(new SimpleMove(logbgcoffset,1),10,"bgc offset");
			scheduler.Register(new SimpleMove(logbgcoffset,0.1),10,"bgc offset");
			scheduler.Register(new SimpleMove(logbgcoffset,0.01),10,"bgc offset");

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

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			Chi->Sample();
			Chi2->Sample();
		}
		chronogram->Sample();

		synsigma->Sample();
		synratetree->Sample();

		ContDiagArray->Sample();
		BGCDiagArray->Sample();

		contsigma->Sample();
		bgcsigma->Sample();
		regarray->Sample();
		phi->Sample();
		// sigma->SetIdentity();

		contprocess->Sample();
		logbgcprocess->Sample();

		logbgcprocess->CutOff(1,0);
		if (alpha)	{
			alpha->Sample();
		}
		if (Reversible())	{
			lambda->Sample();
			at2cg->Sample();
			at2gc->Sample();
			at2ta->Sample();
			gc2cg->Sample();
			if (triplet)	{
				ctlambda->Sample();
				ctat2cg->Sample();
				ctat2gc->Sample();
				ctat2ta->Sample();
				ctgc2cg->Sample();
			}
		}
		else	{
			at2cg->Sample();
			at2gc->Sample();
			at2ta->Sample();
			cg2at->Sample();
			cg2gc->Sample();
			cg2ta->Sample();
			if (triplet)	{
				ctat2cg->Sample();
				ctat2gc->Sample();
				ctat2ta->Sample();
				ctcg2at->Sample();
				ctcg2gc->Sample();
				ctcg2ta->Sample();
			}
		}
		if (triplet)	{
			cpgrate->Sample();
		}
		logbgcoffset->Sample();

		// gcprocess->specialUpdate();

		alpha1->Sample();
		beta1->Sample();
		alpha2->Sample();
		// beta2->Sample();

		for (int gene=0; gene<Ngene; gene++)	{
			rootbgcoffset[gene]->Sample();
			rootbgc[gene]->Sample();

			if (geneprocesstype == 1)	{
				gammageneprocess[gene]->Sample();
			}
			else if (geneprocesstype == 2)	{
				wngeneprocess[gene]->Sample();
			}
			else if (geneprocesstype == 3)	{
				lngeneprocess[gene]->Sample();
			}
			else if (geneprocesstype == 4)	{
				genealpha2->Sample();
				argeneprocess[gene]->Sample();
				wngeneprocess[gene]->Sample();
			}
			else	{
				argeneprocess[gene]->Sample();
			}
		}

		// gctree->specialUpdate()

		if (phyloprocess)	{
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene]->Sample();
			}
		}

		cerr << "ok\n";
	}

	Var<PosReal>* GetScale()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	int GetNgene()	{
		return Ngene;
	}

	bool FlexRho()	{
		return flexrho;
	}

	double GetRho()	{
		return rho->val();
	}

	double GetLambda()	{
		if (!Reversible())	{
			double at = cg2at->val() + cg2ta->val();
			double cg = at2cg->val() + at2gc->val();
			return at / cg;
		}
		return lambda->val();
	}

	double GetContextLambda()	{
		return ctlambda->val();
	}

	double GetCpGRate()	{
		if (! triplet)	{
			cerr << "error : cpg not defined in site-independent model\n";
			exit(1);
		}
		return cpgrate->val();
	}

	double GetAlpha()	{
		if (! alpha)	{
			return 0;
		}
		return alpha->val();
	}

	double GetStatVar()	{
		return genealpha->val() / 2 / rho->val();
	}

	double GetVarCoeff()	{
		return 1.0 / sqrt(genealpha->val());
	}

	LogNormalTreeProcess* GetRhoProcess()	{
		return rhoprocess;
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	LogNormalTreeProcess* GetSynRateTree()	{
		return synratetree;
	}

	LengthTree* GetGeneProcess(int gene)	{
		return geneprocess[gene];
	}

	void UpdateGeneProcess(int gene)	{
		if (geneprocesstype == 4)	{
			argeneprocess[gene]->specialUpdate();
			mixedgeneprocess[gene]->specialUpdate();
		}
		else if (! geneprocesstype || (geneprocesstype == 3))	{
			argeneprocess[gene]->specialUpdate();
		}
	}

	void UpdateGeneProcesses()	{
		for (int gene=0; gene<Ngene; gene++)	{
			UpdateGeneProcess(gene);
		}
	}

	void UpdateBranchProductProcesses()	{
		for (int gene=0; gene<Ngene; gene++)	{
			branchproductprocess[gene]->specialUpdate();
		}
	}

	void UpdateBGC()	{
		GetBGCProcess()->specialUpdate();
		UpdateGeneProcesses();
		UpdateBranchProductProcesses();
	}

	double GetMeanSynRate()	{
		return GetSynRateTree()->GetMeanRate();
	}

	double GetMeanLambda()	{
		return lambdaprocess->GetMean();
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	int GetTotalBGCOverflowCount()	{
		int total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (triplet)	{
				total += GetBGCCpGMatrixTree(gene)->GetBGCOverflowCount();
			}
			else	{
				total += GetBGCMatrixTree(gene)->GetBGCOverflowCount();
			}
		}
		return total;
	}

	double GetPruningTime()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += phyloprocess[gene]->GetPruningTime();
		}
		return total;
	}

	double GetResampleTime()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += phyloprocess[gene]->GetResampleTime();
		}
		return total;
	}

	double GetGeneVarCoeff(int gene)	{
		double mean = 0;
		double var = 0;
		geneprocess[gene]->GetMeanAndVar(mean,var);
		return var / mean / mean;
	}

	double GetMeanFractionAbove(double cutoff, double& meanjoint, double& mean2, double& mean1, double& mean0)	{
		double meann = 0;
		meanjoint = 0;
		mean2 = 0;
		mean1 = 0;
		mean0 = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			int n = 0;
			int totjoint = 0;
			geneprocess[gene]->GetFractionAbove(cutoff,n,totjoint);

			meann += n;
			if (n == 2)	{
				meanjoint += totjoint;
				mean2++;
			}
			if (n == 1)	{
				mean1++;
			}
			if (n == 0)	{
				mean0++;
			}
		}
		meann /= Ngene;
		meanjoint /= Ngene;
		mean2 /= Ngene;
		mean1 /= Ngene;
		mean0 /= Ngene;
		return  meann;
	}

	double GetGeneVarCoeff()	{
		double mean = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			mean += GetGeneVarCoeff(gene);
		}
		mean /= Ngene;
		return mean;
	}

	double GetPropGreaterThan(const Link* from, double threshold)	{
		double n = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			double tmp = branchproductprocess[gene]->GetBranchVal(from->GetBranch())->val();
			if (tmp > threshold)	{
				n++;
			}
		}
		n /= Ngene;
		return n;
	}

	double GetBranchVarCoeff(const Link* from)	{
		double mean = 0;
		double var = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			double tmp = branchproductprocess[gene]->GetBranchVal(from->GetBranch())->val();
			// double tmp = geneprocess[gene]->GetBranchVal(from->GetBranch())->val() * rootbgcoffset[gene]->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= Ngene;
		var /= Ngene;
		var -= mean*mean;
		return sqrt(var / mean / mean);
	}

	double GetSuccessiveBranchesCorrelation()	{
		int n = 0;
		double total = RecursiveGetSuccessiveBranchesCorrelation(GetTree()->GetRoot(),GetTree()->GetRoot(),n);
		return total / n;
	}

	double RecursiveGetSuccessiveBranchesCorrelation(const Link* fromup, const Link* from, int& n)	{

		double total = 0;
		if (!fromup->isRoot())	{
			double time = 0;
			total += GetBranchCorrelation(fromup,from,time);
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetSuccessiveBranchesCorrelation(from,link->Out(),n);
		}
		return total;
	}

	double GetBranchCorrelation(const Link* from1, const Link* from2, double& time)	{
		double mean1 = 0;
		double mean2 = 0;
		double var1 = 0;
		double var2 = 0;
		double covar = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			/*
			double tmp1 = geneprocess[gene]->GetBranchVal(from1->GetBranch())->val();
			double tmp2 = geneprocess[gene]->GetBranchVal(from2->GetBranch())->val();
			*/
			double tmp1 = branchproductprocess[gene]->GetBranchVal(from1->GetBranch())->val();
			double tmp2 = branchproductprocess[gene]->GetBranchVal(from2->GetBranch())->val();
			mean1 += tmp1;
			var1 += tmp1 * tmp1;
			mean2 += tmp2;
			var2 += tmp2 * tmp2;
			covar += tmp1 * tmp2;
		}
		mean1 /= Ngene;
		var1 /= Ngene;
		var1 -= mean1*mean1;
		mean2 /= Ngene;
		var2 /= Ngene;
		var2 -= mean2*mean2;
		covar /= Ngene;
		covar -= mean1*mean2;

		time = chronogram->GetDistance(from1,from2);
		if (time < 0)	{
			cerr << "error: negative time: " << GetTree()->GetLeftMost(from1) << '\t' << GetTree()->GetRightMost(from1) << '\t' << GetTree()->GetLeftMost(from2) << '\t' << GetTree()->GetRightMost(from2) << '\n';
			exit(1);
		}
		// time = fabs(chronogram->GetMidTime(from1) - chronogram->GetMidTime(from2));
		return covar / sqrt(var1 * var2);
	}

	double GetBranchVarCoeff()	{
		int count = 0;
		double totlength = 0;
		double tot = RecursiveGetBranchVarCoeff(lengthtree->GetRoot(),count,totlength);
		// return tot / totlength;
		return tot / count;
	}

	double RecursiveGetBranchVarCoeff(const Link* from, int& count, double& length)	{
		double tot = 0;
		if (! from->isRoot())	{
			count++;
			double tmp = GetLengthTree()->GetBranchVal(from->GetBranch())->val();
			length += tmp;
			tot += GetBranchVarCoeff(from);
			// tot += tmp * GetBranchVarCoeff(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveGetBranchVarCoeff(link->Out(),count,length);
		}
		return tot;
	}

	/*
	double GetMeanGC()	{
		return gctree->GetMeanGCContent();
	}

	double GetVarGC()	{
		return gctree->GetVarGCContent();
	}
	*/

	double GetBGCVar(int i=0)	{
		return (*bgcsigma)[i][i];
	}

	double GetRootBGC(int gene)	{
		return rootbgc[gene]->val();
	}

	double GetMeanRootBGC()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += rootbgc[gene]->val();
		}
		return total / Ngene;
	}

	double GetRootBGCOffset(int gene)	{
		return rootbgcoffset[gene]->val();
	}

	double GetMeanRootBGCOffset()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += rootbgcoffset[gene]->val();
		}
		return total / Ngene;
	}

	double GetRegCoef(int i=0)	{
		return (*regarray->GetVal(0))[i];
	}

	double GetAlphaRegCoef(int i=0)	{
		return (*regarray->GetVal(1))[i];
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tsynrate";

		os << "\tsynsigma";
		if (Reversible())	{
			os << "\tlambda\tat2cg\tat2gc\tat2ta\tgc2cg";
		}
		else	{
			os << "\tat2cg\tat2gc\tat2ta\tcg2at\tcg2gc\tcg2ta";
		}
		if (triplet)	{
			os << "\tcpgrate";
		}

		if (timelambda)	{
			os << "\tmeanlambda";
		}

		os << "\talpha\tbgcoffset";

		for (int i=0; i<Ncont; i++)	{
			os << '\t' << "reg" << i;
		}
		if (withalphaprocess)	{
			for (int i=0; i<Ncont; i++)	{
				os << '\t' << "areg" << i;
			}
		}
		if (timelambda)	{
			for (int i=0; i<Ncont; i++)	{
				os << '\t' << "lreg" << i;
			}
		}


		if (! geneprocesstype)	{
			os << '\t' << "rho";
			if (flexrho)	{
				os << '\t' << "rhomean";
				os << '\t' << "rhovar";
				os << '\t' << "rhosigma";
			}
			os << '\t' << "statvar";
		}
		else if (geneprocesstype == 4)	{
			os << '\t' << "rho";
			os << '\t' << "statvar";
			os << '\t' << "genealpha2";
		}
		else	{
			os << '\t' << "genealpha";
		}
		os << '\t' << "genevar";
		os << '\t' << "branchvar";

		os << '\t' << "rootbgcoffset";
		os << '\t' << "rootbgc";
		os << '\t' << "a1";
		os << '\t' << "b1";
		os << '\t' << "a2";
		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont; l++)	{
				os << '\t' << "cont_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "cont_" << k << '_' << k;
		}
		os << '\t' << "bgc_sigma";
		if (autoregressive)	{
			os << '\t' << "phi";
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}

		os << "\tdim";
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}

		if (triplet)	{
			if (Reversible())	{
				os << "\tctlambda\tat2cg\tat2gc\tat2ta\tgc2cg";
			}
			else	{
				os << "\tctat2cg\tat2gc\tat2ta\tcg2at\tcg2gc\tcg2ta";
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tnumerror";
		}
		os << "\tbgcnumerror";

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os.precision(10);

		// os << GetPruningTime() << '\t' << GetResampleTime() << '\t';
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		os << '\t' << synsigma->val();
		if (Reversible())	{
			os << '\t' << lambda->val();
			os << '\t' << at2cg->val();
			os << '\t' << at2gc->val();
			os << '\t' << at2ta->val();
			os << '\t' << gc2cg->val();
		}
		else	{
			os << '\t' << at2cg->val();
			os << '\t' << at2gc->val();
			os << '\t' << at2ta->val();
			os << '\t' << cg2at->val();
			os << '\t' << cg2gc->val();
			os << '\t' << cg2ta->val();
		}
		if (triplet)	{
			os << '\t' << cpgrate->val();
		}
		if (timelambda)	{
			os << '\t' << GetMeanLambda();
		}

		if (alpha)	{
			os << '\t' << alpha->val();
		}
		else	{
			os << '\t' << 0;
		}
		os << '\t' << logbgcoffset->val();
		for (int i=0; i<Ncont; i++)	{
			os << '\t' << GetRegCoef(i);
		}
		if (withalphaprocess || timelambda)	{
			for (int i=0; i<Ncont; i++)	{
				os << '\t' << GetAlphaRegCoef(i);
			}
		}


		if (! geneprocesstype)	{
			os << '\t' << rho->val();
			if (flexrho)	{
				os << '\t' << rhoprocess->GetMeanLogRate();
				os << '\t' << rhoprocess->GetVarLogRate();
				os << '\t' << rhosigma->val();
			}
			os << '\t' << genealpha->val() / 2 / rho->val();
		}
		else if (geneprocesstype == 4)	{
			os << '\t' << rho->val();
			os << '\t' << genealpha->val() / 2 / rho->val();
			os << '\t' << genealpha2->val();
		}
		else	{
			os << '\t' << genealpha->val();
		}
		os << '\t' << GetGeneVarCoeff();
		os << '\t' << GetBranchVarCoeff();

		os << '\t' << GetMeanRootBGCOffset();
		os << '\t' << GetMeanRootBGC();

		os << '\t' << *alpha1;
		os << '\t' << *beta1;
		os << '\t' << *alpha2;
		// os << '\t' << *beta2;

		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont; l++)	{
				os << '\t' << (*contsigma)[k][l];
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << (*contsigma)[k][k];
		}
		os << '\t' << (*bgcsigma)[0][0];

		if (autoregressive)	{
			os << '\t' << phi->val();
		}
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}
		os << '\t' << GetContProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();

		if (triplet)	{
			if (Reversible())	{
				os << '\t' << ctlambda->val();
				os << '\t' << ctat2cg->val();
				os << '\t' << ctat2gc->val();
				os << '\t' << ctat2ta->val();
				os << '\t' << ctgc2cg->val();
			}
			else	{
				os << '\t' << ctat2cg->val();
				os << '\t' << ctat2gc->val();
				os << '\t' << ctat2ta->val();
				os << '\t' << ctcg2at->val();
				os << '\t' << ctcg2gc->val();
				os << '\t' << ctcg2ta->val();
			}
		}

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\t' << GetTotalBGCOverflowCount();
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		os << *synsigma << '\n';
		os << *synratetree << '\n';
		os << *ContDiagArray << '\n';
		os << *BGCDiagArray << '\n';
		os << *contsigma << '\n';
		os << *bgcsigma << '\n';
		os << *regarray << '\n';
		if (autoregressive)	{
			os << *phi << '\n';
		}
		os << *contprocess << '\n';
		os << *logbgcprocess << '\n';
		if (Reversible())	{
			os << *lambda << '\n';
			os << *at2cg << '\n';
			os << *at2gc << '\n';
			os << *at2ta << '\n';
			os << *gc2cg << '\n';
			if (triplet)	{
				os << *ctlambda << '\n';
				os << *ctat2cg << '\n';
				os << *ctat2gc << '\n';
				os << *ctat2ta << '\n';
				os << *ctgc2cg << '\n';
			}
		}
		else	{
			os << *at2cg << '\n';
			os << *at2gc << '\n';
			os << *at2ta << '\n';
			os << *cg2at << '\n';
			os << *cg2gc << '\n';
			os << *cg2ta << '\n';
			if (triplet)	{
				os << *ctat2cg << '\n';
				os << *ctat2gc << '\n';
				os << *ctat2ta << '\n';
				os << *ctcg2at << '\n';
				os << *ctcg2gc << '\n';
				os << *ctcg2ta << '\n';
			}
		}
		if (triplet)	{
			os << *cpgrate << '\n';
		}
		if (alpha)	{
			os << *alpha << '\n';
		}
		os << *logbgcoffset << '\n';
		os << *rho << '\n';
		if (flexrho)	{
			os << *rhosigma << '\n';
			os << *rhoprocess << '\n';
		}
		os << *genealpha << '\n';
		os << *alpha1 << '\n';
		os << *beta1 << '\n';
		os << *alpha2 << '\n';
		// os << *beta2 << '\n';
		for (int gene=0; gene<Ngene; gene++)	{
			os << *rootbgcoffset[gene] << '\n';
			os << *rootbgc[gene] << '\n';
			if (geneprocesstype == 1)	{
				os << *gammageneprocess[gene] << '\n';
			}
			else if (geneprocesstype == 2)	{
				os << *wngeneprocess[gene] << '\n';
			}
			else if (geneprocesstype == 3)	{
				os << *lngeneprocess[gene] << '\n';
			}
			else if (geneprocesstype == 4)	{
				os << *genealpha2 << '\n';
				os << *argeneprocess[gene] << '\n';
				os << *wngeneprocess[gene] << '\n';
			}
			else	{
				os << *argeneprocess[gene] << '\n';
			}
		}
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			is >> *Chi >> *Chi2;
		}
		is >> *synsigma;
		is >> *synratetree;
		is >> *ContDiagArray;
		is >> *BGCDiagArray;
		is >> *contsigma;
		is >> *bgcsigma;
		is >> *regarray;
		if (autoregressive)	{
			is >>  *phi;
		}
		is >> *contprocess;
		is >> *logbgcprocess;
		if (Reversible())	{
			is >> *lambda;
			is >> *at2cg;
			is >> *at2gc;
			is >> *at2ta;
			is >> *gc2cg;
			if (triplet)	{
				is >> *ctlambda;
				is >> *ctat2cg;
				is >> *ctat2gc;
				is >> *ctat2ta;
				is >> *ctgc2cg;
			}
		}
		else	{
			is >> *at2cg;
			is >> *at2gc;
			is >> *at2ta;
			is >> *cg2at;
			is >> *cg2gc;
			is >> *cg2ta;
			if (triplet)	{
				is >> *ctat2cg;
				is >> *ctat2gc;
				is >> *ctat2ta;
				is >> *ctcg2at;
				is >> *ctcg2gc;
				is >> *ctcg2ta;
			}
		}
		if (triplet)	{
			is >> *cpgrate;
		}
		if (alpha)	{
			is >> *alpha;
		}
		is >> *logbgcoffset;

		is >> *rho;
		if (flexrho)	{
			is >> *rhosigma;
			is >> *rhoprocess;
		}
		is >> *genealpha;
		is >> * alpha1;
		is >> *beta1;
		is >> * alpha2;
		// is >> *beta2;
		for (int gene=0; gene<Ngene; gene++)	{
			is >> *rootbgcoffset[gene];
			is >> *rootbgc[gene];
			if (geneprocesstype == 1)	{
				is >> *gammageneprocess[gene];
			}
			else if (geneprocesstype == 2)	{
				is >> *wngeneprocess[gene];
			}
			else if (geneprocesstype == 3)	{
				is >> *lngeneprocess[gene];
			}
			else if (geneprocesstype == 4)	{
				is >> *genealpha2;
				is >> *argeneprocess[gene];
				is >> *wngeneprocess[gene];
			}
			else	{
				is >> *argeneprocess[gene];
			}
		}
	}

};

#endif
