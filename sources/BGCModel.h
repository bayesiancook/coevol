
#ifndef BGC
#define BGC

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSequenceAlignment.h"
#include "RYSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "CoalCalibratedChronogram.h"

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

class BGCModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;
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
	Const<PosReal>* Ten;
	Const<PosReal>* Tenth;
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
	// 2 : bd with cauchy proper lower bounds
	// 3 : rbd with cauchy proper lower bounds

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
	Jeffreys* phi;

	Jeffreys* alpha;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* contprocess;
	LinRegNormalProcess* logbgcprocess;
	Normal* logbgcoffset;
	BGCProcess* bgcprocess;
	Jeffreys* lambda;
	Jeffreys* at2cg;
	Jeffreys* at2gc;
	Jeffreys* at2ta;
	Jeffreys* gc2cg;
	Jeffreys* rootbgc;
	RandomMutSubMatrix* mutmatrix;
	BGCMatrixTree* nucmatrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	bool meanexp;

	bool conjpath;
	bool priorsampling;

	bool normalise;

	int nrep;

	int df;
	int Ninterpol;

	int clampsuffstat;
	string suffstatfile;

	bool clamptree;
	bool clampdiag;
	bool autoregressive;

	double fixalpha;
	bool withalphaprocess;

	int Nreg;

	MeanExpTreeFromMultiVariate* alphaprocess;

	public:

	SequenceAlignment* GetData()	{
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

	BGCModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, int indf, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inNinterpol, string insuffstatfile, string rootfile, double infixalpha, double inlambda, double inat2cg, double inat2gc, double ingc2cg, bool sample=true)	{

		clamptree = inclamptree;
		fixalpha = infixalpha;
		if (fixalpha == -2)	{
			fixalpha = -1;
			withalphaprocess = true;
			Nreg = 2;
		}
		else	{
			withalphaprocess = false;
			Nreg = 1;
		}


		autoregressive = inautoregressive;
		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		Ninterpol = inNinterpol;
		df = indf;

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		meanexp = inmeanexp;

		// get data from file

		nucdata = new FileSequenceAlignment(datafile);
		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

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
			nrep = conjpath ? 30 : 1;
		}
		cerr << "nrep : " << nrep << '\n';
		cerr << conjpath << '\n';
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
		Ten = new Const<PosReal>(10);
		Tenth = new Const<PosReal>(0.1);

		RootAlpha = 0;
		RootBeta = 0;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		if (calibfile != "None")	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,One,RootAlpha,RootBeta,calibset);
			}
			else {
				cerr << "BD\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,RootAlpha,RootBeta,calibset,chronoprior);
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
		}
		else	{
			phi = 0;
		}
		regarray->SetAtZero();
		/*
		regarray->GetIIDNormal(0)->ClampAt(-1,1);
		regarray->GetIIDNormal(0)->ClampAt(1,2);
		*/
		// regarray->GetIIDNormal(0)->ClampAt(0,1);
		// regarray->GetIIDNormal(0)->ClampAt(0,2);

		cerr << "bgc process\n";

		logbgcprocess = new LinRegNormalProcess(lengthtree,contprocess,regarray,bgcsigma,phi);
		logbgcprocess->CutOff(0.1,0);
		(*logbgcprocess->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode()))[0] = 0;
		if (withalphaprocess)	{
			logbgcprocess->CutOff(0.1,1);
			(*logbgcprocess->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode()))[1] = 0;
		}
		logbgcprocess->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode())->Clamp();

		alpha = 0;
		if (fixalpha != 0)	{
			double min = 0.01;
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

		logbgcoffset = new Normal(Zero,Ten);
		logbgcoffset->setval(0);
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

		rootbgc = new Jeffreys(min,max,Zero);
		rootbgc->setval(1.0);

		cerr << "mutmatrix\n";
		mutmatrix = new RandomMutSubMatrix(lambda,at2cg,at2gc,at2ta,gc2cg,true);

		cerr << "bggcprocess\n";
		bgcprocess = new BGCProcess(logbgcprocess,logbgcoffset);

		if (alphaprocess)	{
			nucmatrixtree = new AlphaBGCMatrixTree(bgcprocess,mutmatrix,alphaprocess,alpha,rootbgc,normalise);
		}
		else	{
			nucmatrixtree = new BGCMatrixTree(bgcprocess,mutmatrix,alpha,rootbgc,normalise,0);
			// nucmatrixtree = new BGCMatrixTree(bgcprocess,mutmatrix,rootbgc,normalise);
		}

		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i,i,contdatatype);
			}
		}

		// make substitution mappings
		if (conjpath)	{
			pathconjtree = new BranchMatrixPathConjugateTree(synratetree, nucmatrixtree, GetData());

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
			phyloprocess = new BranchMatrixPhyloProcess(synratetree, nucmatrixtree, GetData());
		}

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
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(Ten);
		RootRegister(Tenth);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		else if (chronoprior)	{
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
			if (RootAlpha)	{
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
		total += lambda->GetLogProb();
		total += at2cg->GetLogProb();
		total += at2gc->GetLogProb();
		total += at2ta->GetLogProb();
		total += gc2cg->GetLogProb();
		total += logbgcoffset->GetLogProb();
		total += rootbgc->GetLogProb();
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

			if (RootAlpha && (! clamptree))	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
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

			scheduler.Register(new SimpleMove(rootbgc,1),10,"root bgc");
			scheduler.Register(new SimpleMove(rootbgc,0.1),10,"root bgc");
			scheduler.Register(new SimpleMove(rootbgc,0.01),10,"root bgc");

			if (alpha)	{
				scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
				scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
				scheduler.Register(new SimpleMove(alpha,0.01),10,"alpha");
			}

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
		lambda->Sample();
		at2cg->Sample();
		at2gc->Sample();
		at2ta->Sample();
		gc2cg->Sample();
		logbgcoffset->Sample();

		// gcprocess->specialUpdate();

		rootbgc->Sample();

		// gctree->specialUpdate()

		if (phyloprocess)	{
			phyloprocess->Sample();
		}

		cerr << "ok\n";
	}

	Var<PosReal>* GetScale()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	BGCProcess* GetBGCProcess()	{
		return bgcprocess;
	}

	LogNormalTreeProcess* GetSynRateTree()	{
		return synratetree;
	}

	double GetMeanSynRate()	{
		return GetSynRateTree()->GetMeanRate();
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	int GetTotalBGCOverflowCount()	{
		return nucmatrixtree->GetBGCOverflowCount();
	}

	/*
	double GetMeanGC()	{
		return gctree->GetMeanGCContent();
	}

	double GetVarGC()	{
		return gctree->GetVarGCContent();
	}
	*/

	double GetRootBGC()	{
		return rootbgc->val();
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
		os << "\tlambda\tat2cg\tat2gc\tat2ta\tgc2cg\talpha\tbgcoffset";
		os << '\t' << "reg0";
		os << '\t' << "reg1";
		os << '\t' << "reg2";
		if (withalphaprocess)	{
			os << '\t' << "areg0";
			os << '\t' << "areg1";
			os << '\t' << "areg2";
		}
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
		if (RootAlpha)	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}

		os << "\tdim";
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}

		os << '\t' << "rootbgc";
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tnumerror";
		}
		os << "\tbgcnumerror";

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		os << '\t' << synsigma->val();
		os << '\t' << lambda->val();
		os << '\t' << at2cg->val();
		os << '\t' << at2gc->val();
		os << '\t' << at2ta->val();
		os << '\t' << gc2cg->val();
		if (alpha)	{
			os << '\t' << alpha->val();
		}
		else	{
			os << '\t' << 0;
		}
		os << '\t' << logbgcoffset->val();
		os << '\t' << GetRegCoef(0);
		os << '\t' << GetRegCoef(1);
		os << '\t' << GetRegCoef(2);
		if (withalphaprocess)	{
			os << '\t' << GetAlphaRegCoef(0);
			os << '\t' << GetAlphaRegCoef(1);
			os << '\t' << GetAlphaRegCoef(2);
		}

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
		if (RootAlpha)	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}
		os << '\t' << GetContProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();
		os << '\t' << rootbgc->val();

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}
		os << '\t' << GetTotalBGCOverflowCount();

		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (RootAlpha)	{
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
		os << *lambda << '\n';
		os << *at2cg << '\n';
		os << *at2gc << '\n';
		os << *at2ta << '\n';
		os << *gc2cg << '\n';
		if (alpha)	{
			os << *alpha << '\n';
		}
		os << *logbgcoffset << '\n';
		os << *rootbgc << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (RootAlpha)	{
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
		is >> *lambda;
		is >> *at2cg;
		is >> *at2gc;
		is >> *at2ta;
		is >> *gc2cg;
		if (alpha)	{
			is >> *alpha;
		}
		is >> *logbgcoffset;
		is >> *rootbgc;
	}

};

#endif
