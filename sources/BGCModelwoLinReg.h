
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

#include "MeanJitteredLogitTree.h"

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
	InverseWishartMatrix* contsigma;
	// ConjugateInverseWishart* contsigma;
	// Rvar<CovMatrix>* contsigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* contprocess;
	// ConjugateMultiVariateTreeProcess* contprocess;
	Normal* logbgcoffset;
	BGCProcess* bgcprocess;
	Jeffreys* lambda;
	Jeffreys* kappa;
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

	BGCModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, int indf, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inNinterpol, string insuffstatfile, string rootfile, bool sample=true)	{

		clamptree = inclamptree;

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

		ContDiagArray = new JeffreysIIDArray(Ncont+1,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			ContDiagArray->setval(1.0);
		}
		else	{
			ContDiagArray->ClampAt(priorsigma);
		}

		ContSigmaZero = new SigmaZero(ContDiagArray);

		if (clampdiag)	{
			// contsigma = new DiagonalCovMatrix(ContSigmaZero,Ncont+1+df);
		}
		else	{
			// contsigma = new ConjugateInverseWishart(ContSigmaZero,Ncont+1+df);
			contsigma = new InverseWishartMatrix(ContSigmaZero,Ncont+1+df);
		}

		rootmean = 0;
		rootvar = 0;
		cerr << rootfile << '\n';
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(Ncont+1));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont+1));
			ifstream is(rootfile.c_str());
			for (int i=0; i<Ncont+1; i++)	{
				double mean, var;
				is >> mean >> var;
				(*rootmean)[i] = mean;
				(*rootvar)[i] = var;
			}
		}

		if (clampdiag)	{
			// contprocess = new MultiVariateTreeProcess(contsigma,lengthtree,0,0,rootmean,rootvar);
		}
		else	{
			// contprocess = new ConjugateMultiVariateTreeProcess(contsigma,lengthtree,0,0,rootmean, rootvar);
			contprocess = new MultiVariateTreeProcess(contsigma,lengthtree,0,0,rootmean, rootvar);
		}

		double min = 0.001;
		double max = 1000;

		cerr << "bgc process\n";

		logbgcoffset = new Normal(Zero,Ten);
		logbgcoffset->setval(0);
		lambda = new Jeffreys(min,max,Zero);
		lambda->ClampAt(2.0);
		// lambda->ClampAt(4.0);
		kappa = new Jeffreys(min,max,Zero);
		kappa->setval(1.0);

		rootbgc = new Jeffreys(min,max,Zero);
		rootbgc->setval(1.0);

		cerr << "mutmatrix\n";
		mutmatrix = new RandomMutSubMatrix(lambda,kappa,true);

		cerr << "bggcprocess\n";
		bgcprocess = new BGCProcess(contprocess,logbgcoffset,0);

		nucmatrixtree = new BGCMatrixTree(bgcprocess,mutmatrix,rootbgc,normalise);

		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i+1,i,contdatatype);
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
		total += contsigma->GetLogProb();
		total += contprocess->GetLogProb();
		total += lambda->GetLogProb();
		total += kappa->GetLogProb();
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

			/*
			scheduler.Register(new ConjugateMultiVariateMove(contsigma,contprocess,10,10),1,"conj sigma process");
			scheduler.Register(new ConjugateMultiVariateMove(contsigma,contprocess,1,10),1,"conj sigma process");
			scheduler.Register(new ConjugateMultiVariateMove(contsigma,contprocess,0.1,10),1,"conj sigma process");
			scheduler.Register(new ConjugateMultiVariateMove(contsigma,contprocess,0.01,10),1,"conj sigma process");
			*/

			scheduler.Register(new SimpleMove(contprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(contsigma,10),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.01),100,"cont sigma");

			scheduler.Register(new SimpleMove(ContDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,0.1),10,"theta");

			scheduler.Register(new SimpleMove(rootbgc,1),10,"root bgc");
			scheduler.Register(new SimpleMove(rootbgc,0.1),10,"root bgc");
			scheduler.Register(new SimpleMove(rootbgc,0.01),10,"root bgc");

			scheduler.Register(new SimpleMove(kappa,1),10,"kappa");
			scheduler.Register(new SimpleMove(kappa,0.1),10,"kappa");
			scheduler.Register(new SimpleMove(kappa,0.01),10,"kappa");

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

		contsigma->Sample();
		// sigma->SetIdentity();

		contprocess->Sample();
		contprocess->CutOff(1,0);
		lambda->Sample();
		kappa->Sample();
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

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tsynrate";

		os << "\tsynsigma";
		os << "\tlambda\tkappa\tbgcoffset";
		for (int k=0; k<Ncont+1; k++)	{
			for (int l=k+1; l<Ncont+1; l++)	{
				os << '\t' << "cont_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+1; k++)	{
			os << '\t' << "cont_" << k << '_' << k;
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
		os << '\t' << kappa->val();
		os << '\t' << logbgcoffset->val();

		for (int k=0; k<Ncont+1; k++)	{
			for (int l=k+1; l<Ncont+1; l++)	{
				os << '\t' << (*contsigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+1; k++)	{
			os << '\t' << (*contsigma)[k][k];
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
		os << *contsigma << '\n';
		os << *contprocess << '\n';
		os << *lambda << '\n';
		os << *kappa << '\n';
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
		is >> *contsigma;
		is >> *contprocess;
		is >> *lambda;
		is >> *kappa;
		is >> *logbgcoffset;
		is >> *rootbgc;
	}

};

#endif
