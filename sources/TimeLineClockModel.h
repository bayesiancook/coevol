
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "BDCalibratedChronogram.h"
#include "CoalCalibratedChronogram.h"

#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GeneralConjugatePath.h"

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
#include "GenericTimeLine.h"
#include "StatProcess.h"


class LGMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	LGMatrixTree(BranchVarTree<Profile>* instattree, bool innormalise = true) {
		SetWithRoot(true);
		stattree = instattree;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new LGRandomSubMatrix(stattree->GetBranchVal(link->GetBranch()),normalise);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;
	bool normalise;

};

class TimeLineClockModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;

	SequenceAlignment* data;
	// ProteinSequenceAlignment* proteindata;

	ContinuousData* contdata;
	AncestralData* ancdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ncont;
	int Nanc;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;
	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;
	bool calibchrono;

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
	NodeBranchVarTree<PosReal,PosReal>* finalchrono;
	Var<PosReal>* scale;
	Gamma* freescale;

	// if different branches have different scaling factors
	BranchPartition* partition;
	SplitBranchPartition* splitpartition;
	GammaMixTree* gammamixtree;
	string mix;

	int Nmat;
	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	IIDArray<CovMatrix>* sigmaarray;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	LinearlyInterpolatedTimeLine* timeline;
	Normal* alpha;
	int flextimeline;

	MultiVariateTreeProcess* process;

	// is it a copmponent of the multivariate tree process
	// or a separate lognormal process?
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;
	LengthTree* lengthtree;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;

	// for homogeneous model
	Dirichlet* stationary;
	RandomSubMatrix* lgmatrix;

	// for both
	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	StatTree* stattree;
	BranchStatTree* branchstattree;

	// BranchMatrixRASPhyloProcess* phyloprocess;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamproot;
	bool clamptree;

	bool meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int L;
	int M;
	int Kaa;

	bool conjpath;
	bool priorsampling;

	bool normalise;

	int nrep;

	int df;
	int Ninterpol;

	string bounds;
	bool withdrift;
	bool withexpdrift;
	bool withreldrift;

	int clampsuffstat;
	string suffstatfile;

	bool autoregressive;

	bool freestat;

	public:

	SequenceAlignment* GetData()	{
		return data;
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

	TimeLineClockModel() {
	}

	BranchPartition* GetPartition()	{
		if (partition && Split())	{
			return splitpartition;
		}
		return partition;
	}


	TimeLineClockModel(string datafile, string treefile, string contdatafile, string ancdatafile, string timelinefile, int inflextimeline, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, int indf, bool inclampdiag, int inconjpath, int contdatatype, int ancdatatype, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, string inmix, int inNinterpol, string insuffstatfile, string rootfile, bool sample=true, GeneticCodeType type=Universal)	{

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		Ninterpol = inNinterpol;
		mix = inmix;

		df = indf;

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;

		// get data from file

		data = new FileSequenceAlignment(datafile);
		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();
		if (Nstate != 20)	{
			cerr << "error : should be protein alignment\n";
			exit(1);
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
			nrep = conjpath ? 30 : 1;
		}
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		taxonset = data->GetTaxonSet();

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

		// get ancestral data from file
		if (ancdatafile != "None")	{
			cerr << ancdatafile << "\n";
			ancdata = new FileAncestralData(tree,ancdatafile);
			Nanc = ancdata->GetNsite();
			cerr << "ok\n";
		}
		else	{
			ancdata = 0;
			Nanc = 0;
		}

		cerr << "tree and data ok\n";
		cerr << '\n';

		// number of substitution parameters
		L = 0;
		freestat = (!Nanc);
		if (freestat)	{
			L = Nstate;
		}
		// total dimension of the process
		M = L + Ncont + Nanc;
		Kaa = Ncont + Nanc;

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		calibchrono = false;
		RootAlpha = 0;
		RootBeta = 0;
		freescale = 0;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		if (rootage != 0)	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
		}
		if (calibfile != "None")	{
			calibchrono = true;
			/*
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			*/
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
			scale = GetCalibratedChronogram()->GetScale();
		}
		else	{
			chronogram = new Chronogram(tree,One);
			if (rootage)	{
				freescale = new Gamma(RootAlpha,RootBeta);
				scale = freescale;
			}
			else	{
				scale = One;
			}
		}

		cerr << calibfile << '\t' << RootAlpha << '\n';
		if (clamptree)	{
			chronogram->Clamp();
		}

		if (Split())	{
			finalchrono = new SplitChronogram(chronogram,GetSplitTree());
		}
		else	{
			finalchrono = chronogram;
		}

		// ratetree

		// a log normal process on that tree
		sigma = new Gamma(One,One);
		// sigma->ClampAt(1);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		lengthtree = lognormaltree;
		// alternatively:
		// L++;
		// lengthtree = new MeanExpTreeFromMultiVariate(lengthtree,subprocess,0,INTEGRAL,false,meanexp);

		if (mix != "None")	{
			partition = new BranchPartition(tree,mix);
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
			// gammamixtree = new GammaMixTree(partition,One,One,true); // clamp first component to 1
		}
		else	{
			partition = new BranchPartition(tree);
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
		}

		Nmat = GetPartition()->GetNComponent();

		double mindiag = 0.001;
		double maxdiag = 1000;

		DiagArray = new JeffreysIIDArray(M,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}

		sigmaZero = new SigmaZero(DiagArray);
		cerr << M << '\n';

		cerr << "sigmaarray\n";
		if (clampdiag)	{
			sigmaarray = new DiagonalWishartArray(Nmat, df, sigmaZero, sigmaZero);
		}
		else	{
			sigmaarray = new InverseWishartArray(Nmat, df, sigmaZero, sigmaZero);
		}

		cerr << "root\n";
		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(M));
			rootvar = new Const<PosRealVector>(PosRealVector(M));
			ifstream is(rootfile.c_str());
			for (int i=0; i<M; i++)	{
				double mean, var;
				is >> mean >> var;
				(*rootmean)[i] = mean;
				(*rootvar)[i] = var;
			}
		}

		// to be completed
		flextimeline = inflextimeline;
		timeline = 0;
		alpha = 0;
		if (timelinefile != "None")	{
			timeline = new LinearlyInterpolatedTimeLine();
			ifstream is(timelinefile.c_str());
			timeline->FromStream(is);
			alpha = new Normal(Zero,One);
			if (! flextimeline)	{
				alpha->ClampAt(1.0);
			}
		}

		cerr << "process\n";
		process = new PartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),0,0,0,finalchrono, rootmean, rootvar, 0, 0, GetScale(), 1 , timeline, alpha);
		cerr << "processok\n";

		if (contdata)	{
			cerr << "set and clamp continuous data\n";
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,i,i,contdatatype);
			}
		}
		if (ancdata)	{
			cerr << "set and clamp ancestral data\n";
			process->SetAndClamp(ancdata,Ncont,ancdatatype);
			cerr << "ok\n";
		}

		cerr << "cutoff\n";
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		stationary = new Dirichlet(Nstate);
		stationary->setuniform();

		if (!Nanc)	{
			lgmatrix = 0;
			cerr << "stattree\n";
			stattree = new StatTree(stationary,process,Kaa);
			cerr << "branch stat tree\n";
			branchstattree = new BranchStatTree(stattree);
			matrixtree = new LGMatrixTree(branchstattree,normalise);

			// make substitution mappings
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(lengthtree, matrixtree, data);
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
					phyloprocess = new BranchMatrixPhyloProcess(lengthtree, matrixtree, data);
				}
			}
		}
		else	{
			stattree = 0;
			branchstattree = 0;
			matrixtree = 0;
			lgmatrix = new LGRandomSubMatrix(stationary,normalise);
			// make substitution mappings
			if (conjpath)	{
				pathconjtree = new OneMatrixPathConjugateTree(lengthtree, lgmatrix, data);
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
					phyloprocess = new OneMatrixPhyloProcess(lengthtree, lgmatrix, data);
				}
			}
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
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		RootRegister(stationary);
		if (rootmean)	{
			RootRegister(rootmean);
			RootRegister(rootvar);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
			if (calibchrono)	{
				cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
				cerr << "root age : " << GetScale()->val() << '\n';
			}
		}

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~TimeLineClockModel() {}


	Tree* GetTree() {return tree;}
	Tree* GetFineGrainedTree() {return splittree;}

	MultiVariateTreeProcess* GetProcess() {return process;}
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
	AncestralData* GetAncestralData() {return ancdata;}

	int GetL() {return L;}
	int GetNcont() {return Ncont;}
	double GetAlpha()	{
		return alpha->val();
	}

	CovMatrix* GetCovMatrix(int mat = 0) {return sigmaarray->GetVal(mat);}

	CalibratedChronogram* GetCalibratedChronogram()	{
		CalibratedChronogram* tmp = dynamic_cast<CalibratedChronogram*>(chronogram);
		if (!tmp)	{
			cerr << "error in get calibrated chrono\n";
			exit(1);
		}
		return tmp;
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

		if (gammamixtree)	{
			total += gammamixtree->GetLogProb();
		}

		total += DiagArray->GetLogProb();
		total += sigmaarray->GetLogProb();
		total += process->GetLogProb();

		total += stationary->GetLogProb();
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

		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(sigma,1),10,"nu");
		scheduler.Register(new SimpleMove(sigma,0.1),10,"nu");
		scheduler.Register(new SimpleMove(sigma,0.01),10,"nu");

		if (timeline)	{
			scheduler.Register(new SimpleMove(alpha,1),10,"timeline alpha");
			scheduler.Register(new SimpleMove(alpha,0.1),10,"timeline alpha");
			scheduler.Register(new SimpleMove(alpha,0.01),10,"timeline alpha");
		}

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

			if (gammamixtree)	{
				scheduler.Register(new SimpleMove(gammamixtree,10),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,1),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,0.1),10,"gammamixtree");
			}
			if (RootAlpha && (! clamptree))	{
				if (freescale)	{
					scheduler.Register(new SimpleMove(freescale,1),10,"root age");
					scheduler.Register(new SimpleMove(freescale,0.1),10,"root age");
					scheduler.Register(new SimpleMove(freescale,0.01),10,"root age");
				}
				else	{
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				}
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

			scheduler.Register(new SimpleMove(process,10),10,"multinormal");
			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(sigmaarray,10),100," sigma");
			scheduler.Register(new SimpleMove(sigmaarray,1),100," sigma");
			scheduler.Register(new SimpleMove(sigmaarray,0.1),100," sigma");
			scheduler.Register(new SimpleMove(sigmaarray,0.01),100," sigma");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
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
		if (freescale)	{
			freescale->Sample();
		}
		if (gammamixtree)	{
			gammamixtree->Sample();
		}
		DiagArray->Sample();

		sigmaarray->Sample();
		process->Sample();

		if (timeline)	{
			alpha->Sample();
		}
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		lognormaltree->specialUpdate();
		stationary->Sample();

		if (phyloprocess)	{
			phyloprocess->Sample();
		}

		cerr << "ok\n";
	}


	Var<PosReal>* GetScale()	{
		return scale;
		/*
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
		*/
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			if (freescale)	{
				return freescale->val();
			}
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tlength";

		if (stattree)	{
			os << "\tmeanstatent";
			os << "\tvarstatent";
		}
		os << "\trefstatent";
		if (RootAlpha)	{
			os << "\trootage";
		}
		if (timeline)	{
			os << "\talpha";
		}

		for (int mat=0; mat<Nmat; mat++)	{
			/*
			for (int k=0; k<M; k++)	{
				for (int l=k+1; l<M; l++)	{
					os << '\t' << "cov_" << mat << '_' << k << '_' << l;
				}
			}
			*/
			for (int k=0; k<M; k++)	{
				os << '\t' << "cov_" << mat << '_' << k << '_' << k;
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}

		os << "\tdim";
		for (int k=0; k<M; k++)	{
			os << '\t' << "root_" << k;
		}

		if (gammamixtree)	{
			for (int k=0; k<gammamixtree->GetNComponent(); k++)	{
				os << "\tgammamix" << k ;
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tnumerror";
		}

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();

		if (stattree)	{
			os << '\t' << stattree->GetMeanEntropy();
			os << '\t' << stattree->GetVarEntropy();
		}
		os << '\t' << stationary->GetEntropy();
		if (RootAlpha)	{
			os << '\t' << GetRootAge();
		}
		if (timeline)	{
			os << '\t' << alpha->val();
		}

		for (int mat=0; mat<Nmat; mat++)	{
			/*
			for (int k=0; k<M; k++)	{
				for (int l=k+1; l<M; l++)	{
					os << '\t' << (*sigmaarray->GetVal(mat))[k][l];
				}
			}
			*/
			for (int k=0; k<M; k++)	{
				os << '\t' << (*sigmaarray->GetVal(mat))[k][k];
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}
		os << '\t' << GetProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();

		os << '\t' << stationary->val().GetEntropy();

		if (gammamixtree)	{
			for (int k=0; k<gammamixtree->GetNComponent(); k++)	{
				os << '\t' << gammamixtree->GetComponent(k)->val();
			}
		}

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (RootAlpha)	{
			if (freescale)	{
				os << *freescale << '\n';
			}
			else	{
				os << *GetCalibratedChronogram()->GetScale() << '\n';
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		if (gammamixtree)	{
			os << *gammamixtree << '\n';
		}
		if (timeline)	{
			os << *alpha << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigmaarray << '\n';
		os << *process << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (RootAlpha)	{
			if (freescale)	{
				is >> *freescale;
			}
			else	{
				is >> *GetCalibratedChronogram()->GetScale();
			}
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			is >> *Chi >> *Chi2;
		}
		if (gammamixtree)	{
			is >> *gammamixtree;
		}
		if (timeline)	{
			is >> *alpha;
		}
		is >> *DiagArray;
		is >> *sigmaarray;
		is >> *process;
		is >> *stationary;
	}
};

#endif
