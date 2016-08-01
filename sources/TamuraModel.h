

#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "MG3OmegaCodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"
#include "TamuraSubMatrix.h"

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

class TamuraMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	TamuraMatrixTree(BranchVarTree<PosReal>* intstvtree, BranchVarTree<PosReal>* ingctree, BranchVarTree<PosReal>* inattree, Var<UnitReal>* inrootgc, bool innormalise) {
		SetWithRoot(true);
		tstvtree = intstvtree;
		gctree = ingctree;
		attree = inattree;
		rootgc = inrootgc;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~TamuraMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return tstvtree->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new TamuraRandomSubMatrix(0,rootgc,normalise);
		}
		return new TamuraRandomSubMatrix(tstvtree->GetBranchVal(link->GetBranch()), gctree->GetBranchVal(link->GetBranch()), attree->GetBranchVal(link->GetBranch()),normalise);
	}

	private:

	BranchVarTree<PosReal>* tstvtree;
	BranchVarTree<PosReal>* gctree;
	BranchVarTree<PosReal>* attree;
	Var<UnitReal>* rootgc;
	bool normalise;

};

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};


class GCMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GCMatrixTree(CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* innucmatrixtree, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	Var<PosReal>* rootomega;

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


class TamuraModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
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
	Chronogram* chronogram;

	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;
	Gamma* phi;
	IIDUniform* mean;

	MultiVariateTreeProcess* process;

	MeanExpTreeFromMultiVariate* gcsynratetree;
	MeanExpTreeFromMultiVariate* atsynratetree;
	MeanExpTreeFromMultiVariate* tstvtree;
	MeanExpTreeFromMultiVariate* nonsynratetree;
	BranchValPtrTree<Dvar<PosReal> >* omegatree;

	Beta* rootgc;
	TamuraMatrixTree* nucmatrixtree;

	// for both
	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;
	// BranchMatrixPhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	// if true: OUP else BP
	bool autoregressive;

	bool clamproot;

	bool meanexp;

	int omegaratiotree;
	// 0 : rate model
	// 1 : dS omega
	// 2 : dS dN

	// total number of substitution parameters modelled as non homogeneous
	int L;

	bool conjpath;

	bool normalise;

	int nrep;

	public:

	SequenceAlignment* GetData()	{
		if (omegaratiotree)	{
			return codondata;
		}
		return nucdata;
	}

	TamuraModel() {
	}

	TamuraModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, int inomegaratiotree, bool inclamproot, bool inmeanexp, bool innormalise, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		normalise = innormalise;

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		omegaratiotree = inomegaratiotree;

		if (omegaratiotree)	{
			L = 4;
		}
		else	{
			L = 3;
		}
		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		if (omegaratiotree)	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else	{
			conjpath = inconjpath;
		}

		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 30 : 1;
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
			else if (chronoprior == 1)	{
				cerr << "BD\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,RootAlpha,RootBeta,calibset);
			}
			else	{
				cerr << "error : chronoprior\n";
				exit(1);
			}
		}
		else	{
			chronogram = new Chronogram(tree,One);
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

		process->Reset();

		// cut off to avoid numerical errors
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		rootgc = new Beta(One,One);
		gcsynratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
		atsynratetree = new MeanExpTreeFromMultiVariate(process,1,INTEGRAL,false,meanexp);
		tstvtree = new MeanExpTreeFromMultiVariate(process,2,MEAN,false,meanexp);
		nonsynratetree = 0;
		if (omegaratiotree == 0)	{
			omegatree = 0;
		}
		else if (omegaratiotree == 1)	{
			omegatree = new MeanExpTreeFromMultiVariate(process,3,MEAN,false,meanexp);
		}
		else	{
			nonsynratetree = new MeanExpTreeFromMultiVariate(process,3,INTEGRAL,false,meanexp);
			omegatree = new RatioTree(nonsynratetree,gcsynratetree);
		}

		nucmatrixtree = new TamuraMatrixTree(tstvtree,gcsynratetree,atsynratetree,rootgc,normalise);

		if (omegaratiotree)	{
			matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(gcsynratetree, matrixtree, codondata);
				phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(gcsynratetree, matrixtree, codondata);
			}
		}
		else	{
			matrixtree = 0;
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(gcsynratetree, nucmatrixtree, nucdata);
				phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(gcsynratetree, nucmatrixtree, nucdata);
			}
		}

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		// register model
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
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~TamuraModel() {}

	Tree* GetTree() {return tree;}

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

	MeanExpTreeFromMultiVariate* GetSynRateTree() {return gcsynratetree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	Chronogram* GetChronogram() {return chronogram;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double Update(bool check = false)	{
		double ret = ProbModel::Update();
		// phyloprocess->Sample();
		// ret = ProbModel::Update();
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

		total += rootgc->GetLogProb();
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

			scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
			scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
			scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");

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
		scheduler.Cycle(1,1,true,true);
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
		process->CutOff(2,1);
		gcsynratetree->specialUpdate();
		if (omegaratiotree)	{
			UpdateOmegaTree();
		}

		rootgc->Sample();

		phyloprocess->Sample();

		cerr << "ok\n";
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	double GetMaxdS()	{
		if (! omegaratiotree)	{
			cerr << "error : get max dS called under simple rate model\n";
			exit(1);
		}
		return gcsynratetree->GetMax();
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
		return gcsynratetree->GetTotal();
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

	Var<RealVector>* GetMeanVector()	{
		return mean;
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

		os << "\trootgc";

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

		os << '\t' << *rootgc;
		os << '\n';
		os.flush();
	}

	void PrintEntries(ostream& os)	{

		os << "Kgc\n";
		os << "Kat\n";
		os << "Ts/Tv\n";
		if (omegaratiotree == 2)	{
			os << "dN\n";
		}
		else if (omegaratiotree == 1) {
			os << "omega\n";
		}
		for (int k=0; k<Ncont; k++)	{
			os << "character " << k + 1 << '\n';
		}
	}

	void ToStream(ostream& os)	{
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
		os << *rootgc << '\n';
	}

	void FromStream(istream& is)	{
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
		is >> *rootgc;
	}
};

#endif
