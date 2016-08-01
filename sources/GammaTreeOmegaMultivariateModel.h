
#ifndef GAMMATREEOMEGAMULTI
#define GAMMATREEOMEGAMULTI


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "CalibratedChronogram.h"
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

	GCMatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
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
	NucMatrixTree* nucmatrixtree;
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


class GammaTreeOmegaMultivariateModel : public ProbModel {

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

	// tree and branch lengths
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	GammaTree* gamtree;

	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;
	Gamma* phi;
	IIDNormal* mean;

	MultiVariateTreeProcess* process;

	MeanExpTreeFromMultiVariate* omegatree;

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

	// if true: gc model
	bool gc;

	// total number of substitution parameters modelled as non homogeneous
	int L;


	public:

	GammaTreeOmegaMultivariateModel() {
	}

	GammaTreeOmegaMultivariateModel(string datafile, string treefile, string contdatafile, double priorsigma, bool ingc = false, bool inclampdiag = false, bool inautoregressive = false, bool sample=true, GeneticCodeType type=Universal)	{

		bool meanexp = false;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		gc = ingc;
		L = gc ? 2 : 1;

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

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

		PriorLambda = new Const<PosReal>(10);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,One,lambda);

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
			mean = new IIDNormal(Ncont+L,Zero,One);
			process = new AutoRegressiveMultiVariateTreeProcess(sigma,mean,phi,gamtree);
		}
		else	{
			phi = 0;
			mean = 0;
			process = new MultiVariateTreeProcess(sigma,gamtree);
		}

		// process->ClampTransversal(1,0);

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,L+i,i,0);
			}
		}

		// cut off to avoid numerical errors
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		omegatree = new MeanExpTreeFromMultiVariate(process,0,MEAN,false,meanexp);

		// make rate matrices
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (gc)	{
			gctree = new MeanLogitTreeFromMultiVariate(process,2,MEAN,false);
			rootgc = new Beta(One,One);
			stattree = new GCStatTree(gctree,rootgc);
			nucmatrixtree = new NucMatrixTree(relrate,stattree);
			matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);

			nucmatrix = 0;
			stationary = 0;
		}
		else	{
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);
			matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

			gctree = 0;
			rootgc = 0;
			stattree = 0;
			nucmatrixtree = 0;
		}

		// make substitution mappings
		pathconjtree = new BranchMatrixPathConjugateTree(gamtree, matrixtree, codondata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		// register model
		RootRegister(PriorLambda);
		RootRegister(One);
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
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~GammaTreeOmegaMultivariateModel() {}

	Tree* GetTree() {return tree;}

	MeanExpTreeFromMultiVariate* GetOmegaTree() {return omegatree;}
	MeanLogitTreeFromMultiVariate* GetGCTree() {return gctree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	GammaTree* GetGammaTree() {return gamtree;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

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

		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();

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

		// scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<30; i++)	{
			scheduler.Register(new SimpleMove(lambda,1),100,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.1),100,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.01),100,"lambda");
			scheduler.Register(new SimpleMove(gamtree,1),10,"bl");
			scheduler.Register(new SimpleMove(gamtree,0.1),10,"bl");
			scheduler.Register(new SimpleMove(gamtree,0.01),10,"bl");

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
		}
	}

	void drawSample()	{
		cerr << "sample\n";

		lambda->Sample();
		gamtree->Sample();
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

		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		omegatree->specialUpdate();

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

	double GetLength()	{
		return gamtree->GetTotalLength();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
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
		os << "#logprior\tlnL\tlength\tomega";

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
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

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();
		os << '\t' << GetMeanOmega();

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
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

	int GetNcheck()	{
		int tmp =  4 + (Ncont + L) * (Ncont + L+1) / 2;
		if (autoregressive)	{
			tmp += 3;
		}
		return tmp;
	}

	void GetCheck(double* p)	{
		p[0] = GetLength();
		GetOmegaTree()->specialUpdate();
		p[1] = GetMeanOmega();
		int n = 2;
		if (autoregressive)	{
			p[2] = phi->val();
			p[3] = mean->GetMean();
			p[4] = mean->GetVar();
			n = 5;
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
		if (gc)	{
			p[n] =  GetMeanGC();
			n++;
		}
		else	{
			p[n] =  stationary->val().GetEntropy();
			n++;
		}
		p[n] = relrate->val().GetEntropy();
		n++;
		if (n!= GetNcheck())	{
			cerr << "error : non matching Ncheck\n";
			cerr << n << '\t' << GetNcheck() << '\n';
			cerr << Ncont << '\n';
			exit(1);
		}
	}

	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *gamtree << '\n';
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
		is >> *lambda;
		is >> *gamtree;
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

	void Simulate(string name)	{
		lambda->Sample();
		gamtree->Sample();
		sigma->Sample();
		if (autoregressive)	{
			phi->Sample();
			// phi->setval(1);
			mean->Sample();
			/*
			for (int k=0; k<Ncont+L; k++)	{
				(*mean)[k] = 0;
			}
			*/
		}
		for (int i=0; i<Ncont; i++)	{
			process->Unclamp(i+L);
		}
		process->Sample();
		if (gc)	{
			rootgc->Sample();
		}
		else	{
			stationary->Sample();
		}
		relrate->Sample();
		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample();
		phyloprocess->GetLeafData(datacopy);
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
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

};

#endif
