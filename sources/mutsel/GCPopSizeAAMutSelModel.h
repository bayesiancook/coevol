
// a chronogram
//
// 2 independent log normal processes
// - log mutation rate (per unit of time)
// - log popsize (meant to be clamped at two leaf nodes)
// log(omega) = beta + alpha * log(popsize)

#ifndef POPSIZEMODEL_H
#define POPSIZEMODEL_H


#include "Random.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "Chronogram.h"
#include "OneMatrixPhyloProcess.h"
#include "ScatterMatrix.h"
#include "MGMixAAMutSelCodonSubMatrix.h"

#include "GCProcess.h"
#include "BranchMatrixPhyloProcess.h"

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, Var<RealVector>* inMu, VarArray<RealVector>* inMixture, Var<Profile>* inWeight, LengthTree* inpopsizetree) {
		SetWithRoot(true);
		nucmatrixtree = innucmatrixtree;
		mu = inMu;
		mixture = inMixture;
		weight = inWeight;
		popsizetree = inpopsizetree;
		statespace = instatespace;
		rootpopsize = new Const<PosReal>(1.0);
		RecursiveCreate(GetRoot());
	}

	Var<PosReal>* GetRootPopSize() {return rootpopsize;}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	protected:

	void specialUpdate(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->specialUpdate();
			specialUpdate(link->Out());
		}
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGMixAAMutSelCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),mu,mixture,weight,0,rootpopsize);
		}
		return new RandomMGMixAAMutSelCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),mu,mixture,weight,0,popsizetree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return popsizetree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* popsizetree;
	NucMatrixTree* nucmatrixtree;
	Const<PosReal>* rootpopsize;
	Var<RealVector>* mu;
	VarArray<RealVector>* mixture;
	Var<Profile>* weight;
};

class GCPopSizeAAMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	SequenceAlignment* codondata;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;
	Const<Real>* Zero;

	// chronogram 
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;
	
	// autocorrelated process
	Dvar<PosReal>* PriorSigma;
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;

	// variations of population size
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* popsizetree;

	// variations of population size
	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* gctree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;
	
	IIDNormal* MeanAALogFitness;
	int P; // number of degrees of freedom
	Dirichlet* weight;
	NormalIIDArray* iidarray;

	MatrixTree* codonmatrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GCPopSizeAAMutSelModel(string datafile, string treefile, int inP=20, bool sample=true, GeneticCodeType type=Universal)	{
		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------
		
		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		// a tree with all branch lengths iid from an exponential distribution of mean meanlength
		// meanlength is itself endowed with an exponential prior of mean 0.1
		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu); 
		mu->ClampAt(1);
		chronogram = new Chronogram(tree,mu);

		// a log normal process on that tree for the variations of the mutation rate
		PriorSigma = new Const<PosReal>(1);		
		sigma = new Gamma(One,PriorSigma);
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma);
		
		// another log normal process for the variations of pop size
		PriorTau = new Const<PosReal>(1);		
		tau = new Gamma(One,PriorTau);
		popsizetree = new LogNormalTreeProcess(chronogram,tau);
		
		// another log normal process for the variations of gc
		PriorTheta = new Const<PosReal>(1);		
		theta = new Gamma(One,PriorTheta);
		gctree = new LogNormalTreeProcess(chronogram,theta);
		
		// should clamp at the root
		popsizetree->GetRootRate()->ClampAt(0);

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stattree = new GCStatTree(gctree);
		nucmatrixtree = new NucMatrixTree(relrate,stattree);

		MeanAALogFitness = new IIDNormal(Naa,Zero,One);
		MeanAALogFitness->SetAtZero();
		MeanAALogFitness->Clamp();

		cerr << "iidarray\n";
		P = inP;
		weight = new Dirichlet(P);
		iidarray = new NormalIIDArray(Naa,P,Zero,One);
		for (int p=0; p<P; p++)	{
			for (int k=0; k<Naa; k++)	{
				(*iidarray->GetVal(p))[k] = 0.1 * Random::Uniform();
			}
		}

		cerr << "codon matrix\n";
		cerr << ((CodonStateSpace*) codondata->GetStateSpace())->GetNstate() << '\n';
		codonmatrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(),nucmatrixtree,MeanAALogFitness,iidarray,weight,popsizetree);

		cerr << "process\n";
		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(lognormaltree, codonmatrixtree, codondata);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorTau);
		RootRegister(PriorSigma);
		RootRegister(PriorTheta);
		// is this important?
		// RootRegister(chronogram->GetRootAge());
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(popsizetree->GetRootRate());
		RootRegister(codonmatrixtree->GetRootPopSize());
		RootRegister(relrate);
		RootRegister(weight);
		Register();

		MakeScheduler();

		if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
			cerr << '\n';
			GetTree()->Print(cerr);
			cerr << '\n';
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GCPopSizeAAMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetPopSizeTree() {return popsizetree;}
	LogNormalTreeProcess* GetGCTree() {return gctree;}
	LengthTree* GetChronogram() {return chronogram;}
	
	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += mu->GetLogProb();
		total += chronogram->GetLogProb();

		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();

		total += tau->GetLogProb();
		total += popsizetree->GetLogProb();

		total += theta->GetLogProb();
		total += gctree->GetLogProb();

		total += relrate->GetLogProb();

		total += MeanAALogFitness->GetLogProb();
		total += iidarray->GetLogProb();
		total += weight->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(mu,1),10,"mu");
		scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		scheduler.Register(new SimpleMove(chronogram,1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),30,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),30,"chrono");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),30,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),30,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(popsizetree,1),30,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.1),30,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.01),30,"popsize");

		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),100,"theta");
		scheduler.Register(new SimpleMove(gctree,1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.01),30,"gc");

		scheduler.Register(new ProfileMove(relrate,1,1),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.3,1),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.1,2),30,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,3),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.03),30,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),30,"relrates");

		scheduler.Register(new SimpleMove(MeanAALogFitness,1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.01),10,"mean");

		scheduler.Register(new ProfileMove(weight,1,1),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.3,2),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.3,1),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.1,2),20,"weight");
		scheduler.Register(new ProfileMove(weight,0.1,4),20,"weight");
		scheduler.Register(new SimpleMove(weight,0.1),20,"weight");
		scheduler.Register(new SimpleMove(weight,0.01),20,"weight");

		scheduler.Register(new SimpleMove(iidarray,1),4,"iidarray/sigma");
		scheduler.Register(new SimpleMove(iidarray,0.1),4,"iidarray/sigma");
		scheduler.Register(new SimpleMove(iidarray,0.01),4,"iidarray/sigma");
	
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	// Draw a sample from the prior
	void drawSample()	{
		cerr << "sample\n";
		mu->Sample();
		// chronogram->Sample();
		sigma->Sample();
		sigma->setval(10);
		lognormaltree->Sample();
		tau->Sample();
		tau->setval(10);
		popsizetree->Sample();
		// set all log popsizes at 0
		popsizetree->Reset();
		popsizetree->Clamp();
		theta->Sample();
		theta->setval(10);
		gctree->Sample();
		gctree->Reset();
		cerr << "relrate\n";
		relrate->Sample();
		weight->Sample();
		// MeanAALogFitness->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		codonmatrixtree->specialUpdate();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanLogPopSize()	{
		return popsizetree->GetMeanLogRate();
	}
	
	double GetVarLogPopSize()	{
		return popsizetree->GetVarLogRate();
	}

	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}
	
	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tgc\tvar\tlength\tlogpopsize\tvar\tgrandmean\tgrandvar\tweightent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << stattree->GetMeanGCContent();
		os << '\t' << stattree->GetVarGCContent();
		os << '\t' << GetMeanLogPopSize();
		os << '\t' << GetVarLogPopSize();
		os << '\t' << iidarray->GetGrandMean();
		os << '\t' << iidarray->GetGrandVar();
		os << '\t' << weight->GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *tau << '\n';
		os << *popsizetree << '\n';
		os << *theta << '\n';
		os << *gctree << '\n';
		os << *relrate << '\n';
		os << *MeanAALogFitness << '\n';
		os << *weight << '\n';
		os << *iidarray << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *tau;
		is >> *popsizetree;
		is >> *theta;
		is >> *gctree;
		is >> *relrate;
		is >> *MeanAALogFitness;
		is >> *weight;
		is >> *iidarray;
	}

	void GetAASelectionCoefficients(double** tmps)	{
		for (int i=0; i<Naa; i++)	{
			for (int j=0; j<Naa; j++)	{
				tmps[i][j] = 0;
			}
		}

		for (int p=0; p<P; p++)	{
			for (int i=0; i<Naa; i++)	{
				for (int j=0; j<Naa; j++)	{
					if (i!=j)	{
						tmps[i][j] += (*weight)[p] * ((*iidarray->GetVal(p))[j] - (*iidarray->GetVal(p))[i]);
					}
				}
			}
		}

		for (int i=0; i<Naa; i++)	{
			for (int j=0; j<Naa; j++)	{
				tmps[i][j] /= P;
			}
		}
	}
};

#endif

