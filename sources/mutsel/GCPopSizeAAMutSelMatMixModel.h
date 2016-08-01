
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
#include "BranchProcess.h"
// #include "PhyloProcess.h"
#include "GTRSubMatrix.h"
#include "MGAAMutSelCodonSubMatrix.h"

#include "BranchMatrixMixture.h"
#include "GCProcess.h"
#include "PrecisionNormalTreeProcess.h"


class AAMutSelMixtureRandomMatrixTree : public MixtureRandomMatrixTree<RealVector>	{

	public:

	AAMutSelMixtureRandomMatrixTree(Var<Real>* inmean, Var<PosReal>* invar, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* inmatrixtree, LengthTree* inpopsizetree)	{
		mean = inmean;
		var = invar;
		statespace = instatespace;
		nucmatrixtree = inmatrixtree;
		popsizetree = inpopsizetree;
		CreateRandomVariable();
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return nucmatrixtree->GetTree();}

	Rvar<RealVector>* CreateRandomVariable()	{
		rvar = new IIDNormal(Naa,mean,var);
		rvar->SetName("rvar / iid normal");
		return rvar;
	}

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		RandomSubMatrix* matrix = 0;
		if (popsizetree)	{
			matrix = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrixtree->GetBranchVal(link->GetBranch()), rvar,0,0,popsizetree->GetBranchVal(link->GetBranch()));
		}
		else	{
			matrix = new RandomMGAAMutSelCodonSubMatrix(statespace, nucmatrixtree->GetBranchVal(link->GetBranch()), rvar,0,0,0);
		}
		matrix->SetName("matrix");
		return matrix;
	}
	
	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	LengthTree* popsizetree;
};


class AAMutSelBranchMatrixFiniteMixture : public BranchMatrixFiniteMixture<RealVector>	{

	public:

	AAMutSelBranchMatrixFiniteMixture(int insize, int incomponentnumber, Var<Real>* inmean, Var<PosReal>* invar, CodonStateSpace* instatespace, BranchValPtrTree<RandomSubMatrix>* innucmatrixtree, LengthTree* inpopsizetree) :
			BranchMatrixFiniteMixture<RealVector>(innucmatrixtree->GetTree(), insize, incomponentnumber)	{

		mean = inmean;
		var = invar;
		statespace = instatespace;
		nucmatrixtree = innucmatrixtree;
		popsizetree = inpopsizetree;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetIIDNormal(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	Tree* GetTree() {return nucmatrixtree->GetTree();}
	
	protected:

	IIDNormal* GetIIDNormal(int k)	{
		IIDNormal* tmp = dynamic_cast<IIDNormal*>(GetComponent(k));
		if (! tmp)	{
			cerr << "error in AAMutSelMatrixFiniteMixture: component is not iid normal\n";
			exit(1);
		}
		return tmp;
	}

	MixtureRandomMatrixTree<RealVector>* CreateComponent(int k)	{
		return new AAMutSelMixtureRandomMatrixTree(mean,var,statespace,nucmatrixtree,popsizetree);
	}

	private:
	Var<Real>* mean;
	Var<PosReal>* var;
	CodonStateSpace* statespace;
	BranchValPtrTree<RandomSubMatrix>* nucmatrixtree;
	LengthTree* popsizetree;
};


class AAMutSelBranchMatMixValMove : public MCUpdate	{

	public:

	AAMutSelBranchMatMixValMove(AAMutSelBranchMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}	
	
	protected:

	AAMutSelBranchMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
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
	CodonStateSpace* codonstatespace;

	bool withgc;
	bool withpopsize;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ntaxa;

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

	// pop size
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* popsizetree;

	// substitution matrix is relrate * stationary
	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* gctree;

	Dirichlet* relrate;
	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;
	
	int P; // number of degrees of freedom

	public :

	AAMutSelBranchMatrixFiniteMixture* aamutselmix;
	BranchMatrixMixturePhyloProcess<RealVector>* phyloprocess;

	// constructor
	// this is where the entire graph structure of the model is created

	GCPopSizeAAMutSelModel(string datafile, string treefile, int inP, bool inwithgc, bool inwithpopsize, bool sample=true, GeneticCodeType type=Universal)	{

		withgc = inwithgc;
		withpopsize = inwithpopsize;

		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();
		Ntaxa = taxonset->GetNtaxa();

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
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		
		// another log normal process for the variations of pop size
		PriorTau = new Const<PosReal>(1);		
		tau = new Gamma(One,PriorTau);
		popsizetree = new LogNormalTreeProcess(chronogram,tau,MEAN);
		popsizetree->GetRootRate()->ClampAt(0);
		
		// another log normal process for the variations of gc
		PriorTheta = new Const<PosReal>(1);		
		theta = new Gamma(One,PriorTheta);
		gctree = new LogNormalTreeProcess(chronogram,theta,MEAN);
		
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stattree = new GCStatTree(gctree);
		nucmatrixtree = new NucMatrixTree(relrate,stattree);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		if (withpopsize)	{
			aamutselmix = new AAMutSelBranchMatrixFiniteMixture(Nsite,P,Zero,One,codonstatespace,nucmatrixtree,popsizetree);
		}
		else	{
			aamutselmix = new AAMutSelBranchMatrixFiniteMixture(Nsite,P,Zero,One,codonstatespace,nucmatrixtree,0);
		}

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixMixturePhyloProcess<RealVector>(lognormaltree,aamutselmix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorTheta);
		RootRegister(PriorTau);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(aamutselmix->GetWeightVector());
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
	LogNormalTreeProcess* GetGCTree() {return gctree;}
	LogNormalTreeProcess* GetPopSizeTree() {return popsizetree;}
	LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return Ntaxa;}
	
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

		total += aamutselmix->GetLogProb();
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
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.3,10),1,"aamutsel mix simple val");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.1,10),1,"aamutsel mix simple val");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(popsizetree,1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.01),10,"popsize");

		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),100,"theta");
		scheduler.Register(new SimpleMove(gctree,1),10,"gc");
		scheduler.Register(new SimpleMove(gctree,0.1),10,"gc");
		scheduler.Register(new SimpleMove(gctree,0.01),10,"gc");

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.003),10,"relrates");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	/*
	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(mu,1),10,"mu");
		scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.3,10),1,"aamutsel mix simple val");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.1,10),1,"aamutsel mix simple val");
		scheduler.Register(new BranchMatMixValMove<RealVector>(aamutselmix,0.01,10),1,"aamutsel mix simple val");
		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,1,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,1,2,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,3,10),1,"aamutsel mix subset");
		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.3,4,10),1,"aamutsel mix subset");
		scheduler.Register(new AAMutSelBranchMatMixValMove(aamutselmix,0.1,5,1),10,"aamutsel mix subset");
		scheduler.Register(new BranchMatMixWeightAllocMove<RealVector>(aamutselmix,1),1,"aamutsel mix weight alloc");

		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(tau,1),100,"tau");
		scheduler.Register(new SimpleMove(tau,0.1),100,"tau");
		scheduler.Register(new SimpleMove(popsizetree,1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.1),10,"popsize");
		scheduler.Register(new SimpleMove(popsizetree,0.01),10,"popsize");

		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),100,"theta");
		scheduler.Register(new SimpleMove(gctree,1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.1),30,"gc");
		scheduler.Register(new SimpleMove(gctree,0.01),30,"gc");

		scheduler.Register(new ProfileMove(relrate,1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.3,2),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.1,3),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.03),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	*/
	
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
		cerr << "relrate\n";
		relrate->Sample();
		cerr << "gc\n";
		theta->Sample();
		theta->setval(10);
		gctree->Sample();
		gctree->Reset();

		aamutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}
	
	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}

	double GetMeanLogPopSize()	{
		return popsizetree->GetMeanLogRate();
	}
	
	double GetVarLogPopSize()	{
		return popsizetree->GetVarLogRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	double GetGrandMeanLogFitness()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*aamutselmix)[i]->GetMean();
		}
		mean /= Nsite;
		return mean;
	}

	double GetGrandVarLogFitness()	{
		double mean = 0;
		for (int i=0; i<Nsite; i++)	{
			mean += (*aamutselmix)[i]->GetVar();
		}
		mean /= Nsite;
		return mean;
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tgc\tvar\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << aamutselmix->GetEffSize();
		os << '\t' << GetGrandMeanLogFitness();
		os << '\t' << GetGrandVarLogFitness();
		// os << '\t' << aamutselmix->GetWeightEntropy();
		os << '\t' << stattree->GetMeanGCContent();
		os << '\t' << stattree->GetVarGCContent();
		os << '\t' << GetMeanLogPopSize();
		os << '\t' << GetVarLogPopSize();
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
		os << *aamutselmix << '\n';
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
		is >> *aamutselmix;
	}

	double ObservedCompositionalHeterogeneity(ostream& os)	{

		SequenceAlignment* data = codondata;

		double** taxfreq = new double*[Ntaxa];
		for (int j=0; j<Ntaxa; j++)	{
			taxfreq[j] = new double[Nstate];
			for (int k=0; k<Nstate; k++)	{
				taxfreq[j][k] = 0;
			} 
		}

		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Ntaxa; j++)	{
				int state = data->GetState(j,i);
				if (state != unknown)	{
					taxfreq[j][state]++;
				}
			}
		}
				
		// make global freqs out of tax-specific freqs
		double* globalfreq = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			globalfreq[k] = 0;
			for (int j=0; j<Ntaxa; j++)	{
				globalfreq[k] += taxfreq[j][k];
			}
		} 

		// normalise
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			total += globalfreq[k];
		}
		for (int k=0; k<Nstate; k++)	{
			globalfreq[k] /= total;
		}
		for (int j=0; j<Ntaxa; j++)	{
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += taxfreq[j][k];
			}
			for (int k=0; k<Nstate; k++)	{
				taxfreq[j][k] /= total;
				os << taxfreq[j][k] << '\t';
			}
			os << '\n';
		}
		os << '\n';

		// compute max distance
		double maxdist = 0;
		for (int j=0; j<Ntaxa; j++)	{
			double dist = 0;
			for (int k=0; k<Nstate; k++)	{
				double tmp = (taxfreq[j][k] - globalfreq[k]);
				dist += tmp * tmp;
			}
			if (maxdist < dist)	{
				maxdist = dist;
			}
		}
		return maxdist;
	}

	double PostPredCompositionalHeterogeneity(ostream& os)	{
		phyloprocess->PostPredSample(); // Nielsen (recursive accept reject) unclamped
		return phyloprocess->CompositionalHeterogeneityIndex(os);
	}

};

#endif

