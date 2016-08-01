
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
#include "OneMatrixPhyloProcess.h"
#include "MGMixAAMutSelCodonSubMatrix.h"
#include "ConjugateMultiVariateTreeProcess.h"

#include "ACGTProcess.h"

class CodonMatrixTree : public BranchValPtrTree<RandomMGMixAAMutSelCodonSubMatrix>	{


	public:

	CodonMatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, Var<RealVector>* inMu, VarArray<RealVector>* inMixture, Var<Profile>* inWeight, LengthTree* inpopsizetree) {
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

	Tree* GetTree() {return popsizetree->GetTree();}

	protected:

	void specialUpdate(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->specialUpdate();
			specialUpdate(link->Out());
		}
	}

	RandomMGMixAAMutSelCodonSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGMixAAMutSelCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),mu,mixture,weight,0,rootpopsize);
		}
		return new RandomMGMixAAMutSelCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),mu,mixture,weight,0,popsizetree->GetBranchVal(link->GetBranch()));
	}

	private:

	CodonStateSpace* statespace;
	LengthTree* popsizetree;
	NucMatrixTree* nucmatrixtree;
	Const<PosReal>* rootpopsize;
	Var<RealVector>* mu;
	VarArray<RealVector>* mixture;
	Var<Profile>* weight;
};

class AAMatrixTree : public BranchValPtrTree<RandomAminoAcidReducedCodonSubMatrix>	{


	public:

	AAMatrixTree(CodonMatrixTree* incodonmatrixtree)	{
		SetWithRoot(true);
		codonmatrixtree = incodonmatrixtree;
		RecursiveCreate(GetRoot());
	}

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

	RandomAminoAcidReducedCodonSubMatrix* CreateBranchVal(const Link* link)	{
		return new RandomAminoAcidReducedCodonSubMatrix(codonmatrixtree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return codonmatrixtree->GetTree();}

	private:

	CodonMatrixTree* codonmatrixtree;
};

class BranchMatrixPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchMatrixPhyloProcess(LengthTree* intree, BranchValPtrTree<RandomAminoAcidReducedCodonSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:
	BranchValPtrTree<RandomAminoAcidReducedCodonSubMatrix>* matrixtree;
};

class ACGTPopSizeAAReducedMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	// amino acid data
	SequenceAlignment* data;
	TaxonSet* taxonset;
	CodonStateSpace* codonstatespace;

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

	// variations of population size
	Dvar<PosReal>* PriorTau;
	Gamma* tau;
	LogNormalTreeProcess* popsizetree;

	// parameters for the multivariate tree process describing variations of ATGC
	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* Sigma;
	
	// NodeMultiVariateTree
	ConjugateMultiVariateTreeProcess* process;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	StatTree* stattree;
	BranchStatTree* branchstattree;
	NucMatrixTree* nucmatrixtree;
	
	IIDNormal* MeanAALogFitness;
	int P; // number of degrees of freedom
	Dirichlet* weight;
	NormalIIDArray* iidarray;

	CodonMatrixTree* codonmatrixtree;
	AAMatrixTree* aamatrixtree;

	// phylo process
	BranchMatrixPhyloProcess* phyloprocess;
	
	public:

	// constructor
	// this is where the entire graph structure of the model is created

	ACGTPopSizeAAReducedMutSelModel(string datafile, string treefile, int inP=20, bool sample=true, GeneticCodeType type=Universal)	{
		// fetch data from file
		data = new FileSequenceAlignment(datafile);
		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)
		if (Nstate != Naa)	{
			cerr << "error : model works for amino acids only\n";
			exit(1);
		}

		codonstatespace = new CodonStateSpace(type);

		taxonset = data->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		Ntaxa = taxonset->GetNtaxa();

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
		// should clamp at the root
		popsizetree->GetRootRate()->ClampAt(0);
		
		GammaDiag = new Gamma*[Nnuc];
		diag = new Rvar<PosReal>*[Nnuc];
		for (int k=0; k<Nnuc; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(1);
		}
		priorOnSigmaZero = new RvarVec(diag,Nnuc);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		Sigma = new ConjugateInverseWishart(sigmaZero,Nnuc + 3);

		process = new ConjugateMultiVariateTreeProcess(Sigma,chronogram);

		stattree = new StatTree(process,Nnuc,0);
		branchstattree = new BranchStatTree(stattree);
		
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		nucmatrixtree = new NucMatrixTree(relrate,branchstattree);

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

		codonmatrixtree = new CodonMatrixTree(codonstatespace,nucmatrixtree,MeanAALogFitness,iidarray,weight,popsizetree);
		aamatrixtree = new AAMatrixTree(codonmatrixtree);

		cerr << "process\n";
		// a phylogenetic process
		phyloprocess = new BranchMatrixPhyloProcess(lognormaltree, aamatrixtree, data);
		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorTau);
		RootRegister(PriorSigma);
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

	~ACGTPopSizeAAReducedMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetPopSizeTree() {return popsizetree;}
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

		for (int k=0; k<Nnuc; k++)	{
			total += GammaDiag[k]->GetLogProb();
		}
		total += Sigma->GetLogProb();
		total += process->GetLogProb();

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

		/*
		for (int k=0; k<Nnuc; k++)	{
			scheduler.Register(new SimpleMove(GammaDiag[k],10),10,"theta");
			scheduler.Register(new SimpleMove(GammaDiag[k],1),10,"theta");
			scheduler.Register(new SimpleMove(GammaDiag[k],0.1),10,"theta");
		}
		*/

		/*
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,10,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.01,10),1,"conjugate sigma - process");
		*/

		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,10,1),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,1),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.1,1),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.01,1),1,"conjugate sigma - process");

		/*
		scheduler.Register(new SimpleMove(Sigma,10),100,"sigma");
		scheduler.Register(new SimpleMove(Sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(Sigma,0.1),100,"sigma");

		scheduler.Register(new SimpleMove(process,10),40,"multinormal");
		scheduler.Register(new SimpleMove(process,1),40,"multinormal");
		scheduler.Register(new SimpleMove(process,0.1),40,"multinormal");
		scheduler.Register(new SimpleMove(process,0.01),40,"multinormal");
		*/

		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new SimpleMove(MeanAALogFitness,1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.1),10,"mean");
		scheduler.Register(new SimpleMove(MeanAALogFitness,0.01),10,"mean");

		scheduler.Register(new ProfileMove(weight,1,1),10,"weight");
		scheduler.Register(new ProfileMove(weight,1,2),10,"weight");
		scheduler.Register(new ProfileMove(weight,1,4),10,"weight");
		scheduler.Register(new SimpleMove(weight,0.1),10,"weight");
		scheduler.Register(new SimpleMove(weight,0.01),10,"weight");

		scheduler.Register(new SimpleMove(iidarray,1),4,"iidarray");
		scheduler.Register(new SimpleMove(iidarray,0.1),4,"iidarray");
		scheduler.Register(new SimpleMove(iidarray,0.01),4,"iidarray");
	
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
		cerr << "multivariate process\n";
		for (int k=0; k<Nnuc; k++)	{
			GammaDiag[k]->Sample();
			GammaDiag[k]->setval(1);
		}
		Sigma->Sample();
		Sigma->SetIdentity();
		process->Sample();
		process->Reset();
		cerr << "relrate\n";
		relrate->Sample();
		weight->Sample();
		// MeanAALogFitness->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		codonmatrixtree->specialUpdate();
		aamatrixtree->specialUpdate();
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
		os << "#logprior\tlnL\tgc\tvar\tlength\tlogpopsize\tvar\tgrandmean\tgrandvar\tweightent\tstatent\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << stattree->GetMeanEntropy();
		os << '\t' << stattree->GetVarEntropy();
		os << '\t' << GetMeanLogPopSize();
		os << '\t' << GetVarLogPopSize();
		os << '\t' << iidarray->GetGrandMean();
		os << '\t' << iidarray->GetGrandVar();
		os << '\t' << weight->GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		for (int k=0; k<Nnuc; k++)	{
			for (int l=k+1; l<Nnuc; l++)	{
				os << '\t' << (*Sigma)[k][l];
			}
		}
		for (int k=0; k<Nnuc; k++)	{
			os << '\t' << (*Sigma)[k][k];
		}
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
		for (int k=0; k<Nnuc; k++)	{
			os << *GammaDiag[k] << '\n';
		}
		os << *Sigma << '\n';
		os << *process << '\n';
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
		for (int k=0; k<Nnuc; k++)	{
			is >> *GammaDiag[k];
		}
		is >> *Sigma;
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

	double ObservedCompositionalHeterogeneity()	{

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
			}
		}

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

	double PostPredCompositionalHeterogeneity()	{
		phyloprocess->PostPredSample(); // Nielsen (recursive accept reject) unclamped
		return phyloprocess->CompositionalHeterogeneityIndex();
	}
};


