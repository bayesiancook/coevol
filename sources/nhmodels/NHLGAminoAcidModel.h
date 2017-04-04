
#include "Random.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "BranchProcess.h"
#include "GTRSubMatrix.h"
#include "PrecisionNormalTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "StatProcess.h"
// #include "ACGTProcess.h"
#include "OneMatrixPhyloProcess.h"

class LGMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	LGMatrixTree(BranchVarTree<Profile>* instattree) {
		SetWithRoot(true);
		stattree = instattree;
		RecursiveCreate(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new LGRandomSubMatrix(stattree->GetBranchVal(link->GetBranch()),true);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;

};

class BranchMatrixRASPhyloProcess : public PhyloProcess	{

	
	protected:

	public:

	BranchMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
		rate = inrate;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), GetRate(site), matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	

	Var<PosReal>* GetRate(int site)	{
		if (! rate)	{
			return 0;
		}
		return rate->GetVal(site);
	}

	BranchValPtrTree<RandomSubMatrix>* matrixtree;
	VarArray<PosReal>* rate;
};


class NHAminoAcidModel: public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	ProteinSequenceAlignment* proteindata;
	ContinuousData* contdata;
	GCContinuousData* gcdata;
	TaxonSet* taxonset;
	CodonStateSpace* codonstatespace;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ntaxa;

	public:
	int Ncont;
	int K;
	int Kaa;
	int Kgc;
	int Kcont;

	private:
	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;

	// tree and branch lengths
	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	GammaTree* gamtree;
	
	Gamma* alpha;
	GammaIIDArray* rate;

	// parameters for the multivariate tree process describing variations of ATGC
	Gamma** GammaDiag;
	Rvar<PosReal>** diag;
	RvarVec* priorOnSigmaZero;
	SigmaZero* sigmaZero;
	// InverseWishartMatrix* Sigma;
	ConjugateInverseWishart* Sigma;
	
	// NodeMultiVariateTree
	// MultiVariateTreeProcess* process;
	ConjugateMultiVariateTreeProcess* process;

	Dirichlet* refstat;
	StatTree* stattree;
	BranchStatTree* branchstattree;
	LGMatrixTree* aamatrixtree;
	
	BranchMatrixRASPhyloProcess* phyloprocess;

	bool withsep;
	Gamma** theta;
	LogNormalTreeProcess** separateprocess;

	
	public :

	// constructor
	// this is where the entire graph structure of the model is created

	NHAminoAcidModel(string datafile, string treefile, string contdatafile, bool inwithsep, bool sample=true, GeneticCodeType type=Universal)	{

		withsep = inwithsep;

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

		cerr << "protein data\n";
		proteindata = new ProteinSequenceAlignment(codondata);

		cerr << "gc data\n";
		gcdata = new GCContinuousData(codondata,2);

		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------
		
		One = new Const<PosReal>(1);

		PriorLambda = new Const<PosReal>(10);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,One,lambda);
		lambda->setval(10);
		gamtree->Sample();

		alpha = new Gamma(One,One);
		rate = new GammaIIDArray(Nsite,alpha,alpha);
		alpha->Sample();
		rate->Sample();


		// compute number of components of the multivariate process
		// aa : + 20
		// gc : + 1
		// contdata : + Ncont
		Kcont = 0;
		Kaa = 0;
		Kgc = Naa;
		K = Naa + 1;
		if (contdata)	{
			Kcont = Naa + 1;
			K += Ncont;
		}

		cerr << "K : " << K << '\t' << Kaa << '\t' << Kgc << '\t' << Kcont << '\n';

		GammaDiag = new Gamma*[K];
		diag = new Rvar<PosReal>*[K];
		for (int k=0; k<K; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(1);
		}
		priorOnSigmaZero = new RvarVec(diag,K);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		// Sigma = new InverseWishartMatrix(sigmaZero,K + 3);
		Sigma = new ConjugateInverseWishart(sigmaZero,K + 3);
		Sigma->Sample();

		// process = new MultiVariateTreeProcess(Sigma,gamtree);
		process = new ConjugateMultiVariateTreeProcess(Sigma,gamtree);
		process->Sample();
		process->Reset();
		for (int i=0; i<Naa; i++)	{
			process->GetMultiNormal(process->GetRoot())->ClampAt(0,i);
		}

		refstat = new Dirichlet(Naa);
		refstat->Sample();
		refstat->setuniform();
		cerr << "stattree\n";
		stattree = new StatTree(refstat,process,Kaa);
		cerr << "branch stat tree\n";
		branchstattree = new BranchStatTree(stattree);
		// in fact, this is a aa matrix tree
		aamatrixtree = new LGMatrixTree(branchstattree);

		// here, true means : igc is logit transformed...
		process->SetAndClamp(gcdata,Kgc,0,true);

		// ... whereas condata are log transformed
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,Kcont+i,i);
			}
		}

		if (withsep)	{
			theta = new Gamma*[Ncont+1];
			separateprocess = new LogNormalTreeProcess*[Ncont+1];

			for (int i=0; i<Ncont+1; i++)	{
				theta[i] = new Gamma(One,One);
				separateprocess[i] = new LogNormalTreeProcess(gamtree,theta[i],MEAN);
				if (!i)	{
					separateprocess[i]->ClampAt(gcdata,0,true);
				}
				else	{
					separateprocess[i]->ClampAt(contdata,i-1);
				}
			}
		}
				
		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixRASPhyloProcess(gamtree, rate, aamatrixtree, proteindata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(refstat);
		RootRegister(PriorLambda);
		// RootRegister(process->GetMultiNormal(process->GetRoot()));
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

	~NHAminoAcidModel() {}

	Tree* GetTree() {return tree;}
	BranchStatTree* GetBranchStatTree()	{return branchstattree;}
	StatTree* GetStatTree()	{return stattree;}
	LengthTree* GetGamTree() {return gamtree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	CovMatrix* GetCovMatrix() {return Sigma;}
	
	int GetNtaxa() {return Ntaxa;}
	int GetNstate() {return Nstate;}

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

		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();

		total += alpha->GetLogProb();
		total += rate->GetLogProb();

		for (int k=0; k<K; k++)	{
			total += GammaDiag[k]->GetLogProb();
		}
		total += Sigma->GetLogProb();
		total += process->GetLogProb();
	
		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				total += theta[i]->GetLogProb();
				total += separateprocess[i]->GetLogProb();
			}
		}

		total += refstat->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),30,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.1),30,"gamtree");
		scheduler.Register(new SimpleMove(gamtree,0.01),30,"gamtree");

		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(rate,3),10,"rates across sites");
		scheduler.Register(new SimpleMove(rate,0.3),10,"rates across sites");

		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,30),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,1,30),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.1,30),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(Sigma,process,0.01,30),1,"conjugate sigma - process");

		scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,10,20,1,30),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,1,20,1,30),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,0.1,20,1,30),1,"conjugate sigma - process - piecewise translation");
		
		for (int i=0; i<Ncont; i++)	{
			scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,10,Kcont+i,1,30),1,"conjugate sigma - process - piecewise translation");
			scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,1,Kcont+i,1,30),1,"conjugate sigma - process - piecewise translation");
			scheduler.Register(new ConjugateMultiVariatePiecewiseTranslationMove(Sigma,process,0.1,Kcont+i,1,30),1,"conjugate sigma - process - piecewise translation");
		}
		/*
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.1,0,20,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.1,0,20,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.01,0,20,10),1,"conjugate sigma - process - piecewise translation");

		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,1,20,1,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.1,20,1,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.01,20,1,10),1,"conjugate sigma - process - piecewise translation");

		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,1,21,1,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.1,21,1,10),1,"conjugate sigma - process - piecewise translation");
		scheduler.Register(new ConjugateMultiVariateWholeTreePiecewiseTranslationMove(Sigma,process,0.01,21,1,10),1,"conjugate sigma - process - piecewise translation");
		*/

		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				scheduler.Register(new SimpleMove(theta[i],1),100,"theta");
				scheduler.Register(new SimpleMove(theta[i],0.1),100,"theta");
				scheduler.Register(new SimpleMove(separateprocess[i],10),10,"sep");
				scheduler.Register(new SimpleMove(separateprocess[i],1),10,"sep");
				scheduler.Register(new SimpleMove(separateprocess[i],0.1),10,"sep");
			}
		}

		scheduler.Register(new ProfileMove(refstat,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(refstat,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(refstat,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(refstat,0.001),10,"stat");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	
	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
		// lambda->Sample();
		cerr << "gamtree\n";
		// gamtree->Sample();

		alpha->Sample();
		alpha->setval(10);
		rate->Sample();

		cerr << "multivariate process\n";
		for (int k=0; k<K; k++)	{
			GammaDiag[k]->Sample();
			GammaDiag[k]->setval(1);
		}
		Sigma->Sample();
		Sigma->SetIdentity();
		process->Sample();
		// process->Reset();
		// process->SetEmpirical();

		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				theta[i]->Sample();
				separateprocess[i]->Sample();
			}
		}

		cerr << "phyloprocess\n";
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	/*
	double GetMeanGC()	{
	}
	
	double GetVarGC()	{
	}
	*/

	double GetLength()	{
		return gamtree->GetTotalLength();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\talpha\tmeanlogit(gc)";
		for (int k=0; k<Ncont; k++)	{
			os << "\tcont" << k;
		}
		os << "\tmeanaa\tvar\tstatent";
		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				os << "\tsep" << i;
			}
		}
		os << "\tsigma\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << *alpha;
		os << '\t' << process->GetMean(Kgc);
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << process->GetMean(Kcont + k);
		}
		os << '\t' << stattree->GetMeanEntropy();
		os << '\t' << stattree->GetVarEntropy();
		os << '\t' << refstat->GetEntropy();
		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				os << '\t' << separateprocess[i]->GetMeanRate();
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << (*Sigma)[Kgc][Kcont + k];
		}
		for (int k=0; k<Naa; k++)	{
			os << '\t' << (*Sigma)[k][Kgc];
		}
		for (int k=0; k<Ncont; k++)	{
			for (int i=0; i<Naa; i++)	{
				os << '\t' << (*Sigma)[i][Kcont+k];
			}
		}
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *alpha << '\n';
		os << *rate << '\n';
		for (int k=0; k<K; k++)	{
			os << *GammaDiag[k] << '\n';
		}
		os << *Sigma << '\n';
		os << *process << '\n';
		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				os << *theta[i] << '\n';
				os << *separateprocess[i] << '\n';
			}
		}
			
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *alpha;
		is >> *rate;
		for (int k=0; k<K; k++)	{
			is >> *GammaDiag[k];
		}
		is >> *Sigma;
		is >> *process;
		if (withsep)	{
			for (int i=0; i<Ncont+1; i++)	{
				is >>  *theta[i];
				is >> *separateprocess[i];
			}
		}
	}

	/*
	void MakeDataCopy()	{
		datacopy = new CodonSequenceAlignment(codondata);
	}

	void PostPredSample()	{
		phyloprocess->PostPredSample();
		if (! datacopy)	{
			MakeDataCopy();
		}
		phyloprocess->GetLeafData(datacopy);
	}

	CodonSequenceAlignment* GetPostPredCodonData() {return datacopy;}
	*/
};

