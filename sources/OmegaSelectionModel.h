
#ifndef SELECTIONGTR_H
#define SELECTIONGTR_H

#include <omp.h> 

#include <stdio.h>
// #include <gsl/gsl_sf_bessel.h>
// #include <gsl/gsl_sf_psi.h>
 
#include "ProbModel.h"
#include "OneMatrixPhyloProcess.h"
#include "Move.h"
#include "GTRSubMatrix.h"
#include "CodonSubMatrix.h"
#include "CodonStateSpace.h"
#include "CodonSequenceAlignment.h"

#include "IID.h"
#include "IIDGamma.h"
#include "SelectionPhyloProcess.h"
#include "ProfileConjugatePath.h"
#include "Jeffreys.h"
#include "ConjugateInverseWishart.h"
#include "ValArray.h"
#include "ContinuousData.h"

class OmegaSelectionModel : public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	FileSequenceAlignment* data;
	TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;
    ContinuousData* contdata;

	// number of sites 
	int Nsite;

	// number of states (4 for nucleic acids, 20 for amino-acids, 61 for codons)
	int Nstate;

	// number of categories
	int K; 

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* One;
	Const<PosReal>* Ten;

	Exponential* lambda;
	GammaTree* gamtree;

	Dirichlet* relrate;
	Dirichlet* nucstationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;


//different omega for all categories K
	Const<PosReal>* PriorShape;
	Const<PosReal>* PriorScale;

	ExpArray* shape;
	ExpArray* scale;

	GammaArray** omega;

	RandomSubMatrix*** submatrix;

	AllocationTree* allocatetree;

	ProfilePathConjugateArray* patharray;
	SelectionMatrixTree* matrixtree;
	PhyloProcess* phyloprocess;

	int conjugate;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	OmegaSelectionModel(string datafile, string contdatafile, string treefile, int inK, bool sample = true) {

		conjugate =1;
		K = inK;
		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data,true );
        contdata = new FileContinuousData(contdatafile);

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate(); // # states (61 for codons)

		int Nnuc = 4;
//		int Naa = 20;

		cerr << Nsite << '\t' << Nstate << '\n';

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------

		One = new Const<PosReal>(1);
		Ten = new Const<PosReal>(10);
		
		// a tree with all branch lengths iid from an exponential distribution of rate lambda 
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		lambda = new Exponential(Ten,Exponential::MEAN);
		lambda->setval(50);
		gamtree = new GammaTree(tree,One,lambda);

        int offset = 0;
        if (K == 2) {
            offset = -1;
        }
		allocatetree = new AllocationTree(tree, contdata, K, offset);

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		nucstationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate, nucstationary, true);

		PriorShape = new Const<PosReal>(1);
		PriorScale = new Const<PosReal>(1);

		shape = new ExpArray(K,PriorShape, Exponential::MEAN);
		scale = new ExpArray(K,PriorScale, Exponential::MEAN);


		omega = new GammaArray* [K];

		for (int k=0; k<K; k++)	{
			omega[k] = new GammaArray(Nsite,shape->GetExpVal(k),scale->GetExpVal(k));

		}

		cerr << "submatrices\n";

		submatrix = new RandomSubMatrix** [K];
		for (int k=0; k<K; k++)	{
			submatrix[k] = new RandomSubMatrix* [Nsite];
			for (int i=0; i<Nsite;i++){
				submatrix[k][i] = new RandomMGOmegaCodonSubMatrix((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omega[k]->GetGammaVal(i),false);
			}
		}



		if (conjugate)	{
			cerr << "conjugate\n";
			matrixtree = 0;
			patharray = new ProfilePathConjugateArray(Nsite,K,submatrix);
			patharray->InactivateSufficientStatistic();
			phyloprocess = new ProfileConjugateSelectionPhyloProcess(allocatetree,gamtree,codondata,patharray);
		}
		else	{
			matrixtree = new SelectionMatrixTree(allocatetree, submatrix);
			cerr << "create phylo process\n";
			phyloprocess = new SelectionPhyloProcess(gamtree,0,matrixtree,codondata);
		}

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "register\n";
		RootRegister(One);
		RootRegister(Ten);
		RootRegister(relrate);
		RootRegister(nucstationary);
		RootRegister(PriorShape);
		RootRegister(PriorScale);

		Register();

		if (sample)	{
			cerr << "initialise\n";
			// Sample();
			cerr << "sample completed\n";
			Update();
			cerr << "update completed\n" ;
			cerr << "ln L " << GetLogLikelihood() << '\n';
			// cerr << "random calls " << rnd::GetRandom().GetCount() << '\n';
			Trace(cerr);
		}

		cerr << "trace completed\n";

		cerr << "scheduler\n";
		MakeScheduler();

		cerr << "model created\n";
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~OmegaSelectionModel() {}

	Tree* GetTree() {return tree;}


        double Update(bool check = false)       {
                cerr << "update with phyloprocess\n";
                double ret = ProbModel::Update();
                phyloprocess->Move(1);
                ret = ProbModel::Update();
                return ret;
        }

// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();
		total += relrate->GetLogProb();
		total += nucstationary->GetLogProb();

		for(int k=0;k<K;k++){
				total += shape->GetVal(k)->GetLogProb();
				total += scale->GetVal(k)->GetLogProb();
				total += omega[k]->GetLogProb();
		
		}
		return total;
	}

	double GetTreeLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();
		return total;
	}

	double GetNucLogPrior()	{
		double total = 0;
		total += relrate->GetLogProb();
		total += nucstationary->GetLogProb();
		return total;
	}


	double GetOmegaHyperLogPrior()	{
		double total = 0;
		for(int k=0;k<K;k++){
			total += shape->GetVal(k)->GetLogProb();
			total += scale->GetVal(k)->GetLogProb();
		
		}
		return total;
	}

	double GetOmegaLogPrior()	{
		double total = 0;
		for(int k=0;k<K;k++){
			total += omega[k]->GetLogProb();
		
		}
		return total;
	}


	double GetLogLikelihood()	{
		// return 0;
		return phyloprocess->GetLogProb();
	}

	void MakeScheduler()	{

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");

		for( int n=0;n<3;n++){
			// first register the moves on the global variables
			scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");

			scheduler.Register(new SimpleMove(gamtree,3),4,"branch lengths");
			scheduler.Register(new SimpleMove(gamtree,2.3),4,"branch lengths");

			if (conjugate)	{
				scheduler.Register(new ProfileConjugateMove(patharray,1),1,"activate suff stat");
			}

			for (int m=0; m<20; m++)	{

				scheduler.Register(new ProfileMove(relrate,1,1),1,"relrates");
				scheduler.Register(new ProfileMove(relrate,0.3,3),1,"relrates");
				scheduler.Register(new ProfileMove(relrate,0.1,3),1,"relrates");

				scheduler.Register(new ProfileMove(nucstationary,0.1,1),1,"nucstationary");
				scheduler.Register(new ProfileMove(nucstationary,0.03,1),1,"nucstationary");

				for (int omrep=0; omrep<10; omrep++)	{
				for(int k=0;k<K;k++){


					std::stringstream temp;
					std::string indice;
					temp << k;
					temp >> indice;

					scheduler.Register(new GammaArrayMove(omega[k],1),1,"omega" + indice);

					scheduler.Register(new ExpArrayMove(shape,1),10,"shape");
					scheduler.Register(new ExpArrayMove(shape,0.3),10,"shape");

					scheduler.Register(new ExpArrayMove(scale,1),10,"scale");
					scheduler.Register(new ExpArrayMove(scale,0.3),10,"scale");

				}
				}
			}

			if (conjugate)	{
				scheduler.Register(new ProfileConjugateMove(patharray,0),1,"inactivate suff stat");
			}
		}
	}


	double Move(double tuning)	{
		scheduler.Cycle(1,1,false,false);
		return 1;
	}


	void drawSample()	{
		cerr << "in sample\n";
		exit(1);
	}

	int GetSite()	{
		return Nsite;
	}

	int GetaaState()	{
		return Naa;
	}

	int GetCodonState()	{
		return Nstate;
	}

	int GetCategory()	{
		return K;
	}


	double GetLength()	{
		return gamtree->GetTotalLength(); 
	}

	double GetOmega(int cat,int site)	{
		return (*omega[cat]->GetVal(site));
	}

	GTRRandomSubMatrixWithNormRates* Getnucmatrix()	{
		 return nucmatrix;
	}

	RandomSubMatrix* Getsubmatrix(int cat, int site)	{
		int k = cat;
		int i = site;
		return (submatrix[k][i]);
	}


	CodonStateSpace* GetCodonStateSpace()	{
		return codondata->GetCodonStateSpace();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\ttreeprior\tnucprior\tomegaprior\thyperprior";
		os << "\tlnL\tlength\t";
		for(int k=0;k<K;k++){
			os << "shape" << k << '\t';
			os << "scale" << k << '\t';
			os << "omega" << k << '\t';
			os << "possel" << k << '\t';
		}
		os << "rrent\n";
	}

	int GetAboveOne(int k)	{
		int tot = 0;
		for (int i=0; i<Nsite; i++)	{
			if (omega[k]->GetVal(i)->val() > 1)	{
				tot++;
			}
		}
		return tot;
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior();
		os << '\t' << GetTreeLogPrior();
		os << '\t' << GetNucLogPrior();
		os << '\t' << GetOmegaLogPrior();
		os << '\t' << GetOmegaHyperLogPrior();
		os << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();

		for( int k=0;k<K;k++){
			os << '\t' << *shape->GetVal(k);
			os << '\t' << *scale->GetVal(k);
			os << '\t' << omega[k]->GetMean();
			os << '\t' << GetAboveOne(k);
		}
		os << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os.precision(7);
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *relrate << '\n';
		os << *nucstationary << '\n';
		for( int k=0;k<K;k++){
			os << *shape << '\n';
			os << *scale << '\n';
			os << *omega[k] << '\n';
		}
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *relrate;
		is >> *nucstationary;
		for( int k=0;k<K;k++){
			is >> *shape;
			is >> *scale;
			is >> *omega[k];
		}
	}

	void PostPredAli(string name)	{

		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample(true);
		phyloprocess->GetLeafData(datacopy);
		ostringstream s;
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
	}
	
};

#endif
