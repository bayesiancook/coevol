
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
#include "MSCodonSubMatrix.h"

// #include "IIDArray.h"
#include "IID.h"
#include "SelectionPhyloProcess.h"
#include "ProfileConjugatePath.h"
#include "Jeffreys.h"
#include "ConjugateInverseWishart.h"


class ComplexDirichletIIDArrayMove : public MCUpdate	{

	public: 

	ComplexDirichletIIDArrayMove(DirichletIIDArray* inselectarray, double intuning, int innrep) : selectarray(inselectarray), tuning(intuning), nrep(innrep)	{}

	double Move(double tuning_modulator)	{
		

		double* tot = new double[selectarray->GetSize()];
		for (int i=0; i<selectarray->GetSize(); i++)	{
			tot[i] = 0;
		}
		
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif

		for (int i=0; i<selectarray->GetSize(); i++)	{
			for (int rep=0; rep<nrep; rep++)	{
				tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*1,2);
				// tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.3,2);
				// tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*1,5);
				tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.3,5);
				tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.1,10);
				// tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.03,10);
				// tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.03,10);
			}
		}

		double total = 0;
		for (int i=0; i<selectarray->GetSize(); i++)	{
			total+= tot[i];
		}
		delete[] tot;

		return total / selectarray->GetSize() / nrep / 6;
	}	

	private:

	DirichletIIDArray* selectarray;
	double tuning;
	int nrep;
};


class DirichletCodonUsageSelectionModelMS : public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	FileSequenceAlignment* data;
	TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

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

//same center and concentration for all categories K
//	Dirichlet* center;
//	Const<PosReal>* PriorConcentration;
//	Exponential* concentration;

//different center and concentration for all categories K
	Dirichlet** center;
	Const<PosReal>** PriorConcentration;
	Exponential** concentration;

	DirichletIIDArray** selectionprofile;
	
	RandomSubMatrix*** submatrix;

	AllocationTree* allocatetree;

	ProfilePathConjugateArray* patharray;
	SelectionMatrixTree* matrixtree;
	PhyloProcess* phyloprocess;

	Dirichlet* codonusageselection;

	int conjugate;
	string type;
	string mechanism;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	DirichletCodonUsageSelectionModelMS(string datafile, string treefile, int inK, string intype, int inconjugate, string inmechanism, bool sample = true) {

		conjugate = inconjugate;
		type = intype;
		mechanism = inmechanism;
		K = inK;

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data,true );

		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate(); // # states (61 for codons)

		int Nnuc = 4;
		int Naa = 20;

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

		allocatetree = new AllocationTree(tree, K);

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		nucstationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate, nucstationary, true);

		codonusageselection = new Dirichlet(Nstate);
		codonusageselection->setuniform();
		codonusageselection->Clamp();


		center = new Dirichlet* [K];
		PriorConcentration = new Const<PosReal>* [K];
		concentration = new Exponential* [K];

		for (int k=0; k<K; k++)	{

				center[k] = new Dirichlet(Naa);
				PriorConcentration[k] = new Const<PosReal>(Naa);
				concentration[k] = new Exponential(PriorConcentration[k],Exponential::MEAN);

			if(type == "clamp_MCMC")	{
				center[k]->setuniform();
				center[k]->Clamp();

				concentration[k]->setval(20);
				concentration[k]->Clamp();
			}
		}


		selectionprofile = new DirichletIIDArray*[K];
		for (int k=0; k<K; k++)	{
			selectionprofile[k] = new DirichletIIDArray(Nsite,center[k],concentration[k]);
			// selectionprofile[k]->SetUniform();
		}

		cerr << "submatrices\n";

//Square Root //phenimenological
		if(mechanism == "SR")	{
			submatrix = new RandomSubMatrix** [K];
			for (int k=0; k<K; k++)	{
				submatrix[k] = new RandomSubMatrix* [Nsite];
				for (int i=0; i<Nsite;i++){
					submatrix[k][i] = new RandomMGSRFitnessCodonUsageSubMatrix((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, selectionprofile[k]->GetVal(i), codonusageselection,false);
				}
			}
		}


//Mutation Selection // mechanistic
		if(mechanism == "MS")	{
			submatrix = new RandomSubMatrix** [K];
			for (int k=0; k<K; k++)	{
				submatrix[k] = new RandomSubMatrix* [Nsite];
				for (int i=0; i<Nsite;i++){
					submatrix[k][i] = new RandomMGMSFitnessCodonUsageSubMatrix((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, selectionprofile[k]->GetVal(i), codonusageselection,false);
				}
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
		RootRegister(codonusageselection);
		for (int i=0; i<K; i++)	{
			RootRegister(PriorConcentration[i]);
			RootRegister(center[i]);
		}
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

	~DirichletCodonUsageSelectionModelMS() {}

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

		for(int i=0;i<K;i++){
			total += concentration[i]->GetLogProb();
			total += center[i]->GetLogProb();

			for (int j=0; j<Nsite;j++){
				total += selectionprofile[i]->GetVal(j)->GetLogProb();
			}
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

				scheduler.Register(new ProfileMove(relrate,0.1,1),1,"relrates");
				scheduler.Register(new ProfileMove(relrate,0.03,3),1,"relrates");
				scheduler.Register(new ProfileMove(relrate,0.01,3),1,"relrates");

				scheduler.Register(new ProfileMove(nucstationary,0.1,1),1,"nucstationary");
				scheduler.Register(new ProfileMove(nucstationary,0.01,1),1,"nucstationary");

				scheduler.Register(new ProfileMove(codonusageselection,0.3,1),1,"codonusageselection");
				scheduler.Register(new ProfileMove(codonusageselection,0.1,3),1,"codonusageselection");
				scheduler.Register(new SimpleMove(codonusageselection,0.01),1,"codonusageselection");

				for(int i=0;i<K;i++){


					std::stringstream temp;
					std::string indice;
					temp << i;
					temp >> indice;


					scheduler.Register(new ProfileMove(center[i],0.1,1),10,"center");
					scheduler.Register(new ProfileMove(center[i],0.01,5),10,"center");

					scheduler.Register(new SimpleMove(concentration[i],0.1),10,"conc");
					scheduler.Register(new SimpleMove(concentration[i],0.03),10,"conc");

					scheduler.Register(new ComplexDirichletIIDArrayMove(selectionprofile[i],0.15,10),1,"stat" + indice);
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

	double GetSelectionProfile(int cat,int site, int state)	{
		return (*selectionprofile[cat]->GetVal(site))[state];
	}

	double GetStationary(int cat,int site, int state)	{
		return GetSelectionProfile(cat,site,state);
	}

	double GetLength()	{
		return gamtree->GetTotalLength(); 
	}

	double GetCenterEntropy(int k)	{

		double tot = 0;
		for (int i=0; i<center[k]->GetDim(); i++)	{
			double tmp = (*center[k])[i];
			tot -= tmp * log (tmp);
		}
		return tot;
	}

	GTRRandomSubMatrixWithNormRates* Getnucmatrix()	{
		 return nucmatrix;
	}

	RandomSubMatrix* Getsubmatrix(int cat, int site)	{
		int k = cat;
		int i = site;
		return (submatrix[k][i]);
	}

	double GetCodonUsage(int state)	{
		int i = state;
		return (*codonusageselection)[i];
	}

	CodonStateSpace* GetCodonStateSpace()	{
		return codondata->GetCodonStateSpace();
	}

	// creates the header of the <model_name>.trace file 
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\t";
		for(int i=0;i<K;i++){
			os << "center" << i << '\t';
			os << "conc" << i << '\t';
			os << "selentropy" << i << '\t';
		}
		os << "rrent\n";
		// os << "rrent\tcount\tmax\n";
	}


	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior();
		os << '\t' << GetLogLikelihood();
		os << '\t' << GetLength();
		for( int i=0;i<K;i++){
			os << '\t' << GetCenterEntropy(i);
			os << '\t' << concentration[i]->val();
			os << '\t' << selectionprofile[i]->GetMeanEntropy();
		}
		os << '\t' << relrate->val().GetEntropy();
        /*
		os << '\t' << normerrorcount;
		os << '\t' << normerrormax;
        */
		os << '\n';
		os.flush();
	}
	
	void ToStream(ostream& os)	{
		os.precision(7);
		os << *lambda << '\n';
		os << *gamtree << '\n';
		os << *relrate << '\n';
		os << *nucstationary << '\n';
		os << *codonusageselection << '\n';
		for( int i=0;i<K;i++){
			os << *center[i] << '\n';
			os << *concentration[i] << '\n';
			os << *selectionprofile[i] << '\n';
		}
	}

	void FromStream(istream& is)	{
		is >> *lambda;
		is >> *gamtree;
		is >> *relrate;
		is >> *nucstationary;
		is >> *codonusageselection;
		for( int i=0;i<K;i++){
			is >> *center[i];
			is >> *concentration[i];
			is >> *selectionprofile[i];
		}
	}

	void PostPredAli(string name)	{

		CodonSequenceAlignment* datacopy = new CodonSequenceAlignment(codondata);
		phyloprocess->PostPredSample(false);
		phyloprocess->GetLeafData(datacopy);
		ostringstream s;
		ofstream os((name + ".ali").c_str());
		datacopy->ToStream(os);
	}

};

#endif
