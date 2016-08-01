#ifndef GTRMODEL_H
#define GTRMODEL_H

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "GTRSubMatrix.h"
#include "OneMatrixPhyloProcess.h"
#include "Move.h"
#include "PhyloProcessMHMove.h"

class MHOneMatrixPhyloProcess : public OneMatrixPhyloProcess	{

	protected:

	public:

	MHOneMatrixPhyloProcess(LengthTree* intree, RandomSubMatrix* randmatrix,  SequenceAlignment* indata) :
		OneMatrixPhyloProcess(intree,randmatrix,indata), propmatrix(0) {
	}

	// EmpiricalSubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
	SubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
		if (! propmatrix)	{
			CreateProposalMatrix();
		}
		return propmatrix;
	}

	protected:

	void Unfold()	{
		PhyloProcess::Unfold();
		SetProposalMatrices();
	}

	void CreateProposalMatrix()	{
		cerr << "create proposal matrix\n";
		int nstate = GetMatrix()->GetNstate();
		rr = new double[nstate * (nstate-1) / 2];
		stat = new double[nstate];
		GetData()->GetEmpiricalFreq(stat);
		for (int i=0; i<nstate; i++)	{
			cerr << stat[i] << '\t';
		}
		cerr << '\n';
		if (nstate == Nnuc)	{
			// transition / transversion = 2
			rr[0] = 1;
			rr[1] = 2;
			rr[2] = 1;
			rr[3] = 1;
			rr[4] = 2;
			rr[5] = 1;
			propmatrix = new GTRSubMatrix(nstate,rr,stat,GetMatrix()->isNormalised());
		}
		else if (nstate == Naa)	{
			propmatrix = new LGSubMatrix(stat,GetMatrix()->isNormalised());
		}
		else	{
			cerr << "error in MHOneMatrixPhyloProcess: unrecognized alphabet\n";
			exit(1);
		}
		propmatrix->CorruptMatrix();
		propmatrix->UpdateMatrix();
	}

	// void CreateProposalMatrix()	{
	// 	propmatrix = new EmpiricalGTRSubMatrix(GetMatrix()->GetNstate(),GetMatrix()->isNormalised());
	// 	RecursiveRegisterProposalMatrix(GetRoot());
	// 	propmatrix->RefreshStatistics();
	// }

	// void RecursiveRegisterProposalMatrix(const Link* from)	{
	// 	for (int i=0; i<GetNsite(); i++)	{
	// 		propmatrix->AddPath(GetPath(from->GetBranch(),i));
	// 	}
	// 	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
	// 		RecursiveRegisterProposalMatrix(link->Out());
	// 	}
	// }

	// EmpiricalSubMatrix* propmatrix;
	SubMatrix* propmatrix;
	double* rr;
	double* stat;
};

class MHOneMatrixRASPhyloProcess : public OneMatrixRASPhyloProcess	{

	protected:

	public:

	MHOneMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate, RandomSubMatrix* randmatrix,  SequenceAlignment* indata) :
		OneMatrixRASPhyloProcess(intree,inrate,randmatrix,indata)	{
	}

	// EmpiricalSubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
	SubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
		if (! propmatrix)	{
			CreateProposalMatrix();
		}
		return propmatrix;
	}

	protected:

	void Unfold()	{
		PhyloProcess::Unfold();
		SetProposalMatrices();
	}

	void CreateProposalMatrix()	{
		cerr << "create proposal matrix\n";
		int nstate = GetMatrix()->GetNstate();
		rr = new double[nstate * (nstate-1) / 2];
		stat = new double[nstate];
		GetData()->GetEmpiricalFreq(stat);
		for (int i=0; i<nstate; i++)	{
			cerr << stat[i] << '\t';
		}
		cerr << '\n';
		if (nstate == Nnuc)	{
			// transition / transversion = 2
			rr[0] = 1;
			rr[1] = 2;
			rr[2] = 1;
			rr[3] = 1;
			rr[4] = 2;
			rr[5] = 1;
			propmatrix = new GTRSubMatrix(nstate,rr,stat,GetMatrix()->isNormalised());
		}
		else if (nstate == Naa)	{
			propmatrix = new LGSubMatrix(stat,GetMatrix()->isNormalised());
		}
		else	{
			cerr << "error in MHOneMatrixPhyloProcess: unrecognized alphabet\n";
			exit(1);
		}
		propmatrix->CorruptMatrix();
		propmatrix->UpdateMatrix();
	}

	SubMatrix* propmatrix;
	double* rr;
	double* stat;
};


class GTRModel : public ProbModel	{

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	// tree and branch lengths
	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	GammaTree* gamtree;

	// relative exchange rates of the matrix
	Dvar<PosReal>* PriorMeanRelRate;
	IIDExp* relrate;

	// equilibrium frequencies of the matrix
	Dirichlet* stationary;

	// substitution matrix is relrate * stationary
	GTRRandomSubMatrix* matrix;

	bool mh;

	// rates across sites
	Dvar<PosReal>* PriorVarRate;
	Exponential* alpha;
	GammaIIDArray* rate;
	PhyloProcess* phyloprocess;
	// OneMatrixRASPhyloProcess* phyloprocess;
	// MHOneMatrixPhyloProcess* phyloprocess;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GTRModel(string datafile, string treefile, bool inmh) {
		// fetch data from file
		data = new FileSequenceAlignment(datafile);
		Nsite = data->GetNsite();	// # columns
		Nstate = data->GetNstate();	// # states (20 for amino acids)

		taxonset = data->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";

		// ----------
		// construction of the graph
		// ----------

		mh = inmh;
		// a tree with all branch lengths iid from an exponential distribution of rate lambda
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,mu,lambda);
		/*
		gamtree->SetBranchLengths();
		gamtree->Clamp();
		*/

		// a collection of Nsite rates, i.i.d. from a gamma distribution of mean 1 and variance 1/alpha
		// alpha is itself from an exponential prior of mean 1
		PriorVarRate = new Const<PosReal>(PosReal(1));
		alpha = new Exponential(PriorVarRate,Exponential::MEAN);
		rate = new GammaIIDArray(Nsite,alpha,alpha);


		//
		// a general time reversible (GTR) substitution matrix
		// for any pair of states l,m (distinct)
		//
		// Q_lm = rho_lm pi_m
		//
		// where rho_lm are the relative rates
		// and pi_m is the vector of equilibrium frequencies of the substitution process
		//
		// the proces is reversible, thus rho_lm = rho_ml, and there are Nstate (Nstate-1) / 2 distinct relative rates
		// all iid from an exponential of mean 1
		// the stationary probabilities (pi) are from a uniform distribution (Dirichlet of weights 1,1,1,1,...)
		relrate = new IIDExp(Nstate*(Nstate-1)/2);
		stationary = new Dirichlet(Nstate);
		matrix = new GTRRandomSubMatrix(relrate,stationary);

		/*
		alpha->ClampAt(1);
		for (int i=0; i<Nsite; i++)	{
			(*rate)[i]->ClampAt(1);
		}

		relrate->setall(1);
		relrate->Clamp();

		double* empfreq = new double[Nstate];
		data->GetEmpiricalFreq(empfreq);
		stationary->setarray(empfreq);
		stationary->Clamp();
		*/

		// a phylogenetic process
		// this is where everything comes together
		// the process is of type GTRRAS
		// which means that it gives the same GTR matrix for every site and every branch
		// but gives each site its own rate from the IID vector "rate"
		//
		cerr << "create phylo process\n";
		// phyloprocess = new OneMatrixPhyloProcess(gamtree,matrix,data);
		if (mh)	{
			phyloprocess = new MHOneMatrixRASPhyloProcess(gamtree,rate,matrix,data);
			// phyloprocess = new MHOneMatrixPhyloProcess(gamtree,matrix,data);
		}
		else	{
			phyloprocess = new OneMatrixRASPhyloProcess(gamtree,rate,matrix,data);
			// phyloprocess = new OneMatrixPhyloProcess(gamtree,matrix,data);
		}
		// phyloprocess = new OneMatrixRASPhyloProcess(gamtree,rate,matrix,data);
		// a phyloprocess needs to "unfold"
		// (recursive construction of all branch/site substitution processes)
		cerr << "unfold\n";
		cerr << "one matrix model\n";
		phyloprocess->Unfold();
		/*
		if (mh)	{
			phyloprocess->SetProposalMatrices();
		}
		*/
		cerr << "register\n";

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		RootRegister(PriorVarRate);
		RootRegister(PriorLambda);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		cerr << "initialise\n";
		Sample();
		cerr << "ok\n";
		Trace(cerr);

		cerr << "update\n";
		Update();
		cerr << "update ok\n";
		cerr << "scheduler\n";
		MakeScheduler();

		cerr << "model created\n";

		phyloprocess->SetName("phyloprocess");
		rate->SetName("rate");
		alpha->SetName("alpha");
		gamtree->SetName("bl");
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GTRModel() {}

	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += alpha->GetLogProb();
		total += rate->GetLogProb();
		total += gamtree->GetLogProb();
		total += stationary->GetLogProb();
		total += relrate->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void MakeScheduler()	{

		if (mh)	{
			scheduler.Register(new PhyloProcessMHMove(phyloprocess,500,-1,-1),1,"mapping");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,0.3),1,"mapping");
		}

		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(rate,1),10,"rate");
		scheduler.Register(new SimpleMove(rate,0.1),10,"rate");
		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"branch lengths");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new ProfileMove(stationary,1,2),10,"stat4");
		scheduler.Register(new SimpleMove(stationary,0.1),10,"stat");
		scheduler.Register(new SimpleMove(stationary,0.01),10,"stat");

		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,0.1),1,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,0.1),1,"lengthrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,relrate,1),1,"lengthrelrate");
		scheduler.Register(new MultiplicativeCompensatoryMove(gamtree,rate,1),1,"lengthrate");
	}

	double OldFashionedMove(double tuning)	{
		double total2 = 0;
		for (int i=0; i<10; i++)	{
			for (int j=0; j<10; j++)	{
				alpha->Move(tuning);
				rate->Move(tuning);
				lambda->Move(tuning);
				total2 += relrate->Move(tuning/10);
			}

			for (int k=0; k<2; k++)	{
				stationary->Move(tuning,1);
				stationary->Move(tuning/10,1);
				stationary->Move(tuning/50,2);
				stationary->Move(0.1*tuning);
				stationary->Move(0.01*tuning);
				stationary->Move(0.001*tuning);
			}
			gamtree->Move(tuning);
		}
		phyloprocess->Move(tuning);
		return 1;
	}

	// Draw a sample from the prior

	void drawSample()	{
		cerr << "alpha\n";
		alpha->Sample();
		cerr << "rate\n";
		rate->Sample();
		cerr << "ok\n";
		lambda->Sample();
		cerr << "gamtree\n";
		gamtree->Sample();
		stationary->Sample();
		relrate->Sample();
		cerr << "phyloprocess sample\n";
		phyloprocess->Sample();
		cerr << "ok\n";
	}

	// various summary statistics
	// used to check mcmc convergence

	double GetMeanRate()	{
		double total = 0;
		for (int i=0; i<Nsite; i++)	{
			total += (*rate)[i]->val();
		}
		total /= Nsite;
		return total;
	}

	double GetVarRate()	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<Nsite; i++)	{
			double tmp = (*rate)[i]->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= Nsite;
		var /= Nsite;
		var -= mean * mean;
		return var;
	}

	double GetLength()	{
		return gamtree->GetTotalLength() * matrix->GetRate() * GetMeanRate();
	}


	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tlambda\talpha\tmeanrate\tvarrate\tstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << GetLength();
		os << '\t' << lambda->val();
		os << '\t' << alpha->val();
		os << '\t' << GetMeanRate();
		os << '\t' << GetVarRate();
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
	}

	void ToStream(ostream& os)	{
		os << *alpha << '\n';
		os << *rate << '\n';
		os << *lambda<< '\n';
		os << *gamtree << '\n';
		os << *stationary << '\n';
		os << *relrate << '\n';
	}

	void FromStream(istream& is)	{
		is >> *alpha;
		is >> *rate;
		is >> *lambda;
		is >> *gamtree;
		is >> *stationary;
		is >> *relrate;
	}

	// various accessory functions to print tree with branch lengths
	// and make stochastic mappinga of substitution events along the tree
	void PrintTree(ostream& os)	{
		// gamtree->Print(os);
	}

	void PostMapping(ostream& os, bool red)	{
		phyloprocess->ResampleSub();
		for (int i=0; i<Nsite; i++)	{
			phyloprocess->GetSiteMapping(i)->Print(os, red);
		}
	}

	void PostMapping(int site, ostream& os, bool red)	{
		phyloprocess->ResampleSub(site); // Nielsen (recursive accept reject) clamped
		phyloprocess->GetSiteMapping(site)->Print(os, red);
	}

	void PostPredMapping(ostream& os, bool red)	{
		phyloprocess->PostPredSample(); // Nielsen (recursive accept reject) unclamped
		for (int i=0; i<Nsite; i++)	{
			phyloprocess->GetSiteMapping(i)->Print(os, red);
		}
	}

	void PostPredMapping(int site, ostream& os, bool red)	{
		phyloprocess->PostPredSample(site);
		phyloprocess->GetSiteMapping(site)->Print(os, red);
	}

};

#endif // GTRMODEL_H
