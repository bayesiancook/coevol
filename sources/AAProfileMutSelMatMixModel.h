
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
#include "MGAAProfileMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"
#include "AAProfileMutSelMixtureRandomMatrix.h"

class AAProfileMutSelMatrixFiniteMixture : public MatrixFiniteMixture<Profile>	{

	public:

	AAProfileMutSelMatrixFiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, Var<PosReal>* inNeff, string inprofilefile) :
			MatrixFiniteMixture<Profile>(insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		Neff = inNeff;
		profilefile = inprofilefile;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};

	void SetProfileValues()	{
		if (profilefile == "none")	{
			cerr << "no profile file given\n";
			exit(1);
		}
		else	{
			ifstream is((profilefile).c_str());
			if (!is)	{
				cerr << "Profile file not found\n";
				exit(1);
			}
			for (int k=0; k<GetComponentNumber(); k++)	{
				if (is.eof())	{
					cerr << "Profile file ends too soon.  Does the initial number of components given match the number in the file?\n";
					exit(1);
				}
				is >> *(GetDirichlet(k));
				GetDirichlet(k)->Clamp();  // This will prevent any moves on the profiles... maybe change this...
				//cerr << *(GetDirichlet(k)) << "\n";
			}
			string temp;
			is >> temp;
			if (!is.eof())	{
				cerr << "Are there more profiles given than the number of components specified?";
				cerr << " ...using only first " << GetComponentNumber() << " to preset mixture.\n";
			}
		}
	}


	protected:

	Dirichlet* GetDirichlet(int k)    {
		Dirichlet* temp = dynamic_cast<Dirichlet*>(GetComponent(k));
		if (!temp)	{
			cerr << "null pointer...\n";
			exit(1);
		}
                return temp;
	}

	MixtureRandomMatrix<Profile>* CreateComponent(int k)	{
		return new AAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<PosReal>* Neff;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
	string profilefile;
};

class AAProfileMutSelMatMixValMove : public MCUpdate	{

	public:

	AAProfileMutSelMatMixValMove(AAProfileMutSelMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	AAProfileMutSelMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};


class AAProfileMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	CodonSequenceAlignment* datacopy;
	ContinuousData* contdata;
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

	/*
	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	// autocorrelated process
	Dvar<PosReal>* PriorSigma;
	Gamma* sigma;
	LogNormalTreeProcess* lognormaltree;
	*/

	// tree and branch lengths
	Dvar<PosReal>* PriorLambda;
	Exponential* lambda;
	Dvar<PosReal>* mu;
	GammaTree* gamtree;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;

	// effective population size
	Const<PosReal>* neff;

	GTRRandomSubMatrixWithNormRates* nucmatrix;

	int P; // number of degrees of freedom
	Dirichlet* AAProfileCenter;
        Exponential* AAProfileConcentration;
        Const<PosReal>* AAProfileConcentrationPrior;

	public :

	AAProfileMutSelMatrixFiniteMixture* aaprofilemutselmix;
	MatrixMixturePhyloProcess<Profile>* phyloprocess;

	string profilefile;
	// constructor
	// this is where the entire graph structure of the model is created

	AAProfileMutSelModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal, string inprofilefile="none")	{
		// fetch data from file
		//
		cout << "in constructor\n";
		cout.flush();

		cout << "datafile: " << datafile << "\n";
		cout.flush();

		nucdata = new FileSequenceAlignment(datafile);
		cout << "One\n";
		cout.flush();
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		datacopy = 0;
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (e.g., 20 for amino acids)

		profilefile = inprofilefile;

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
		/*
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
		sigma->setval(10);
		cout << "\n" << "sigma: " << *sigma << "\n";

		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		lognormaltree->Reset();
		*/

		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,mu,lambda);
		//*
		//gamtree->SetBranchLengths();
		cout << "\n" << "length: " << GetLength() << "\n\n";

		// another log normal process for the variations of pop size

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		neff = new Const<PosReal>(1);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';
		AAProfileCenter = new Dirichlet(Naa);
		AAProfileCenter->setuniform();
                AAProfileConcentrationPrior = new Const<PosReal>((double)(Naa));
                AAProfileConcentration = new Exponential(AAProfileConcentrationPrior, Exponential::MEAN);
		AAProfileConcentration->setval(20);
		aaprofilemutselmix = new AAProfileMutSelMatrixFiniteMixture(Nsite,P,AAProfileCenter,AAProfileConcentration,codonstatespace,nucmatrix, neff, profilefile);
		if (profilefile != "none")	{
			aaprofilemutselmix->SetProfileValues();
			AAProfileCenter->Clamp();
			AAProfileConcentration->Clamp();
		}

		// update before phyloprocess
		Update();

		cerr << "create phylo process\n";
		//phyloprocess = new MatrixMixturePhyloProcess<Profile>(lognormaltree,aaprofilemutselmix,codondata);
		phyloprocess = new MatrixMixturePhyloProcess<Profile>(gamtree,aaprofilemutselmix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		phyloprocess->Sample();

		cerr << "register\n";
		/*
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		*/
		RootRegister(PriorLambda);
		RootRegister(mu);
		RootRegister(neff);
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(AAProfileCenter);
		RootRegister(AAProfileConcentrationPrior);
		RootRegister(aaprofilemutselmix->GetWeightVector());
		Register();

		MakeScheduler();

		if (sample)	{
			//cerr << "sample model\n";
			//Sample();

			cerr << "update\n";
			Update();

			//TraceHeader(cerr);
			//Trace(cerr);
			//cerr << '\n';
			//GetTree()->Print(cerr);
			//cerr << '\n';
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~AAProfileMutSelModel() {}

	Tree* GetTree() {return tree;}
	//LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	//LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return  Ntaxa;}
	int GetNstate() {return Nstate;}
	CodonSequenceAlignment* GetCodonData() {return codondata;}

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

		/*
		total += mu->GetLogProb();
		total += chronogram->GetLogProb();

		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();
		*/

		total += lambda->GetLogProb();
		total += gamtree->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();

		total += aaprofilemutselmix->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	double GetSiteLogLikelihood(int i)	{
		double ret = phyloprocess->SiteLogLikelihood(i);
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		//scheduler.Register(new SimpleMove(mu,1),10,"mu");
		//scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		/*
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
		*/
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aaprofilemutselmix,1),1,"aaprofilemutsel mix weight alloc");
		scheduler.Register(new MatMixValMove<Profile>(aaprofilemutselmix,0.3,10),1,"aaprofilemutsel mix simple val");
		scheduler.Register(new MatMixValMove<Profile>(aaprofilemutselmix,0.1,10),1,"aaprofilemutsel mix simple val");
		scheduler.Register(new MatMixValMove<Profile>(aaprofilemutselmix,0.01,10),1,"aaprofilemutsel mix simple val");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aaprofilemutselmix,1),1,"aaprofilemutsel mix weight alloc");
		scheduler.Register(new AAProfileMutSelMatMixValMove(aaprofilemutselmix,1,1,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelMatMixValMove(aaprofilemutselmix,1,2,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelMatMixValMove(aaprofilemutselmix,0.3,3,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aaprofilemutselmix,1),1,"aaprofilemutsel mix weight alloc");
		scheduler.Register(new AAProfileMutSelMatMixValMove(aaprofilemutselmix,0.3,4,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelMatMixValMove(aaprofilemutselmix,0.1,5,1),10,"aaprofilemutsel mix subset");
		scheduler.Register(new MatMixWeightAllocMove<Profile>(aaprofilemutselmix,1),1,"aaprofilemutsel mix weight alloc");

		scheduler.Register(new ProfileMove(AAProfileCenter,0.05,5),10,"aaprofile center profile move");
		scheduler.Register(new ProfileMove(AAProfileCenter,0.06,5),10,"aaprofile center profile move");
		scheduler.Register(new ProfileMove(AAProfileCenter,0.07,5),10,"aaprofile center profile move");
		scheduler.Register(new ProfileMove(AAProfileCenter,0.08,5),10,"aaprofile center profile move");
		scheduler.Register(new ProfileMove(AAProfileCenter,0.09,5),10,"aaprofile center profile move");
		scheduler.Register(new ProfileMove(AAProfileCenter,0.1,5),10,"aaprofile center profile move");

		scheduler.Register(new SimpleMove(AAProfileConcentration,0.1),10,"aaprofile concentration");
		scheduler.Register(new SimpleMove(AAProfileConcentration,0.2),10,"aaprofile concentration");
		scheduler.Register(new SimpleMove(AAProfileConcentration,0.3),10,"aaprofile concentration");
		scheduler.Register(new SimpleMove(AAProfileConcentration,0.4),10,"aaprofile concentration");
		scheduler.Register(new SimpleMove(AAProfileConcentration,0.5),10,"aaprofile concentration");
		scheduler.Register(new SimpleMove(AAProfileConcentration,1.0),10,"aaprofile concentration");
		/*
		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");
		*/

		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"branch lengths");

		scheduler.Register(new ProfileMove(relrate,1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.3,2),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.1,3),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.03),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}

	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample...\n";
		/*
		mu->Sample();
		// chronogram->Sample();
		sigma->Sample();
		sigma->setval(10);
		lognormaltree->Sample();
		*/
		cerr << "stat\n";
		stationary->Sample();
		cerr << "relrate\n";
		relrate->Sample();
		AAProfileConcentration->Sample();
		AAProfileCenter->Sample();
		aaprofilemutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";
	}


	// various summary statistics
	// used to check mcmc convergence

	/*
	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}

	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}
	*/

	double GetLength()	{
		//return lognormaltree->GetTotalLength();
		return gamtree->GetTotalLength();
	}

	//double GetGrandMeanLogFitness()	{
	//	double mean = 0;
	//	for (int i=0; i<Nsite; i++)	{
	//		mean += (*aaprofilemutselmix)[i]->GetMean();
	//	}
	//	mean /= Nsite;
	//	return mean;
	//}

	//double GetGrandVarLogFitness()	{
	//	double mean = 0;
	//	for (int i=0; i<Nsite; i++)	{
	//		mean += (*aaprofilemutselmix)[i]->GetVar();
	//	}
	//	mean /= Nsite;
	//	return mean;
	//}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		//os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tstatent\trrent\n";
		os << "#lnpri\tlnL\t\tncomp\teffsize\tcenterent\tcenter\tconcentration\tlength\tstatent\tstat\trrent\trr\n";
		//Trace(os);
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	//
	//
	// ***CHECK THESE***
	//
	//
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << aaprofilemutselmix->GetComponentNumber() << '\t';
		os << aaprofilemutselmix->GetEffSize() << '\t';
		os << AAProfileCenter->val().GetEntropy() << '\t';
		os << *AAProfileCenter << '\t';
		os << *AAProfileConcentration << '\t';
		os << GetLength() << '\t';
		os << stationary->val().GetEntropy() << '\t';
		os << *stationary << '\t';
		os << relrate->val().GetEntropy() << '\t';
		os << *relrate << '\t';
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		/*
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		*/

		os << *lambda<< '\n';
		os << *gamtree << '\n';

		os << *relrate << '\n';
		os << *stationary << '\n';
		os << *AAProfileCenter << '\n';
		os << *AAProfileConcentration << '\n';
		os << *aaprofilemutselmix << '\n';
	}

	void FromStreamWithUpdate(istream& is)	{
		/*
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		*/

		is >> *lambda;
		is >> *gamtree;

		is >> *relrate;
		is >> *stationary;
		is >> *AAProfileCenter;
		is >> *AAProfileConcentration;
		is >> *aaprofilemutselmix;
		Update();
		phyloprocess->Sample();
	}
	void FromStream(istream& is)	{
		/*
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		*/

		is >> *lambda;
		is >> *gamtree;

		is >> *relrate;
		is >> *stationary;
		is >> *AAProfileCenter;
		is >> *AAProfileConcentration;
		is >> *aaprofilemutselmix;

		//Update();
		//phyloprocess->Sample();

	}

	void MakeDataCopy()	{
		datacopy = new CodonSequenceAlignment(codondata);
	}

	void SamplePosteriorMapping()	{
		phyloprocess->Sample();
	}

	void PostPredSample()	{
		phyloprocess->PostPredSample();
		if (! datacopy)	{
			MakeDataCopy();
		}
		phyloprocess->GetLeafData(datacopy);
	}

	CodonSequenceAlignment* GetPostPredCodonData() {return datacopy;}

};

