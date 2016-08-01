
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
#include "GTRSubMatrix.h"
#include "MGAAProfileMutSelCodonSubMatrix.h"

#include "BranchMatrixMixture.h"
//#include "GCProcess.h"
#include "PrecisionNormalTreeProcess.h"
#include "AAProfileMutSelMixtureRandomMatrixTree.h"


class AAProfileMutSelBranchMatrixFiniteMixture : public BranchMatrixFiniteMixture<Profile>	{

	public:

	AAProfileMutSelBranchMatrixFiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, BranchVarTree<PosReal>* innefftree, Var<PosReal>* inrootneff, string inprofilefile) :
			BranchMatrixFiniteMixture<Profile>(innefftree->GetTree(), insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		nefftree = innefftree;
		rootneff = inrootneff;
		profilefile = inprofilefile;
		Create();
	}

	double MoveValues(double tuning, int n) {
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)      {
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	}

	void SetProfileValues()	{
		if (profilefile == "none")	{
			cerr << "no profile file given\n";
			exit(1);
		}
		else	{
			cout << "reading profiles\n";
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

	Tree* GetTree() {return nefftree->GetTree();}

	protected:

	Dirichlet* GetDirichlet(int k)    {
		Dirichlet* temp = dynamic_cast<Dirichlet*>(GetComponent(k));
		if (!temp)      {
			cerr << "null pointer...\n";
			exit(1);
		}
		return temp;
	}

	MixtureRandomMatrixTree<Profile>* CreateComponent(int k)	{
		return new AAProfileMutSelMixtureRandomMatrixTree(center,concentration,statespace, nucmatrix, nefftree, rootneff);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
	BranchVarTree<PosReal>* nefftree;
	Var<PosReal>* rootneff;
	string profilefile;
};


class AAProfileMutSelBranchMatMixValMove : public MCUpdate	{

	public:

	AAProfileMutSelBranchMatMixValMove(AAProfileMutSelBranchMatrixFiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	AAProfileMutSelBranchMatrixFiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};

class NeffAAProfileMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	//ContinuousData* contdata;
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

	// autocorrelated process
	Dvar<PosReal>* PriorTheta;
	Gamma* theta;
	LogNormalTreeProcess* nefftree;
	Const<PosReal>* rootneff;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	//GCStatTree* stattree;
	//NucMatrixTree* nucmatrixtree;
	GTRRandomSubMatrixWithNormRates* nucmatrix;


	int P; // number of degrees of freedom of mixture
	Dirichlet* AAProfileCenter;
        Exponential* AAProfileConcentration;
        Const<PosReal>* AAProfileConcentrationPrior;

	public :

	AAProfileMutSelBranchMatrixFiniteMixture* aaprofilemutselmix;
	BranchMatrixMixturePhyloProcess<Profile>* phyloprocess;

	string profilefile;
	// constructor
	// this is where the entire graph structure of the model is created

	NeffAAProfileMutSelModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal, string inprofilefile="none")	{

		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

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

		// another log normal process for the variations of neff
		PriorTheta = new Const<PosReal>(1);
		theta = new Gamma(One,PriorTheta);
		nefftree = new LogNormalTreeProcess(chronogram,theta,MEAN);
		nefftree->GetRootRate()->ClampAt(0); // NOTE:  clamp the root of nefftree at 0...
		nefftree->Reset();
		nefftree->specialUpdate();

		rootneff = new Const<PosReal>(1);

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);
		//relrate->Sample();
		//stationary->Sample();

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "mixture\n";
		P = inP;
		cerr << P << '\n';

		AAProfileCenter = new Dirichlet(Naa);
		AAProfileCenter->setuniform();
		//AAProfileCenter->Clamp();
		AAProfileConcentrationPrior = new Const<PosReal>((double)(Naa));
		AAProfileConcentration = new Exponential(AAProfileConcentrationPrior, Exponential::MEAN);
		AAProfileConcentration->setval(20);
		//AAProfileConcentration->ClampAt(20);
		aaprofilemutselmix = new AAProfileMutSelBranchMatrixFiniteMixture(Nsite,P,AAProfileCenter,AAProfileConcentration,codonstatespace,nucmatrix, nefftree, rootneff, profilefile);
		if (profilefile != "none")	{
			aaprofilemutselmix->SetProfileValues();
			AAProfileCenter->Clamp();
			AAProfileConcentration->Clamp();
		}

		if (P==Nsite)	{
			aaprofilemutselmix->SetHalpernBrunoAllocations();
		}
		//AAProfileConcentration->Sample();
		//AAProfileCenter->Sample();
		// aaprofilemutselmix->Sample();

		// update before phyloprocess

		Update();

		cerr << "create phylo process\n";
		phyloprocess = new BranchMatrixMixturePhyloProcess<Profile>(lognormaltree,aaprofilemutselmix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "sample\n";
		phyloprocess->Sample();

		cerr << "register\n";
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorTheta);
		RootRegister(rootneff);
		RootRegister(PriorMu);
		RootRegister(PriorSigma);
		RootRegister(lognormaltree->GetRootRate());
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(nefftree->GetRootRate()); // check this...
		RootRegister(AAProfileCenter);
		RootRegister(AAProfileConcentrationPrior);
		RootRegister(aaprofilemutselmix->GetWeightVector());
		Register();

		MakeScheduler();

		if (sample)	{
			//cerr << "sample model\n";
			//Sample();

			//cerr << "update\n";
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

	~NeffAAProfileMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetNeffTree() {return nefftree;}
	LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return Ntaxa;}


	//double Move(double tuning = 1)	{
	//	scheduler.Cycle(1,1,true,true);
	//	return 1;
	//}


	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += mu->GetLogProb();
		total += chronogram->GetLogProb();

		total += sigma->GetLogProb();
		total += lognormaltree->GetLogProb();

		total += theta->GetLogProb();
		total += nefftree->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();

		total += aaprofilemutselmix->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{

		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		//scheduler.Register(new SimpleMove(mu,1),10,"mu");
		//scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		//scheduler.Register(new BranchMatMixWeightAllocMove<Profile>(aaprofilemutselmix,1),1,"aaprofilemutsel mix weight alloc");
		scheduler.Register(new BranchMatMixValMove<Profile>(aaprofilemutselmix,0.3,10),1,"aaprofilemutsel mix simple val");
		scheduler.Register(new BranchMatMixValMove<Profile>(aaprofilemutselmix,0.1,10),1,"aaprofilemutsel mix simple val");
		scheduler.Register(new AAProfileMutSelBranchMatMixValMove(aaprofilemutselmix,1,1,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatMixValMove(aaprofilemutselmix,1,2,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatMixValMove(aaprofilemutselmix,0.3,3,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatMixValMove(aaprofilemutselmix,0.3,4,10),1,"aaprofilemutsel mix subset");

		/*
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
		*/


		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");

		scheduler.Register(new SimpleMove(theta,1),100,"theta");
		scheduler.Register(new SimpleMove(theta,0.1),100,"theta");
		scheduler.Register(new SimpleMove(nefftree,1),10,"neff");
		scheduler.Register(new SimpleMove(nefftree,0.1),10,"neff");
		scheduler.Register(new SimpleMove(nefftree,0.01),10,"neff");

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
		scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.003),10,"relrates");

		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
		scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

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
		cerr << "stat\n";
		stationary->Sample();
		cerr << "relrate\n";
		relrate->Sample();
		cerr << "theta\n";
		theta->Sample();
		theta->setval(10);
		cerr << "nefftree\n";
		nefftree->Sample();
		nefftree->Reset();
		cerr << "aaprofilecenter\n";
		AAProfileCenter->Sample();
		cerr << "aaprofileconcentration\n";
		AAProfileConcentration->Sample();
		cerr << "aaprofilemutselmix\n";
		aaprofilemutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";

		cerr << "drawSample() called...exiting\n";
		exit(1);
	}


	// various summary statistics
	// used to check mcmc convergence

	double GetMeanRho()	{
		return lognormaltree->GetMeanRate();
	}

	double GetVarRho()	{
		return lognormaltree->GetVarRate();
	}

	double GetLength()	{
		return lognormaltree->GetTotalLength();
	}

	//double GetGrandMeanLogFitness()	{
	//	double mean = 0;
	//	for (int i=0; i<Nsite; i++)	{
	//		mean += (*aamutselmix)[i]->GetMean();
	//	}
	//	mean /= Nsite;
	//	return mean;
	//}

	//double GetGrandVarLogFitness()	{
	//	double mean = 0;
	//	for (int i=0; i<Nsite; i++)	{
	//		mean += (*aamutselmix)[i]->GetVar();
	//	}
	//	mean /= Nsite;
	//	return mean;
	//}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#lnpri\tlnL\t\tncomp\teffsize\tcenterent\tcenter\tconcentration\tlength\tstatent\tstat\trrent\trr\n";
		//os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tgc\tvar\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
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
		os << *mu << '\n';
		os << *chronogram << '\n';
		os << *sigma << '\n';
		os << *lognormaltree << '\n';
		os << *theta << '\n';
		os << *nefftree << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
		os << *AAProfileCenter << '\n';
		os << *AAProfileConcentration << '\n';
		os << *aaprofilemutselmix << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		is >> *sigma;
		is >> *lognormaltree;
		is >> *theta;
		is >> *nefftree;
		is >> *relrate;
		is >> *stationary;
		is >> *AAProfileCenter;
		is >> *AAProfileConcentration;
		is >> *aaprofilemutselmix;
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


