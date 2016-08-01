
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


class AAProfileMutSelBranchMatrixInfiniteMixture : public BranchMatrixInfiniteMixture<Profile>	{

	public:

	AAProfileMutSelBranchMatrixInfiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, BranchVarTree<PosReal>* innefftree, Var<PosReal>* inrootneff) :
			BranchMatrixInfiniteMixture<Profile>(innefftree->GetTree(), insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		nefftree = innefftree;
		rootneff = inrootneff;
		Create();
	}

	double MoveValues(double tuning, int n) {
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)      {
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
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
};


class AAProfileMutSelBranchMatInfMixValMove : public MCUpdate	{

	public:

	AAProfileMutSelBranchMatInfMixValMove(AAProfileMutSelBranchMatrixInfiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	AAProfileMutSelBranchMatrixInfiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};



class CodonBranchSitePath : public BranchSitePath	{


};



class ProfileBranchInfiniteMixturePhyloProcess : public BranchMatrixInfiniteMixturePhyloProcess<Profile>	{

	public:

	ProfileBranchInfiniteMixturePhyloProcess(LengthTree* intree, BranchMatrixInfiniteMixture<Profile>* inmatmix,  CodonSequenceAlignment* indata) : BranchMatrixInfiniteMixturePhyloProcess<Profile>(intree,inmatmix,indata)	{data = indata;}

	double	MyStatistic(ostream& os)	{
		double temp = check();
		os << temp;
		return temp;
	}

	double check()	{
		return 1.0;
	}

	int NonsynSubCount()	{
		int count = 0;
		for (int site=0; site<GetNsite(); site++)	{
			count += BranchSiteNonsynSubCount(GetRoot(), site);
		}
		return count;
	}

	int BranchSiteNonsynSubCount(Link* from, int site)	{
		int count = 0;
		for (Link* blink=from->Next(); blink!=from; blink=blink->Next())	{
			Plink* plink = GetPath(blink->GetBranch(), site)->Init();
			while (plink != GetPath(blink->GetBranch(), site)->Last())	{
				if (! data->GetCodonStateSpace()->Synonymous(plink->GetState(), plink->Next()->GetState()))	{
					count++;
				}
				plink=plink->Next();
			}
			count += BranchSiteNonsynSubCount(blink->Out(), site);
		}
		return count;
	}

	private:
	CodonSequenceAlignment* data;

};




class NeffAAProfileInfMixMutSelModel : public ProbModel {

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	CodonSequenceAlignment* dataunknown;
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

	public :

	int P; // number of degrees of freedom
	Dirichlet* AAProfileCenter;
        Exponential* AAProfileConcentration;
        Const<PosReal>* AAProfileConcentrationPrior;

	//public :

	AAProfileMutSelBranchMatrixInfiniteMixture* aaprofilemutselmix;
	//BranchMatrixInfiniteMixturePhyloProcess<Profile>* phyloprocess;
	ProfileBranchInfiniteMixturePhyloProcess* phyloprocess;

	// constructor
	// this is where the entire graph structure of the model is created

	NeffAAProfileInfMixMutSelModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal)	{

		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		dataunknown = new CodonSequenceAlignment(nucdata, true, type);
		dataunknown->Unclamp();

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

		// another log normal process for the variations of neff
		PriorTheta = new Const<PosReal>(1);
		theta = new Gamma(One,PriorTheta);
		nefftree = new LogNormalTreeProcess(chronogram,theta,MEAN);
		nefftree->GetRootRate()->ClampAt(0); // NOTE:  clamp the root of nefftree at 0...
		//nefftree->GetNodeVal(nefftree->GetRoot())->ClampAt(0); // NOTE:  clamp the root of nefftree at 0...
		nefftree->Reset(); // sets all log Neff values to 0, or Neff values to 1,
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
		AAProfileConcentrationPrior = new Const<PosReal>((double)(Naa));
		AAProfileConcentration = new Exponential(AAProfileConcentrationPrior, Exponential::MEAN);
		AAProfileConcentration->setval(20);
		//AAProfileConcentration->ClampAt(20);
		aaprofilemutselmix = new AAProfileMutSelBranchMatrixInfiniteMixture(Nsite,P,AAProfileCenter,AAProfileConcentration,codonstatespace,nucmatrix, nefftree, rootneff);
		//AAProfileConcentration->Sample();
		//AAProfileCenter->Sample();
		// aaprofilemutselmix->Sample();

		cerr << "create phylo process\n";
		phyloprocess = new ProfileBranchInfiniteMixturePhyloProcess(lognormaltree,aaprofilemutselmix,codondata);

		cerr << "unfold\n";
		phyloprocess->Unfold();
		//phyloprocess->Sample();


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
		RootRegister(nefftree->GetRootRate());
		RootRegister(AAProfileCenter);
		RootRegister(AAProfileConcentrationPrior);
		//RootRegister(aaprofilemutselmix->GetWeightVector());
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
			phyloprocess->Sample();
			Update();
		}

		/*if (sample)	{
			cerr << "sample model\n";
			Sample();

			cerr << "update\n";
			Update();

			TraceHeader(cerr);
			Trace(cerr);
			cerr << '\n';
			GetTree()->Print(cerr);
			cerr << '\n';
		}*/
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~NeffAAProfileInfMixMutSelModel() {}

	Tree* GetTree() {return tree;}
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LogNormalTreeProcess* GetNeffTree() {return nefftree;}
	LengthTree* GetChronogram() {return chronogram;}
	int GetNtaxa() {return Ntaxa;}

	//Exponential* GetAAProfileConcentration() {return AAProfileConcentration;}


	//double Move(double tuning = 1)	{
	//	scheduler.Cycle(1,1,true,true);
	//	return 1;
	//}


	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		//cout << "Entering GetLogPrior\n";
		//cout.flush();
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
		//cout << "Done GetLogPrior\n";
		//cout.flush();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// scheduler is empty. instead, we use the old fashioned move function (below)
	void MakeScheduler()	{

		//scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

		scheduler.Register(new BranchMatInfMixAllocMove<Profile>(aaprofilemutselmix,5,1),1,"aaprofilemutsel mix alloc");
		scheduler.Register(new AAProfileMutSelBranchMatInfMixValMove(aaprofilemutselmix,1,1,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatInfMixValMove(aaprofilemutselmix,0.5,2,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatInfMixValMove(aaprofilemutselmix,0.1,3,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new AAProfileMutSelBranchMatInfMixValMove(aaprofilemutselmix,0.01,4,10),1,"aaprofilemutsel mix subset");
		scheduler.Register(new BranchMatInfMixAlphaMove<Profile>(aaprofilemutselmix,0.1,1),1,"aaprofilemutsel alpha");
		scheduler.Register(new BranchMatInfMixAlphaMove<Profile>(aaprofilemutselmix,0.2,1),1,"aaprofilemutsel alpha");
		scheduler.Register(new BranchMatInfMixAlphaMove<Profile>(aaprofilemutselmix,0.3,1),1,"aaprofilemutsel alpha");
		scheduler.Register(new BranchMatInfMixAlphaMove<Profile>(aaprofilemutselmix,0.4,1),1,"aaprofilemutsel alpha");
		scheduler.Register(new BranchMatInfMixAlphaMove<Profile>(aaprofilemutselmix,0.5,1),1,"aaprofilemutsel alpha");
		scheduler.Register(new BranchMatInfMixAllocMove<Profile>(aaprofilemutselmix,5,1),1,"aaprofilemutsel mix alloc");


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

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates profile move");
		scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates profile move");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates simple move");
		scheduler.Register(new SimpleMove(relrate,0.003),10,"relrates simple move");

		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
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
		AAProfileCenter->setuniform();
		cerr << "aaprofileconcentration\n";
		AAProfileConcentration->Sample();
		cerr << "aaprofilemutselmix\n";
		aaprofilemutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		phyloprocess->Sample();
		cerr << "ok\n";

		cerr << "drawSample() called... exiting...\n";
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
		//os << "#logprior\tlnL\tlength\teffsize\tgrandmean\tgrandvar\tgc\tvar\trrent\n";
		//os << "#lnpri\tlnL\t\tncomp\talpha\tlength\tstatent\trrent\n";
		os << "#lnpri\tlnL\t\tncomp\talpha\tcenterent\tcenter\tconcentration\tlength\tstatent\tstat\trrent\trr\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << aaprofilemutselmix->GetComponentNumber() << '\t';
		os << aaprofilemutselmix->GetAlpha() << '\t';
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
		//cout << "Done Trace\n";
		//cout.flush();
	}

	void ToStream(ostream& os)	{
		//cout << "Entering ToStream\n";
		//cout.flush();
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
		//cout << "Done ToStream\n";
		//cout.flush();
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

	int ObservedNonsynSubCount()	{
		phyloprocess->Sample();
		return phyloprocess->NonsynSubCount();
	}

	int PostPredNonsynSubCount()	{
		//phyloprocess->UnclampData();

		phyloprocess->SetData(dataunknown);
		phyloprocess->Sample();
		//phyloprocess->PostPredSample();
		int count = phyloprocess->NonsynSubCount();
		phyloprocess->SetData(codondata);

		//phyloprocess->ClampData();
		return count;
	}

};


