
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
#include "MGCodonProfileMutSelCodonSubMatrix.h"
#include "MatrixMixture.h"
#include "CodonProfileMutSelMixtureRandomMatrix.h"
#include "ProfileInfiniteMixturePhyloProcess.h"

class CodonProfileMutSelMatrixInfiniteMixture : public MatrixInfiniteMixture<Profile>	{

	public:

	CodonProfileMutSelMatrixInfiniteMixture(int insize, int incomponentnumber, Var<Profile>* incenter, Var<PosReal>* inconcentration, CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, Var<PosReal>* inNeff) :
			MatrixInfiniteMixture<Profile>(insize, incomponentnumber)	{

		center = incenter;
		concentration = inconcentration;
		statespace = instatespace;
		nucmatrix = innucmatrix;
		Neff = inNeff;
		Create();
	}

	double MoveValues(double tuning, int n)	{
		double total = 0;
		for (int k=0; k<GetComponentNumber(); k++)	{
			total += GetDirichlet(k)->Move(tuning,n);
		}
		return total / GetComponentNumber();
	};


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
		CodonProfileMutSelMixtureRandomMatrix* temp = new CodonProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
		return temp;
		//return new AAProfileMutSelMixtureRandomMatrix(center,concentration,statespace,nucmatrix, Neff);
	}

	private:
	Var<Profile>* center;
	Var<PosReal>* concentration;
	Var<PosReal>* Neff;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};

class CodonProfileMutSelMatInfMixValMove : public MCUpdate	{

	public:

	CodonProfileMutSelMatInfMixValMove(CodonProfileMutSelMatrixInfiniteMixture* inmix, double intuning, int inN, int innrep) : mix(inmix), tuning(intuning), N(inN), nrep(innrep) {}

	double Move(double tuning_modulator = 1)	{
		double total = 0;
		for (int rep=0; rep<nrep; rep++)	{
			total += mix->MoveValues(tuning * tuning_modulator,N);
		}
		return total / nrep;
	}

	protected:

	CodonProfileMutSelMatrixInfiniteMixture* mix;
	double tuning;
	int N;
	int nrep;
};

/*
class ProfileInfiniteMixturePhyloProcess : public MatrixInfiniteMixturePhyloProcess<Profile>	{

	public:

	ProfileInfiniteMixturePhyloProcess(LengthTree* intree, MatrixInfiniteMixture<Profile>* inmatmix,  CodonSequenceAlignment* indata) : MatrixInfiniteMixturePhyloProcess<Profile>(intree,inmatmix,indata)	{data=indata;}

	int BranchSiteNonsynSubCount(Link* from, int site)      {
		int count = 0;
		for (Link* blink=from->Next(); blink!=from; blink=blink->Next())        {
			Plink* plink = GetPath(blink->GetBranch(), site)->Init();
			while (plink != GetPath(blink->GetBranch(), site)->Last())      {
				if (! data->GetCodonStateSpace()->Synonymous(plink->GetState(), plink->Next()->GetState()))     {
					count++;
				}
				plink=plink->Next();
			}
			count += BranchSiteNonsynSubCount(blink->Out(), site);
		}
		return count;
        }

	double NonsynSubMean()    {
		int count = 0;
		for (int site=0; site<GetNsite(); site++)       {
			count += BranchSiteNonsynSubCount(GetRoot(), site);
		}
		return ( (double)(count)/GetNsite() );
	}

	double NonsynSubVariance()	{
		int count = 0;
		int sum = 0;
		double var = 0;
		for (int site=0; site<GetNsite(); site++)	{
			count = BranchSiteNonsynSubCount(GetRoot(), site);
			var += (count * count);
			sum += count;
		}
		double mean = (double)(sum)/GetNsite();
		var /= GetNsite();
		var -= (mean * mean);
		return var;
	}

	private:
	CodonSequenceAlignment* data;

};
*/


class CodonProfileMutSelMatInfMixModel : public ProbModel {

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
	Dirichlet* CodonProfileCenter;
        Exponential* CodonProfileConcentration;
        Const<PosReal>* CodonProfileConcentrationPrior;

	public :

	CodonProfileMutSelMatrixInfiniteMixture* codonprofilemutselmix;
	//MatrixInfiniteMixturePhyloProcess<Profile>* phyloprocess;
	ProfileInfiniteMixturePhyloProcess* phyloprocess;

	// constructor
	// this is where the entire graph structure of the model is created
	CodonProfileMutSelMatInfMixModel(string datafile, string treefile, int inP, bool sample=true, GeneticCodeType type=Universal)	{
		// fetch data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		//datacopy = 0;
		datacopy = new CodonSequenceAlignment(nucdata, true, type);
		datacopy->Unclamp();
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)


		//ofstream osaadata("tempAA");
		//codondata->ToStream(osaadata);

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
		lognormaltree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
		lognormaltree->Reset();
		*/

		// a tree with all branch lengths iid from an exponential distribution of rate lambda
		// this is a gamma distribution of shape mu=1 and scale lambda
		// lambda is itself endowed with an exponential prior of mean 10
		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,mu,lambda);
		//*
		gamtree->SetBranchLengths();

		cout << "Tree length: " << GetLength() << "\n";
		cout.flush();

		/*
		gamtree->Clamp();
		*/

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		relrate->setuniform();
		stationary->setuniform();
		neff = new Const<PosReal>(1);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary);

		CodonStateSpace* codonstatespace = (CodonStateSpace*) codondata->GetStateSpace();

		cerr << "infinite mixture\n";
		P = inP;
		cerr << P << '\n';
		CodonProfileCenter = new Dirichlet(Nstate);
		CodonProfileCenter->setuniform();
                CodonProfileConcentrationPrior = new Const<PosReal>((double)(Nstate));
                CodonProfileConcentration = new Exponential(CodonProfileConcentrationPrior, Exponential::MEAN);
		CodonProfileConcentration->setval(Nstate);
		//CodonProfileConcentration->ClampAt(Nstate);

		codonprofilemutselmix = new CodonProfileMutSelMatrixInfiniteMixture(Nsite,P,CodonProfileCenter,CodonProfileConcentration,codonstatespace,nucmatrix, neff);

		// Update before phyloprocess
		Update();

		cerr << "create phylo process\n";
		//phyloprocess = new MatrixInfiniteMixturePhyloProcess<Profile>(lognormaltree,aaprofilemutselmix,codondata);
		//phyloprocess = new ProfileInfiniteMixturePhyloProcess(lognormaltree,aaprofilemutselmix,codondata);
		phyloprocess = new ProfileInfiniteMixturePhyloProcess(gamtree,codonprofilemutselmix,codondata);
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
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(CodonProfileCenter);
		RootRegister(CodonProfileConcentrationPrior);
		RootRegister(neff);
		//RootRegister(aaprofilemutselmix->GetWeightVector());
		Register();

		cout << "after Register\n";
		cout.flush();

		MakeScheduler();

		cout << "after MakeScheduler()\n";
		cout.flush();

		if (sample)	{
			Update();
		}

		cout << "after Update()\n";
		cout.flush();

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
	~CodonProfileMutSelMatInfMixModel() {}

	Tree* GetTree() {return tree;}
	/*
	LogNormalTreeProcess* GetLogNormalTree() {return lognormaltree;}
	LengthTree* GetChronogram() {return chronogram;}
	*/
	int GetNtaxa() {return  Ntaxa;}
	int GetNstate() {return Nstate;}
	CodonSequenceAlignment* GetCodonData() {return codondata;}


	//double Move(double tuning = 1)	{
	//	scheduler.Cycle(1,1,true,true);
	//	return 1;
	//}


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

		//total += codonprofilemutselmix->GetLogProb();
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
		/*
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
		*/

		scheduler.Register(new MatInfMixAllocMove<Profile>(codonprofilemutselmix,5,1),1,"codonprofilemutsel mix weight alloc"); // comment out to fix allocations
		scheduler.Register(new CodonProfileMutSelMatInfMixValMove(codonprofilemutselmix,1,1,10),1,"codonprofilemutsel mix subset");
		scheduler.Register(new CodonProfileMutSelMatInfMixValMove(codonprofilemutselmix,1,2,10),1,"codonprofilemutsel mix subset");
		scheduler.Register(new CodonProfileMutSelMatInfMixValMove(codonprofilemutselmix,0.3,3,10),1,"codonprofilemutsel mix subset");
		scheduler.Register(new MatInfMixAlphaMove<Profile>(codonprofilemutselmix,0.1,1),10,"codonprofilemutsel alpha");
		scheduler.Register(new MatInfMixAlphaMove<Profile>(codonprofilemutselmix,0.2,1),10,"codonprofilemutsel alpha");
		scheduler.Register(new MatInfMixAlphaMove<Profile>(codonprofilemutselmix,0.3,1),10,"codonprofilemutsel alpha");
		scheduler.Register(new MatInfMixAlphaMove<Profile>(codonprofilemutselmix,0.4,1),10,"codonprofilemutsel alpha");
		scheduler.Register(new MatInfMixAlphaMove<Profile>(codonprofilemutselmix,0.5,1),10,"codonprofilemutsel alpha");
		scheduler.Register(new MatInfMixAllocMove<Profile>(codonprofilemutselmix,5,1),1,"codonprofilemutsel mix weight alloc"); // comment out
		scheduler.Register(new CodonProfileMutSelMatInfMixValMove(codonprofilemutselmix,0.3,4,10),1,"codonprofilemutsel mix subset");
		scheduler.Register(new CodonProfileMutSelMatInfMixValMove(codonprofilemutselmix,0.1,5,1),10,"codonprofilemutsel mix subset");


		scheduler.Register(new ProfileMove(CodonProfileCenter,0.05,5),10,"codonprofile center profile move");
		scheduler.Register(new ProfileMove(CodonProfileCenter,0.06,5),10,"codonprofile center profile move");
		scheduler.Register(new ProfileMove(CodonProfileCenter,0.07,5),10,"codonprofile center profile move");
		scheduler.Register(new ProfileMove(CodonProfileCenter,0.08,5),10,"codonprofile center profile move");
		scheduler.Register(new ProfileMove(CodonProfileCenter,0.09,5),10,"codonprofile center profile move");
		scheduler.Register(new ProfileMove(CodonProfileCenter,0.1,5),10,"codonprofile center profile move");

		scheduler.Register(new SimpleMove(CodonProfileConcentration,0.1),10,"codonprofile concentration");
		scheduler.Register(new SimpleMove(CodonProfileConcentration,0.2),10,"codonprofile concentration");
		scheduler.Register(new SimpleMove(CodonProfileConcentration,0.3),10,"codonprofile concentration");
		scheduler.Register(new SimpleMove(CodonProfileConcentration,0.4),10,"codonprofile concentration");
		scheduler.Register(new SimpleMove(CodonProfileConcentration,0.5),10,"codonprofile concentration");
		scheduler.Register(new SimpleMove(CodonProfileConcentration,1.0),10,"codonprofile concentration");

		/*
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.2),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.3),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.4),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.5),100,"sigma");
		scheduler.Register(new SimpleMove(lognormaltree,1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.1),10,"lognormal");
		scheduler.Register(new SimpleMove(lognormaltree,0.01),10,"lognormal");
		*/


		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");
		scheduler.Register(new SimpleMove(gamtree,1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.1),10,"branch lengths");
		scheduler.Register(new SimpleMove(gamtree,0.01),10,"branch lengths");

		//scheduler.Register(new ProfileMove(relrate,1,1),10,"relrates");
		//scheduler.Register(new ProfileMove(relrate,0.3,2),10,"relrates");
		//scheduler.Register(new ProfileMove(relrate,0.1,3),10,"relrates");
		//scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		//scheduler.Register(new SimpleMove(relrate,0.03),10,"relrates");
		//scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates2");
		scheduler.Register(new ProfileMove(relrate,0.1,2),10,"relrates4");
		scheduler.Register(new ProfileMove(relrate,0.1,3),10,"relrates6");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relratesS");
		scheduler.Register(new SimpleMove(relrate,0.02),10,"relratesS");
		scheduler.Register(new SimpleMove(relrate,0.03),10,"relratesS");


		scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.02,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
		scheduler.Register(new ProfileMove(stationary,0.04,2),10,"stat4");
		scheduler.Register(new SimpleMove(stationary,0.01),10,"statS");
		scheduler.Register(new SimpleMove(stationary,0.02),10,"statS");
		scheduler.Register(new SimpleMove(stationary,0.03),10,"statS");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");

	}

	// Draw a sample from the prior

	void drawSample()	{
		cerr << "sample\n";
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
		CodonProfileConcentration->Sample();
		CodonProfileCenter->Sample();
		codonprofilemutselmix->Sample();
		cerr << "iid\n";
		// iidarray->Sample();
		cout << "in drawSample before exit\n";
		exit(1);
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
		//os << "#logprior\tlnL\tlength\tstatent\trrent\tcompnum\n";
		os << "#lnpri\tlnL\t\tncomp\talpha\tcenterent\tcenter\tconcentration\tlength\tstatent\tstat\trrent\trr\n";
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
		os << codonprofilemutselmix->GetComponentNumber() << '\t';
		os << codonprofilemutselmix->GetAlpha() << '\t';
		os << CodonProfileCenter->val().GetEntropy() << '\t';
		os << *CodonProfileCenter << '\t';
		os << *CodonProfileConcentration << '\t';
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
		os << *CodonProfileCenter << '\n';
		os << *CodonProfileConcentration << '\n';
		os << *codonprofilemutselmix << '\n';
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
		is >> *CodonProfileCenter;
		is >> *CodonProfileConcentration;
		is >> *codonprofilemutselmix;

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
		is >> *CodonProfileCenter;
		is >> *CodonProfileConcentration;
		is >> *codonprofilemutselmix;

		//cout << "FromStream\n";
		//cout.flush();

		for (int i=0; i<codonprofilemutselmix->GetSize(); i++)	{
			codonprofilemutselmix->SetAllocation(i, codonprofilemutselmix->GetAllocation(i));
		}

		Update();
		phyloprocess->Sample();
	}

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

	void SamplePosteriorMapping()	{
		phyloprocess->Sample();
	}

	void SamplePosteriorPredictiveMapping()	{
		phyloprocess->SetData(datacopy);
		phyloprocess->Sample();
		phyloprocess->SetData(codondata);
	}

	double NonsynSubMean()	{
		//phyloprocess->Sample();
		return phyloprocess->NonsynSubMean();
	}

	double NonsynSubVariance()	{
		return phyloprocess->NonsynSubVariance();
	}

	//double PostPredNonsynSubMean()	{

		//phyloprocess->SetData(datacopy);
		//phyloprocess->Sample();
		//double meanns = phyloprocess->NonsynSubMean();
		//phyloprocess->SetData(codondata);
		//return meanns;
	//}

	double MeanObservedAADiversity()	{

		CodonSequenceAlignment* data = codondata;
		int* AAPresence = new int[Naa];
		for (int aa=0; aa<Naa; aa++)	{
			AAPresence[aa] = 0;
		}

		int AADiversity;
		int AADiversitySum = 0;
		for (int site=0; site<Nsite; site++)	{
			for (int taxa=0; taxa<Ntaxa; taxa++)	{
				int state = data->GetState(taxa,site);
				if (state != unknown)	{
					AAPresence[data->GetCodonStateSpace()->Translation(state)] = 1;
				}
			}

			AADiversity = 0;
			for (int aa=0; aa<Naa; aa++)	{
				AADiversity += AAPresence[aa];
				AAPresence[aa] = 0;
			}
			AADiversitySum += AADiversity;
		}

		delete[] AAPresence;
		return ((double)(AADiversitySum)/Nsite);
	}

	double MeanPredictiveAADiversity()	{

		CodonSequenceAlignment* data = datacopy;
		int* AAPresence = new int[Naa];
		for (int aa=0; aa<Naa; aa++)	{
			AAPresence[aa] = 0;
		}

		int AADiversity;
		int AADiversitySum = 0;
		for (int site=0; site<Nsite; site++)	{
			for (int taxa=0; taxa<Ntaxa; taxa++)	{
				int state = data->GetState(taxa,site);
				if (state != unknown)	{
					AAPresence[data->GetCodonStateSpace()->Translation(state)] = 1;
				}
			}

			AADiversity = 0;
			for (int aa=0; aa<Naa; aa++)	{
				AADiversity += AAPresence[aa];
				AAPresence[aa] = 0;
			}
			AADiversitySum += AADiversity;
		}
		delete[] AAPresence;
		return ((double)(AADiversitySum)/Nsite);
	}


	void EmpiricalSiteSpecificFrequencies(ostream& os)	{


		os << Nsite << '\n';

		CodonSequenceAlignment* data = codondata;
		double* AAProfile = new double[Naa];
		for (int aa=0; aa<Naa; aa++)	{
			AAProfile[aa] = 0;
		}
		for (int site=0; site<Nsite; site++)	{

			for (int aa=0; aa<Naa; aa++)	{
				AAProfile[aa] = 0;
			}

			double sum = 0;

			for (int taxa=0; taxa<Ntaxa; taxa++)	{
				int state = data->GetState(taxa,site);

				if (state != unknown)	{
					AAProfile[data->GetCodonStateSpace()->Translation(state)] += 1.0;
					sum += 1.0;
				}
			}

			for (int aa=0; aa<Naa; aa++)	{
				os << AAProfile[aa]/sum << '\t';
			}
			os << '\n';
		}

	}
};

