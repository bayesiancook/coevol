
#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include "MCMC.h"
#include "BranchProcess.h"
#include "SiteMapping.h"
#include "RandomBranchSitePath.h"
#include "SequenceAlignment.h"
#include "BranchSiteSubstitutionProcess.h"
#include <map>

#include "Chrono.h"

// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// of the relevant type, and using the relevant parameters (on a model specific basis)
// for each branch/site pair
// this creation of branch/site paths is done through a pure virtual method: CreateRandomBranchSitePath
//
// PhyloProcess and RandomBranchSitePath are two tightly interconnected abstract classes
// to build new phylogenetic model:
// - make a new class deriving from RandomBranchSitePath (say MyRandomBranchSitePath)) with the relevant constructor (where the connections with parameters will be made)
// - make a new class deriving from PhyloProcess (say MyPhyloProcess), implementing CreateRandomBranchSitePath
//   this implementation of CreateRandomBranchSitePath should essentially call the constructor of MyRandomBranchSitePath, with the parameters that are relevant for that particular branch/site pair:
//
//   class MyRandomBranchSitePath : public RandomBranchSitePath	{
//   
//   		MyRandomBranchSitePath( ... model specific list of parameters...)	{
//
//   			// your code here
//   		}
//   };
//
//
//   class MyPhyloProcess : public PhyloProcess	{
//
//   		RandomBranchSitePath(Link* link, int site)	{
//   			return new MyRandomBranchSitePath( ... here model-specific list of arguments ... );
//   		}
//   };
//
//
//
class PhyloProcess : public MCMC	{

	friend class ProbModel;

	public:

	// pure virtual method: this is where all the model-specific things have to be implemented
	virtual RandomBranchSitePath* CreateRandomBranchSitePath(const Link* link, int site) = 0;


				PhyloProcess(LengthTree* intree, SequenceAlignment* indata, bool inbranchmap = true);
	virtual 		~PhyloProcess();

	double 			GetLogProb(); 		// likelihood Felsenstein 1981
	double 			GetFastLogProb(); 		// likelihood Felsenstein 1981
	double 			GetPathLogProb(); 	// probability of the entire mapping
	double 			GetPathLogProb(const Link* from, int site); 	// probability of the entire mapping
	double 			SiteLogLikelihood(int site);
	double 			FastSiteLogLikelihood(int site);

	void			drawSample(); 		// ResampleSub Nielsen: accept reject algorithm for stochastic mappings (clamped)
	double			Move(double tuning); 	// ResampleSub (clamped)

	StateSpace*		GetStateSpace() {return data->GetStateSpace();}
	void			SetStateSpace();
	void			RecursiveSetStateSpace(const Link* from);
	SequenceAlignment*	GetData() {return data;}

	double PrintRootLikelihood(ostream& os);

	int 			GetNsite();
	int 			GetNtaxa();
	int 			GetNstate(int site);
	int			GetMaxNstate();

	int			GetMaxTrial() { return maxtrial;}
	void			SetMaxTrial(int i) { maxtrial = i;}

	void 			ResampleSub();	// clamped Nielsen
	void 			ResampleSub(int site);

	void			PostPredSample(bool rootprior = false); // unclamped Nielsen
	void			PostPredSample(int site, bool rootprior = false);
				// rootprior == true : root state drawn from stationary probability of the process
				// rootprior == false: root state drawn from posterior distribution

	// computes the frequencies of states in each taxon
	// the global frequencies
	// and returns the chi-square score
	double			CompositionalHeterogeneityIndex(ostream& os);
	void			GetLeafFreqs(const Link* from, double** taxfreq);

	void			GetLeafData(SequenceAlignment* data);
	void			RecursiveGetLeafData(const Link* from, SequenceAlignment* data);

	void 			ClampData() {clampdata = true;}
	void			UnclampData() {clampdata = false;}

	// various accessors

	BranchSiteSubstitutionProcess* 		GetBranchSiteSubstitutionProcess(const Branch* branch, int site)	{
					if (isMissing(branch,site))	{
						cerr << "as bssub\n";
						exit(1);
					}
					return (BranchSiteSubstitutionProcess*) GetPath(branch,site);
				}

	public:
	// bool			NonMissingPath(const Branch* branch, int site);
	int&			GetState(const Node* node, int site);
	int			GetData(int taxon, int site);
	bool			isDataCompatible(int taxon, int site, int state)	{
					return GetStateSpace()->isCompatible(GetData(taxon,site),state);
				}

	Tree*			GetTree();
	Link*			GetRoot();
	TaxonSet*		GetTaxonSet();

	void			SetData(SequenceAlignment* data);
	virtual void		Unfold();

	void			SetName(string inname);
	void			RecursiveSetName(const Link* from, string inname);

	int 			GetTotMissing(const Node* node)	{
					return totmissingmap[node];
				}

	bool			isMissing(const Node* node, int site)	{
					return missingmap[node][site];
				}

	bool			isMissing(const Link* link, int site)	{
					if ((! missingmap[link->GetNode()][site]) && (! missingmap[link->Out()->GetNode()][site]) && (pathmap[link->GetBranch()]) && (! pathmap[link->GetBranch()][site]))	{
						cerr << "error in is missing\n";
						exit(1);
					}
					return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
				}

	bool			isMissing(const Branch* branch, int site)	{
					return (! pathmap[branch][site]);
				}

	void 			CreateMissingMap();
	void			RecursiveCreateMissingMap(const Link* from);
	bool			FillMissingMap(const Link* from, int site);
	void			ComputeTotalMissingPerSite(const Link* from);

	double*			GetCondLikelihood(const Link* from)	{
					return condlmap[from];
	}

	/*
	double			GetTotalTime()	{
		return chrono.GetTime();
	}
	*/

	double			GetPruningTime()	{
		return pruningchrono.GetTime();
	}

	/*
	double			GetPruningAncestrsalTime()	{
		return prunancchrono.GetTime();

	}
	*/

	double			GetResampleTime()	{
		return resamplechrono.GetTime();
	}

	protected:


	RandomBranchSitePath**  	CreateRandomBranchSitePath(const Link* link);
	void			DeletePath(const Link* link);

	void			Cleanup();

	void			RecursiveCreate(const Link* link);
	void			RecursiveDelete(const Link* link);

	void 			RecursiveCreateTBL(const Link* link, int innstate);
	void 			RecursiveDeleteTBL(const Link* link);

	void 			Pruning(const Link* link, int site);
	void 			ResampleSub(const Link* link, int site);
	void 			ResampleState();
	void			ResampleState(int site);
	void			PruningAncestral(const Link* link, int site);
	double			RecordPruningAncestralLogProb(const Link* link, int site);
	void			PriorSample(const Link* link, int site, bool rootprior);
	void			PriorSample();
	void			RootPosteriorDraw(int site);

	public:
	// mappings
	class PhyloProcessSiteMapping : public SiteMapping	{

		public:
					PhyloProcessSiteMapping(PhyloProcess* inprocess, int i);
					Tree* GetTree();
					BranchSitePath* GetPath(const Branch* branch);

		private:
		PhyloProcess* myprocess;
		int site;
	};

	// to access the current site mapping at site <site>
	PhyloProcessSiteMapping*	GetSiteMapping(int site);

	bool SampleBranchMapping()	{
		return branchmap;
	}


	// data fields

	LengthTree* tree;
	SequenceAlignment* data;

	int* sitearray;
	double* sitelnL;
	map<const Link* , double*> condlmap;

	int MaxNstate;

	bool clampdata;

	double MHMove(int nrep, double s01, double s10);
	void SetProposalMatrices();

	RandomBranchSitePath*	GetPath(const Branch* branch, int site);

	private:
	map<const Branch*,RandomBranchSitePath**> pathmap;
	map<const Node*,int*> statemap;
	map<const Node*, int> bkstatemap;
	map<const Node*,bool*> missingmap;
	map<const Node*,int>  totmissingmap;

	void RecursiveSetProposalMatrices(const Link* from);

	// to be overriden in Metropolis Hastings PhyloProcess classes
	// EmpiricalSubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
	virtual SubMatrix* GetProposalMatrix(const Branch* branch, int site)	{
		return 0;
	}

	void BackupNodeStates(const Link* from, int site);
	void RestoreNodeStates(const Link* from, int site);

	double ProposeResampleSub(const Link* from, int site);
	void SwapMatrices(const Link* from, int site);
	void RecursiveRegister(const Link* from, int site, Mnode* mnode);
	void ResetFlagMap(const Link* from, bool in);
	const Link* ChooseNodeSet(const Link* from, double s01, double s10, int sw);
	map<const Node*,bool> flagmap;

	int maxtrial;
	PhyloProcessSiteMapping** sitemapping;
	static const int unknown = -1;

	static const int DEFAULTMAXTRIAL = 100;

	// Chrono chrono;
	Chrono pruningchrono;
	// Chrono prunancchrono;
	Chrono resamplechrono;

	bool branchmap;
};


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* PhyloProcessSiteMapping
//-------------------------------------------------------------------------

inline PhyloProcess::PhyloProcessSiteMapping::PhyloProcessSiteMapping(PhyloProcess* inprocess, int i) : myprocess(inprocess), site(i) {}
inline Tree* PhyloProcess::PhyloProcessSiteMapping::GetTree() {return myprocess->GetTree();}
inline BranchSitePath* PhyloProcess::PhyloProcessSiteMapping::GetPath(const Branch* branch) {
		if (myprocess->isMissing(branch,site))	{
			cerr << "PhyloProcessSiteMapping::GetPath called on null path\n";
			exit(1);
		}
		return myprocess->GetPath(branch,site);
}

//-------------------------------------------------------------------------
//	* PhyloProcess
//-------------------------------------------------------------------------

inline int PhyloProcess::GetNsite()	{return data->GetNsite();}
inline int PhyloProcess::GetNtaxa() {return data->GetNtaxa();}
inline int PhyloProcess::GetNstate(int site) {
	if (isMissing(GetRoot(),site))	{
		return 0;
	}
	return GetBranchSiteSubstitutionProcess(0,site)->GetNstate();
}

inline RandomBranchSitePath* PhyloProcess::GetPath(const Branch* branch, int site) {
	if (! pathmap[branch][site])	{
		cerr << "error in phyloprocess::getpath: null path\n";
		exit(1);
	}
	return pathmap[branch][site];
}

/*
inline bool PhyloProcess::NonMissingPath(const Branch* branch, int site)	{
	return (bool) pathmap[branch][site];
}
*/

inline int& PhyloProcess::GetState(const Node* node, int site) {return statemap[node][site];}
inline int PhyloProcess::GetData(int taxon, int site)	{return data->GetState(taxon,site);}
inline Tree* PhyloProcess::GetTree() {return tree->GetTree();}
inline Link* PhyloProcess::GetRoot() {return GetTree()->GetRoot();}
inline TaxonSet* PhyloProcess::GetTaxonSet() {return data->GetTaxonSet();}

inline PhyloProcess::PhyloProcessSiteMapping* PhyloProcess::GetSiteMapping(int site) {return sitemapping[site];}


#endif // PHYLOPROCESS_H
