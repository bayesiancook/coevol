
#ifndef GENERALSEGMENTEDCONJUGATEPATH_H
#define	GENERALSEGMENTEDCONJUGATEPATH_H


#include "RandomSegmentedBranchSitePath.h"
#include "Conjugate.h"
#include "RandomSubMatrix.h"
#include "PhyloProcess.h"
#include "BrownianMove.h"

#include <utility>
#include <map>

//#include "CodonSubMatrix.h"



class SegmentPair {
public:
	int first;
	int second;
	int segment;
	SegmentPair(int infirst, int insecond, int insegment) {first=infirst; second=insecond;segment=insegment;}
	bool operator<(const SegmentPair &o) const {
		return (segment<o.segment)
			  || (segment==o.segment && (first <o.first
									|| (first == o.first && second < o.second)));
	} 
};

class SegmentedPathConjugate : public DSemiConjugatePrior<void>	{

	public:

	SegmentedPathConjugate(Var<PosReal>* inlength, EvolutionSegmentedBranch* inevol, bool inisroot)	{
		length = inlength;
		evol = inevol;
				matrix = evol->GetP();
		isroot = inisroot;
		Register(evol);
				Register(evol->GetGlobalTransMatrix());

				Register(length);
		CreateSuffStat();
	}

	~SegmentedPathConjugate() {
		DeleteSuffStat();
	}

	void specialUpdate()	{
	}

	// accessors?

	// PhyloProcess* GetPhyloProcess() {return myprocess;}
	EvolutionSegmentedMatrix** GetRandomMatrix() {return matrix;}

	//virtual AbstractTransitionMatrix* 	GetMatrix()		{return matrix;}
	int GetNstate() {return evol->GetNstate();}

	virtual const double* 		GetStationary()		{return evol->GetStationary();}

	virtual void CreateSuffStat()	{
	}

	virtual void DeleteSuffStat()	{
	}

	virtual void ResetSufficientStatistic()	{
		rootcount.clear();
		paircount.clear();
	}

	void SaveSufficientStatistic()	{
		cerr << "save sufficient\n";
		exit(1);
	}

	void RestoreSufficientStatistic()	{
		cerr << "restore sufficient\n";
		exit(1);
	}

	void IncrementRootCount(int state)	{
		rootcount[state]++;
	}

	void AddRootCount(int state, int count)	{
		rootcount[state] += count;
	}

	void IncrementPairCount(int state1, int state2, int segment)	{
		paircount[SegmentPair(state1,state2,segment)]++;
	}

	void AddPairCount(int state1, int state2, int segment, int count)	{
		paircount[SegmentPair(state1,state2,segment)] += count;
	}

	int GetPairCount(int state1, int state2, int segment)	{
		return  paircount[SegmentPair(state1,state2,segment)];
	}

	void PrintSuffStat(ostream& os)	{

		// if is root
		if (isroot)	{
			os << rootcount.size() << '\t';
			for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
				os << i->first << '\t' << i->second << '\t';
			}
			os << '\t';
		}
		else	{
			os << paircount.size() << '\t';
			for (map<SegmentPair, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				os << i->first.first << '\t' << i->first.second << '\t' << i->first.segment << '\t' << i->second << '\t';
			}
			os << '\t';
		}
	}

	void AddSuffStatFromStream(istream& is)	{

		// if is root
		if (isroot)	{
			int n;
			is >> n;
			for (int i=0; i<n; i++)	{
				int f, s;
				is >> f >> s;
				AddRootCount(f,s);
			}
		}
		else	{
			int n;
			is >> n;
			for (int i=0; i<n; i++)	{
				int f,s,seg,t;
				is >> f >> s >> seg >> t;
				AddPairCount(f,s,seg,t);
			}
		}
	}

	virtual double SuffStatLogProb()	{

		double total = 0;

		// root part
		if (isroot)	{
			const double* rootstat = GetStationary();
			for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
				total += i->second * log(rootstat[i->first]);
			}
		}

		// non root part
		if (! isroot)	{
			int totnsub = 0;
			for (map<SegmentPair, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				total += i->second * log( (*matrix[i->first.segment]) (i->first.first, i->first.second));
				totnsub += i->second;
			}
						if(std::isnan(total) || std::isinf(total)) {
							cerr << "suff stat log prob : " << total << endl;
							exit(0);
						}

		}

		return total;
  
	}

	double logProb()	{

		if (isActive())	{
					return SuffStatLogProb();
			return SuffStatLogProb();
		}
		return 0;
	}

	Var<PosReal>* GetLength() {return length;}
		EvolutionSegmentedBranch* GetEvol() {return evol;}


	protected:

	map<int,int> rootcount;
	map<SegmentPair, int> paircount;

		EvolutionSegmentedBranch* evol;
	Var<PosReal>* length;
	EvolutionSegmentedMatrix** matrix;
	bool isroot;
};


class SegmentedConjugateRandomBranchSitePath : public virtual ConjugateSampling<void>, public virtual RandomSegmentedBranchSitePath	{

	public:

	SegmentedConjugateRandomBranchSitePath(PhyloProcess* inprocess, SegmentedPathConjugate* inpathconj) : RandomSegmentedBranchSitePath(inprocess)	{
			evol = inpathconj->GetEvol();
		 
			interstate = new int[GetNsegment()+1];
			bkinterstate = new int[GetNsegment()+1];
			transitionmatrix = evol->GetGlobalTransMatrix();
			pathconj = inpathconj;

			if (!pathconj)	{
					cerr << "error in ConjugateRandomBranchSitePath\n";
					exit(1);
			}

			length = pathconj->GetLength();
			active_flag = true;

			conjugate_up.insert(pathconj);
			Register(pathconj);
	}

	protected:

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (parent != pathconj)	{
			cerr << "error in ConjugateRandomBranchSitePath::AddSufficientStatistic\n";
			exit(1);
		}

		if (isRoot())	{
			pathconj->IncrementRootCount(stateup);
		}
		else	{
					for(int i=0; i<GetNsegment(); i++)
			pathconj->IncrementPairCount(interstate[i],interstate[i+1], i);
		}
	}

	private:

	SegmentedPathConjugate* pathconj;
};


class SegmentedPathConjugateTree : public BranchValPtrTree<SegmentedPathConjugate>, public AbstractPathConj 	{

	public:

	SegmentedPathConjugateTree(LengthTree* inlengthtree, SequenceAlignment* indata)	{
		lengthtree = inlengthtree;
		data = indata;
	}

	~SegmentedPathConjugateTree()	{
	}

	LengthTree* GetLengthTree() {return lengthtree;}
	SequenceAlignment* GetData() {return data;}

	Tree* GetTree() {return lengthtree->GetTree();}

	void ActivateSufficientStatistic()	{
		ActivateSufficientStatistic(GetRoot());
	}

	void InactivateSufficientStatistic()	{
		InactivateSufficientStatistic(GetRoot());
	}

	void PrintSuffStat(ostream& os)	{
		RecursivePrintSuffStat(os,GetRoot());
		os << '\n';
	}

	double GetLogProb()	{
				return RecursiveGetLogProb(GetRoot());
	}

	void ReadFromFile(string infile)	{
		ifstream is(infile.c_str());
		ActivateSufficientStatistic();
		int N;
		is >> N;
		for (int i=0; i<N; i++)	{
			RecursiveAddSuffStatFromStream(is,GetRoot());
		}
		/*
		ofstream os("check");
		PrintSuffStat(os);
		*/
	}

	protected:

	double RecursiveGetLogProb(const Link* from)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogProb(link->Out());
		}
		total += GetBranchVal(from->GetBranch())->logProb();
		return total;
	}

	void RecursiveAddSuffStatFromStream(istream& is, const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSuffStatFromStream(is,link->Out());
		}
		GetBranchVal(from->GetBranch())->AddSuffStatFromStream(is);
	}

	void RecursivePrintSuffStat(ostream& os, const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivePrintSuffStat(os,link->Out());
		}
		GetBranchVal(from->GetBranch())->PrintSuffStat(os);
	}

	SegmentedPathConjugate* GetSegmentedPathConjugate(const Branch* branch)	{
		return dynamic_cast<SegmentedPathConjugate*>(GetBranchVal(branch));
	}

	void ActivateSufficientStatistic(const Link* from)	{
		GetSegmentedPathConjugate(from->GetBranch())->ActivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ActivateSufficientStatistic(link->Out());
		}
	}

	void InactivateSufficientStatistic(const Link* from)	{
		GetSegmentedPathConjugate(from->GetBranch())->InactivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			InactivateSufficientStatistic(link->Out());
		}
	}

	LengthTree* lengthtree;
	SequenceAlignment* data;

};

class BranchMatrixSegmentedPathConjugateTree : public SegmentedPathConjugateTree	{

	public:

	BranchMatrixSegmentedPathConjugateTree(LengthTree* inlengthtree, EvolutionSegmentedProcess* inevoltree,  SequenceAlignment* indata) : SegmentedPathConjugateTree(inlengthtree, indata) {
		evoltree = inevoltree;
		SetWithRoot(true);
		if (! evoltree->WithRoot())	{
			cerr << "error in SegmentedPathConjugateTree: evoltree does not have root value\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~BranchMatrixSegmentedPathConjugateTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	SegmentedPathConjugate* CreateBranchVal(const Link* link)	{
		return new SegmentedPathConjugate(lengthtree->GetBranchVal(link->GetBranch()),evoltree->GetBranchVal(link->GetBranch()), link->isRoot());
	}

	EvolutionSegmentedProcess* evoltree;

};

class SegmentedPathConjugatePhyloProcess : public PhyloProcess	{


	protected:

	public:

	SegmentedPathConjugatePhyloProcess(SegmentedPathConjugateTree* inpathconjtree) : PhyloProcess(inpathconjtree->GetLengthTree(), inpathconjtree->GetData(), false)	{
		pathconjtree = inpathconjtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		if (! pathconjtree->GetBranchVal(link->GetBranch()))	{
			cerr << "error in create branch val\n";
			exit(1);
		}
		return  new SegmentedConjugateRandomBranchSitePath(this,pathconjtree->GetBranchVal(link->GetBranch()));
	}

	protected:
	SegmentedPathConjugateTree* pathconjtree;

};







#endif	/* GENERALSEGMENTEDCONJUGATEPATH_H */

