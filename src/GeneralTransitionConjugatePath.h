#ifndef CONJUGATETRANSPATH_H
#define CONJUGATETRANSPATH_H

#include "RandomBranchSitePath.h"
#include "Conjugate.h"
#include "RandomSubMatrix.h"
#include "PhyloProcess.h"

#include <utility>
#include <map>

#include "CodonSubMatrix.h"

class AbstractPathConj {
public:
	virtual void ActivateSufficientStatistic() = 0;
	virtual void InactivateSufficientStatistic() = 0;
	virtual ~AbstractPathConj(){}
};

class TransitionPathConjugate : public DSemiConjugatePrior<void>	{

	public:

	TransitionPathConjugate(Var<PosReal>* inlength, RandomTransitionMatrix* inmatrix, bool inisroot)	{
		length = inlength;
		matrix = inmatrix;
		isroot = inisroot;
		Register(matrix);
		CreateSuffStat();
	}

	~TransitionPathConjugate() {
		DeleteSuffStat();
	}

	void specialUpdate()	{
	}

	// accessors?

	// PhyloProcess* GetPhyloProcess() {return myprocess;}
	RandomTransitionMatrix* GetRandomMatrix() {return matrix;}

	virtual AbstractTransitionMatrix* 	GetMatrix()		{return matrix;}
	int GetNstate() {return matrix->GetNstate();}

	virtual const double* 		GetStationary()		{return matrix->GetStationary();}

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

	void IncrementPairCount(int state1, int state2)	{
		paircount[ pair<int,int>(state1,state2)]++;
	}

	void AddPairCount(int state1, int state2, int count)	{
		paircount[ pair<int,int>(state1,state2)] += count;
	}

	int GetPairCount(int state1, int state2)	{
		return  paircount[ pair<int,int>(state1,state2)];
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
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				os << i->first.first << '\t' << i->first.second << '\t' << i->second << '\t';
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
				int f,s,t;
				is >> f >> s >> t;
				AddPairCount(f,s,t);
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
			AbstractTransitionMatrix& mat = *GetMatrix();
			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				total += i->second * log(mat(i->first.first, i->first.second));
				totnsub += i->second;
			}
		}
		return total;
	}

	double logProb()	{
		if (isActive())	{
			return SuffStatLogProb();
		}
		return 0;
	}

	Var<PosReal>* GetLength() {return length;}

	protected:

	map<int,int> rootcount;
	map< pair<int,int>, int> paircount;

	Var<PosReal>* length;
	RandomTransitionMatrix* matrix;
	bool isroot;
};


class TransitionConjugateRandomBranchSitePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath	{

	public:

	TransitionConjugateRandomBranchSitePath(PhyloProcess* inprocess, TransitionPathConjugate* inpathconj) : RandomBranchSitePath(inprocess)	{

		pathconj = inpathconj;

		if (!pathconj)	{
			cerr << "error in ConjugateRandomBranchSitePath\n";
			exit(1);
		}

		transitionmatrix = pathconj->GetRandomMatrix();
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
			pathconj->IncrementPairCount(stateup,statedown);
		}
	}

	private:

	TransitionPathConjugate* pathconj;
};


class TransitionPathConjugateTree : public BranchValPtrTree<TransitionPathConjugate>, public AbstractPathConj {

	public:

	TransitionPathConjugateTree(LengthTree* inlengthtree, SequenceAlignment* indata)	{
		lengthtree = inlengthtree;
		data = indata;
	}

	~TransitionPathConjugateTree()	{
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

	TransitionPathConjugate* GetTransitionPathConjugate(const Branch* branch)	{
		return dynamic_cast<TransitionPathConjugate*>(GetBranchVal(branch));
	}

	void ActivateSufficientStatistic(const Link* from)	{
		GetTransitionPathConjugate(from->GetBranch())->ActivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ActivateSufficientStatistic(link->Out());
		}
	}

	void InactivateSufficientStatistic(const Link* from)	{
		GetTransitionPathConjugate(from->GetBranch())->InactivateSufficientStatistic();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			InactivateSufficientStatistic(link->Out());
		}
	}

	LengthTree* lengthtree;
	SequenceAlignment* data;

};

class BranchMatrixTransitionPathConjugateTree : public TransitionPathConjugateTree	{

	public:

	BranchMatrixTransitionPathConjugateTree(LengthTree* inlengthtree, BranchValPtrTree<RandomTransitionMatrix>* inmatrixtree,  SequenceAlignment* indata) : TransitionPathConjugateTree(inlengthtree, indata) {
		matrixtree = inmatrixtree;
		SetWithRoot(true);
		if (! matrixtree->WithRoot())	{
			cerr << "error in TransitionPathConjugateTree: matrixtree does not have root value\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~BranchMatrixTransitionPathConjugateTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	TransitionPathConjugate* CreateBranchVal(const Link* link)	{
		return new TransitionPathConjugate(lengthtree->GetBranchVal(link->GetBranch()),matrixtree->GetBranchVal(link->GetBranch()), link->isRoot());
	}

	BranchValPtrTree<RandomTransitionMatrix>* matrixtree;

};

class TransitionPathConjugatePhyloProcess : public PhyloProcess	{


	protected:

	public:

	TransitionPathConjugatePhyloProcess(TransitionPathConjugateTree* inpathconjtree) : PhyloProcess(inpathconjtree->GetLengthTree(), inpathconjtree->GetData(), false)	{
		pathconjtree = inpathconjtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		if (! pathconjtree->GetBranchVal(link->GetBranch()))	{
			cerr << "error in create branch val\n";
			exit(1);
		}
		return  new TransitionConjugateRandomBranchSitePath(this,pathconjtree->GetBranchVal(link->GetBranch()));
	}

	protected:
	TransitionPathConjugateTree* pathconjtree;

};

template <class R, class D> class DSemiConjugateMove : public MCUpdate	{

	public:

	DSemiConjugateMove(R* inrandom, D* indsemi, double intuning, int inn) : random(inrandom), dsemi(indsemi), tuning(intuning), n(inn) {}

	double Move(double tuning_modulator = 1)	{
		dsemi->ActivateSufficientStatistic();
		double total = 0;
		for (int i=0; i<n; i++)	{
			total += random->Move(tuning* tuning_modulator);
		}
		total /= n;
		dsemi->InactivateSufficientStatistic();
		return total;
	}

	protected:

	R* random;
	D* dsemi;
	double tuning;
	int n;
};

class DSemiConjugateMappingMove : public MCUpdate	{

	public:

	DSemiConjugateMappingMove(PhyloProcess* inprocess, AbstractPathConj* inpathconj, double intuning = 1) : process(inprocess), pathconj(inpathconj) {
		tuning  = intuning;
		initialized = 0;
	}

	double Move(double tuning_modulator=1)	{
		pathconj->InactivateSufficientStatistic();
		if (! initialized)	{
			process->Move(1.0);
			initialized = 1;
		}
		else	{
			process->Move(tuning);
		}
		pathconj->ActivateSufficientStatistic();
		return 1;
	}

	protected:

	PhyloProcess* process;
	AbstractPathConj* pathconj;
	double tuning;
	int initialized;
};

#endif


