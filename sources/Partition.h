
#ifndef NODEPARTITION_H
#define NODEPARTITION_H


#include "ValArray.h"
#include "ValTree.h"
#include "IID.h"

#include <map>

using namespace std;

// associate an index in [0..K-1] to each branch of a tree
class BranchPartition	{

	public:

	// BranchPartition() {}

	BranchPartition(Tree* intree, bool all = false)	{
		tree = intree;
		if (all)	{
			NComponent = 0;
			AllReset(tree->GetRoot());
		}
		else	{
			NComponent  = 1;
			Reset(tree->GetRoot());
		}
	}

	BranchPartition(Tree* intree, string filename, bool recursive  = false)	{
		tree = intree;
		Reset(tree->GetRoot());
		ifstream is(filename.c_str());
		ReadFromStream(is);
		if (recursive)	{
			RecursiveSetAlloc(tree->GetRoot(),0);
		}
	}

	Tree* GetTree()	{
		return tree;
	}

	void RecursiveSetAlloc(const Link* from, int i)	{
		if (i)	{
			SetAlloc(from->GetBranch(),i);
			// cerr << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' << i << '\n';
		}
		int j = 0;
		if (! from->isRoot())	{
			j = GetAlloc(from->GetBranch());
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAlloc(link->Out(),j);
		}
	}

	void ReadFromStream(istream& is)	{
		is >> NComponent;
		int N;
		is >> N;
		for (int i=0; i<N; i++)	{
			string name1;
			string name2;
			int k;
			is >> name1 >> name2 >> k;
			cerr << name1 << '\t' << name2 << '\t' << k << '\n';
			SetAlloc(name1,name2,k);
		}
	}

	void Reset(const Link* from)	{
		allocmap[from->GetBranch()] = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			Reset(link->Out());
		}
	}

	void AllReset(const Link* from)	{
		allocmap[from->GetBranch()] = NComponent;
		NComponent++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			AllReset(link->Out());
		}
	}

	int GetAlloc(const Branch* branch)	{
		int tmp = allocmap[branch];
		if ((tmp < 0) || (tmp >= NComponent))	{
			cerr << "error in allocation " << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	void SetAlloc(const Branch* branch, int k)	{
		if ((k < 0) || (k >= NComponent))	{
			cerr << "error in allocation " << k << '\n';
			exit(1);
		}
		allocmap[branch] = k;
	}

	void SetAlloc(string name1, string name2, int k)	{
		const Link* link= tree->GetLCA(name1,name2);
		if (! link)	{
			cerr << "error : did not find last common ancestor of " << name1 << " and " << name2 << '\n';
			exit(1);
		}
		SetAlloc(link->GetBranch(),k);
	}

	int GetNComponent()	{
		return NComponent;
	}

	void Check(const Link* from)	{
		if (! from->isLeaf())	{
			cerr << '(';
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				Check(link->Out());
				if (link->Next() != from)	{
					cerr << ",";
				}
			}
			cerr << ")";
		}
		else	{
			cerr << from->GetNode()->GetName();
		}
		if (from->isRoot())	{
			cerr << ";\n";
		}
		else	{
			cerr  << ':' << GetAlloc(from->GetBranch());
		}
	}

	protected:

	map<const Branch*,int> allocmap;
	Tree* tree;
	int NComponent;

};


class ScaleMixTree : public BranchVarTree<PosReal> {

	public:

	ScaleMixTree(BranchPartition* inpartition, VarArray<PosReal>* incomplist)	{
		partition = inpartition;
		complist = incomplist;
		if (partition->GetNComponent() != complist->GetSize())	{
			cerr << "error in ScaleMixTree: non matching number of components : " << partition->GetNComponent() << '\t' << complist->GetSize() << '\n';
			exit(1);
		}
	}

	Var<PosReal>* GetComponentList(int k)	{
		return complist->GetVal(k);
	}

	Tree* GetTree()	{
		return partition->GetTree();
	}

	Var<PosReal>* GetBranchVal(const Branch* branch)	{
		return complist->GetVal(partition->GetAlloc(branch));
	}

	private:

	BranchPartition* partition;
	VarArray<PosReal>* complist;

};

class GammaMixTree : public LengthTree 	, public virtual MCMC {
// BranchVarTree<PosReal> {

	public:

	GammaMixTree(BranchPartition* inpartition, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool fixfirst = true)	{
		partition = inpartition;
		alpha = inalpha;
		beta = inbeta;
		componentlist = new GammaIIDArray(partition->GetNComponent(),alpha,beta);
		if (fixfirst)	{
			componentlist->GetVal(0)->ClampAt(1.0);
		}
	}

	Tree* GetTree()	{
		return partition->GetTree();
	}

	Var<PosReal>* GetBranchVal(const Branch* branch)	{
		return componentlist->GetVal(partition->GetAlloc(branch));
	}

	void Sample()	{
		componentlist->Sample();
	}

	void drawSample()	{
		componentlist->Sample();
	}

	double Move(double tuning = 1)	{
		return componentlist->Move(tuning);
	}

	double GetLogProb()	{
		return componentlist->GetLogProb();
	}

	int GetNComponent()	{
		return componentlist->GetSize();
	}

	Var<PosReal>* GetComponent(int k)	{
		return componentlist->GetVal(k);
	}

	friend ostream& operator<<(ostream& os, const GammaMixTree& a)	{
		os << *a.componentlist;
		return os;
	}

	friend istream& operator>>(istream& is, GammaMixTree& a)  {
		is >> *a.componentlist;
		return is;
	}

	private:

	BranchPartition* partition;
	GammaIIDArray* componentlist;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};


#endif
