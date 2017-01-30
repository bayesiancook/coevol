
#ifndef ABSTRACTTREE_H
#define ABSTRACTTREE_H

#include "Tree.h"

class AbstractTree	{

	public:
	virtual ~AbstractTree() {}

	virtual Tree* GetTree() = 0;
	Link* GetRoot() {return GetTree()->GetRoot();}
};

template <class V> class BranchVarTree : public virtual AbstractTree	{

	public:
	virtual Var<V>* GetBranchVal(const Branch* branch) = 0;

	void SetName(string inname)	{
		RecursiveSetName(this->GetRoot(),inname);
	}

	void PrintLengths()	{
		RecursivePrintLengths(this->GetRoot());
	}

	/*
	int GetNBranchVals()	{
		int index = 0;
		RecursiveGetNBranchVals(GetRoot(),index);
		return index;
	}

	void GetBranchVals(V* array, int& index)	{
		RecursiveGetBranchVals(GetRoot(),array,index);
	}

	void SetBranchVals(V* array, int& index)	{
		RecursiveSetBranchVals(GetRoot(),array,index);
	}
	*/

	protected:

	/*
	void RecursiveGetNBranchVals(const Link* from, int& index);
		if (WithRoot() || (! from->isRoot()))	{
			index++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetNBranchVals(link->Out(),index);
		}
	}

	void RecursiveGetBranchVals(const Link* from, V* array, int& index)	{
		if (WithRoot() || (! from->isRoot()))	{
			array[index++] = GetBranchVal(from->GetBranch())->val();
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetBranchVals(link->Out(),array,index);
		}
	}

	void RecursiveSetBranchVals(const Link* from, V* array, int& index)	{
		if (WithRoot() || (! from->isRoot()))	{
			GetBranchVal(from->GetBranch())->setval(array[index++]);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetBranchVals(link->Out(),array,index);
		}
	}
	*/

	void RecursiveSetName(Link* from, string inname)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->SetName(inname);
			RecursiveSetName(link->Out(),inname);
		}
	}

	void RegisterBranchTree(Link* from, DAGnode* in)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->Register(in);
			RegisterBranchTree(link->Out(),in);
		}
	}

	void RecursivePrintLengths(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			ostringstream s;
			s << GetBranchVal(link->GetBranch())->val();
			link->GetBranch()->SetName(s.str());
			RecursivePrintLengths(link->Out());
		}
	}
};

template <class V> class NodeVarTree : public virtual AbstractTree	{

	public:
	virtual Var<V>* GetNodeVal(const Node* node) = 0;

	Var<V>* GetRootVal() {return GetNodeVal(GetRoot()->GetNode());}

	virtual void RegisterNodeTree(DAGnode* in)	{
		RegisterNodeTree(GetRoot(),in);
	}

	/*
	int GetNNodeVals()	{
		int index = 0;
		RecursiveGetNNodeVals(GetRoot(),index);
		return index;
	}

	void GetNodeVals(V* array, int& index)	{
		RecursiveGetNodeVals(GetRoot(),array,index);
	}

	void SetNodeVals(V* array, int& index)	{
		RecursiveSetBranchVals(GetRoot(),array,index);
	}
	*/

	protected:

	/*
	void RecursiveGetNNodeVals(const Link* from, int& index);
		index++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetNNodeVals(link->Out(),index);
		}
	}

	void RecursiveGetNodeVals(const Link* from, V* array, int& index)	{
		GetNodeVal(from->GetNode())->GetVals(array,index);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetNodeVals(link->Out(),array,index);
		}
	}

	void RecursiveSetNodeVals(const Link* from, V* array, int& index)	{
		GetNodeVal(from->GetNode())->setval(array[index++]);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetNodeVals(link->Out(),array,index);
		}
	}
	*/

	void RegisterNodeTree(Link* from, DAGnode* in)	{
		GetNodeVal(from->GetNode())->Register(in);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RegisterNodeTree(link->Out(),in);
		}
	}
};
//
// describes a tree with a Var<PosReal> (a branch length)
// associated to each branch
//
class LengthTree	: public virtual BranchVarTree<PosReal> {

	public:

	~LengthTree() {}

	virtual Var<PosReal>* GetBranchLength(const Branch* branch)	{
		return GetBranchVal(branch);
	};

	void SetBranchLengths()	{
		SetBranchLengths(GetRoot());
	}

	double GetTotalLength()	{
		return GetTotalLength(GetRoot());
	}

	protected:

	double GetTotalLength(const Link* from)	{
		double tot = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += GetBranchLength(link->GetBranch())->val();
			tot += GetTotalLength(link->Out());
		}
		return tot;
	}

	void SetBranchLengths(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchLength(link->GetBranch())->setval(atof((link->GetBranch()->GetName()).c_str()));
			SetBranchLengths(link->Out());
		}
	}

	public:

	void GetMeanAndVar(double& mean, double& var)	{
		mean = 0;
		var = 0;
		int count = 0;
		RecursiveGetMeanAndVar(GetRoot(),mean,var,count);
		mean/=count;
		var/=count;
		var-=mean*mean;
	}

	double GetMin()	{
		return RecursiveGetMin(GetRoot());
	}

	double GetMax()	{
		return RecursiveGetMax(GetRoot());
	}

	void GetFractionAbove(double cutoff, int& n, int& totjoint)	{

		double mean = 0;
		double var = 0;
		GetMeanAndVar(mean,var);
		n = 0;
		totjoint = 0;
		RecursiveGetFractionAbove(GetRoot(), GetRoot(), cutoff*mean, n, totjoint);
	}

	protected:

	double RecursiveGetMin(const Link* from)	{

		double min = -1;
		if (! from->isRoot())	{
			min = GetBranchVal(from->GetBranch())->val();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMin(link->Out());
			if ((min == -1) || (min > tmp))	{
				min = tmp;
			}
		}
		return min;
	}

	double RecursiveGetMax(const Link* from)	{

		double max = 0;
		if (! from->isRoot())	{
			max = GetBranchVal(from->GetBranch())->val();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMax(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		return max;
	}

	void RecursiveGetMeanAndVar(const Link* from, double& mean, double& var, int& count)	{
		if (! from->isRoot())	{
			double tmp = GetBranchVal(from->GetBranch())->val();
			mean += tmp;
			var += tmp*tmp;
			count++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetMeanAndVar(link->Out(),mean,var,count);
		}
	}

	void RecursiveGetFractionAbove(const Link* fromup, const Link* from, double value, int& n, int& totjoint)	{
		if (! from->isRoot())	{
			double tmp = GetBranchVal(from->GetBranch())->val();
			if (tmp > value)	{
				n++;
				if (! fromup->isRoot())	{
					double tmpup = GetBranchVal(fromup->GetBranch())->val();
					if (tmpup > value)	{
						totjoint ++;
					}
				}
			}

		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetFractionAbove(from,link->Out(),value,n,totjoint);
		}
	}

	public:

	void CheckRootToTip(const Link* from, double p = 0)	{
		if (from->isLeaf())	{
			cerr << p  << '\t';
		}
		else	{
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				CheckRootToTip(link->Out(), p + GetBranchLength(link->GetBranch())->val());
			}
		}
	}
};

class RandomLengthTree	: public virtual LengthTree , public virtual Multiplicative {

	public:

	virtual void Register(DAGnode* in)	{
		RegisterBranchTree(GetRoot(),in);
	}

	int ScalarMultiplication(double d)	{
		return ScalarMultiplication(GetRoot(),d);
	}

	protected:

	int ScalarMultiplication(const Link* from, double d)	{
		int tot = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += GetBranchLength(link->GetBranch())->ScalarMultiplication(d);
			tot += ScalarMultiplication(link->Out(),d);
		}
		return tot;
	}
};


template <class U, class V> class NodeBranchVarTree : public virtual NodeVarTree<U>, public virtual BranchVarTree<V> {};

template <> class NodeBranchVarTree<PosReal,PosReal> : public virtual NodeVarTree<PosReal>, public virtual LengthTree {};
template <class U> class NodeBranchVarTree<U,PosReal> : public virtual NodeVarTree<U>, public virtual LengthTree {};

#endif

