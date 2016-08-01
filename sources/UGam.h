
#ifndef UGAMJITTER_H
#define UGAMJITTER_H

#include "MultiVariateTreeProcess.h"


class MeanJitteredLogit: public Dvar<UnitReal>{

	private:

		MultiNormal* up;
		MultiNormal* down;
		int index;
		Var<Real>* whitenoisemean;

	public:

	MeanJitteredLogit(MultiNormal* inup, MultiNormal* indown, int inindex, Var<Real>* inwhitenoisemean){
		up =inup;
		down = indown;
		index = inindex;
		Register(up);
		Register(down);
		if (whitenoisemean)	{
			Register(whitenoisemean);
		}
		specialUpdate();
	}

	void specialUpdate(){
		double x = 0.5 * ((*up)[index] + (*down)[index]);
		if (whitenoisemean)	{
			x += whitenoisemean->val();
		}
		double y = exp(x) / (1.0 + exp(x));
		setval(y);
	}
};

class MeanJitteredLogitTree: public BranchValPtrTree<Dvar<UnitReal> >	{

	public:

	MeanJitteredLogitTree(MultiVariateTreeProcess* inprocess, int inindex, BranchProcess<Real>* inwnprocess, bool withroot = true)	{
		process = inprocess;
		wnprocess = inwnprocess;
		index = inindex;
		SetWithRoot(withroot);
		RecursiveCreate(GetRoot());
	}

	~MeanJitteredLogitTree()	{
		RecursiveDelete(GetRoot());
	}

	Var<UnitReal>* GetRootRate() {return GetBranchVal(0);}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	double GetVar()	{
		int n =0;
		double tmp1 = RecursiveGetMeanSquare(GetRoot(),n);
		int m = 0;
		double tmp2 = RecursiveGetMean(GetRoot(),m);
		tmp1 /= n;
		tmp2 /= n;
		tmp1 -= tmp2 * tmp2;
		return tmp1;
	}

	double GetMean()	{
		int n =0;
		double tmp = RecursiveGetMean(GetRoot(),n);
		return tmp / n;
	}

	double GetTotal()	{
		int n =0;
		double tmp = RecursiveGetMean(GetRoot(),n);
		return tmp;
	}

	protected:

	double RecursiveGetMean(Link* from, int& n)	{
		double tmp = 0;
		if ((! from->isRoot()) || WithRoot())	{
			tmp += GetBranchVal(from->GetBranch())->val();
			n++;
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += RecursiveGetMean(link->Out(),n);
		}
		return tmp;
	}

	double RecursiveGetMeanSquare(Link* from, int& n)	{
		double tmp = 0;
		if ((! from->isRoot()) || WithRoot())	{
			double t = GetBranchVal(from->GetBranch())->val();
			tmp += t * t;
			n++;
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += RecursiveGetMeanSquare(link->Out(),n);
		}
		return tmp;
	}

	void specialUpdate(Link* from)	{
		if ((! from->isRoot()) || WithRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<UnitReal>* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new MeanJitteredLogit(process->GetMultiNormal(link), process->GetMultiNormal(link), index, 0);
		}
		return new MeanJitteredLogit(process->GetMultiNormal(link), process->GetMultiNormal(link->Out()), index, wnprocess->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return process->GetTree();}

	private:

	MultiVariateTreeProcess* process;
	BranchProcess<Real>* wnprocess;
	int index;
};

#endif
