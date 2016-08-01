#ifndef FIXEDLENGTHTREE_H
#define FIXEDLENGTHTREE_H

#include "ValTree.h"


class FixedLength : public Const<PosReal>	{

	public:

	FixedLength(Var<PosReal>* inmu, double l) : Const<PosReal>(l)	{
		mu = inmu;
		Register(mu);
	}

	protected:

	Var<PosReal>* mu;
};

class FixedLengthTree : public BranchValPtrTree <Dvar<PosReal> >	{

	public:

	FixedLengthTree(Tree* intree, Var<PosReal>* inmu) : tree(intree) {
		mu = inmu;
		SetWithRoot(false);
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree() {return tree;}

	/*
	bool comp (double a, double b) {return (a>b);}

	vector<double> GetDateList()	{
		vector<double> date;
		PushAges(GetRoot(),date,0);
		sort(date.begin(), date.end(), comp);
		date.push(0);
		for (unsigned int i=0; i<date.size(); i++)	{
			cerr << date[i] << '\n';
		}
		exit(1);
		return date;
	}

	void PushAges(const Link* from, vector<double>* date, double currentdate)	{
		if (! from->isLeaf())	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				PushAges(link->Out(), currentdate + GetBranchVal(link->GetBranch())->val());
			}
		}
	}
	*/

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{

		double l =  atof((link->GetBranch()->GetName()).c_str());
		if (l <= 0)	{
			cerr << "error in fixedlengthtree::createbranchval : " << l << '\n';
			exit(1);
		}
		return new FixedLength(mu, l);
	}

	protected:

	Var<PosReal>* mu;
	Tree* tree;
};

#endif
