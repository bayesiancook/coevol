
#ifndef SPLITCHRONO_H
#define SPLITCHRONO_H

#include "SplitLengthTree.h"

class SplitNodeDate : public Dvar<PosReal>	{

	public:

	SplitNodeDate(Var<PosReal>* inup, Var<PosReal>* indown, double inf)	{
		down = indown;
		up = inup;
		f = inf;
		Register(down);
		Register(up);
		specialUpdate();
	}

	void specialUpdate()	{
		double tmp = down->val() + (up->val() - down->val()) * f;
		if (isnan(tmp))	{
			cerr << "error in SplitNodeDate:: special update: nan\n";
			exit(1);
		}
		setval(tmp);
	}

	private:

	Var<PosReal>* up;
	Var<PosReal>* down;
	double f;
};



class SplitChronogram : public virtual SplitLengthTree, public virtual NodeBranchValPtrTree<Dvar<PosReal>, Dvar<PosReal> >	{

// NodeBranchValPtrTree<Dvar<PosReal> > {

	public:

	SplitChronogram(Chronogram* infrom, SplitTree* insplittree) : SplitLengthTree(infrom,insplittree) , chrono(infrom) {
		NodeValPtrTree<Dvar<PosReal> > ::RecursiveCreate(GetRoot());
	}

	/*
	virtual void Check()	{
		RecursiveCheck(GetRoot());
		cerr << "\n";
	}

	virtual void RecursiveCheck(const Link* from)	{

		if (from->isLeaf())	{
			cerr << "@" << from->GetNode()->GetName() << "@";
			cerr << GetNodeVal(from->GetNode())->val() ;
		}
		else	{
			if (from->Next()->Next() != from)	{
				cerr << '(';
			}
			for (const Link* link = from->Next(); link!=from; link=link->Next())	{
				RecursiveCheck(link->Out());
				if (link->Next() != from)	{
					cerr << ',';
				}
			}
			if (from->Next()->Next() != from)	{
				cerr << ')';
			}
		}
		cerr << ":" << GetNodeVal(from->GetNode())->val();
	}
	*/

	protected:

	void RecursiveUpdate(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveUpdate(link->Out());
		}
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		GetNodeVal(from->GetNode())->specialUpdate();
	}

	/*
	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		SplitBranch* splitbranch = splittree->GetSplitBranch(link->GetBranch());
		return new SplitLength(from->GetBranchVal(splitbranch->GetMother()), splitbranch->GetSplitFraction());
	}
	*/

	Dvar<PosReal>* CreateNodeVal(const Link* link)	{
		SplitNode* snode = splittree->GetSplitNode(link->GetNode());
		const Node* node1 = snode->GetUp();
		const Node* node2 = snode->GetDown();
		double f = snode->GetSplitFraction();
		SplitNodeDate* ret = new SplitNodeDate(chrono->GetNodeDate(node1),chrono->GetNodeDate(node2),f);
		return ret;
	}

	private:
	Chronogram* chrono;
};


#endif

