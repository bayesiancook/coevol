

#ifndef SPLITLENGTHTREE_H
#define SPLITLENGTHTREE_H

#include "SplitTree.h"
#include "ValTree.h"

class SplitLength : public Dvar<PosReal>	{

	public:

	SplitLength(Var<PosReal>* infrom, double inf)	{
		from = infrom;
		f = inf;
		Register(from);
		specialUpdate();
	}

	void specialUpdate()	{
		double tmp = from->val() * f;
		if (std::isnan(tmp))	{
			cerr << "error in special update: nan\n";
			cerr << from->val() << '\t' << f << '\n';
			exit(1);
		}
		setval(from->val() * f);
	}

	private:

	Var<PosReal>* from;
	double f;
};

class SplitLengthTree : public virtual BranchValPtrTree<Dvar<PosReal> > {

	public:


	SplitLengthTree(LengthTree* infrom, SplitTree* insplittree)	{
		cerr << "in split length tree\n";
		SetWithRoot(false);
		splittree = insplittree;
		from = infrom;
		RecursiveCreate(GetRoot());
	}

	virtual ~SplitLengthTree() {}

	Tree* GetTree() {return splittree;}

	void specialUpdate()	{
		RecursiveUpdate(GetRoot());
	}

	virtual void Check()	{
		RecursiveCheck(GetRoot());
		cerr << "\n";
	}

	virtual void RecursiveCheck(const Link* from)	{

		if (from->isLeaf())	{
			cerr << "@" << from->GetNode()->GetName() << "@";
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
		if (! from->isRoot())	{
			cerr << ":" << GetBranchVal(from->GetBranch())->val();
		}
		else	{
			cerr << ";";
		}
	}
	protected:

	virtual void RecursiveUpdate(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveUpdate(link->Out());
		}
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		SplitBranch* splitbranch = splittree->GetSplitBranch(link->GetBranch());
		return new SplitLength(from->GetBranchVal(splitbranch->GetMother()), splitbranch->GetSplitFraction());
	}

	SplitTree* splittree;
	LengthTree* from;
};

#endif
