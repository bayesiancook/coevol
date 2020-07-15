#pragma once

#include "Chronogram.h"

class BranchLength : public Dvar<PosReal>{

	private:
	Var<PosReal>* Before;
	Var<PosReal>* After;
	Var<PosReal>* Rate;


	public:
	BranchLength(Var<PosReal>* inBefore, Var<PosReal>* inAfter, Var<PosReal>* inRate){
		Before = inBefore;
		After = inAfter;
		Rate = inRate;
		Register(inRate);
		Register(inBefore);
		Register(inAfter);
		SetName("time");
		specialUpdate();
	}


	~BranchLength(){};

	void specialUpdate(){
		double newval = (Before->val() - After->val()) * Rate->val();
        if (newval < 0) {
            cerr << "in branch time tree: negative time\n";
            cerr << newval << '\n';
            exit(1);
        }
		if (newval < 1e-7)	{
			newval = 1e-7;
		}
		setval(newval);
	}

};


class BranchTimeTree : public BranchValPtrTree<Dvar<PosReal> >  {

    public:

    BranchTimeTree(NodeVarTree<PosReal>* inchronogram, Var<PosReal>* inrate)    {
        SetWithRoot(false);
        chronogram = inchronogram;
        rate = inrate;
        RecursiveCreate(GetRoot());
        specialUpdate();
    }

    ~BranchTimeTree()   {
        RecursiveDelete(GetRoot());
    }

    Tree* GetTree() {
        return chronogram->GetTree();
    }

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	double GetTotalTime()	{
		return RecursiveGetTotalTime(GetRoot());
	}

    protected:

	double RecursiveGetTotalTime(const Link* from)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotalTime(link->Out());
		}
		if (! from->isRoot())	{
			total += GetBranchVal(from->GetBranch())->val();
		}
		return total;
	}

	void specialUpdate(Link* from)	{
		if ((! from->isRoot()) || WithRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new BranchLength(chronogram->GetNodeVal(link->GetNode()),
                chronogram->GetNodeVal(link->Out()->GetNode()),
                rate);
	}

	private:

	NodeVarTree<PosReal>* chronogram;
	Var<PosReal>* rate;
};
