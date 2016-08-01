
#ifndef BRANCHPROD_H
#define BRANCHPROD_H

class BranchProduct : public Dvar<PosReal>	{

	public:

	BranchProduct(Var<PosReal>* inval1, Var<PosReal>* inval2, Var<PosReal>* inval3 = 0, bool ininactivate = false)	{
		inactivate = ininactivate;
		val1 = inval1;
		val2 = inval2;
		val3 = inval3;
		Register(val1);
		Register(val2);
		Register(val3);
		specialUpdate();
	}

	~BranchProduct() {}

	void specialUpdate()	{
		if (inactivate)	{
			setval(0);
		}
		else	{
			if (val3)	{
				setval(val1->val() * val2->val() * val3->val());
			}
			else	{
				setval(val1->val() * val2->val());
			}
		}
	}

	protected:

	Var<PosReal>* val1;
	Var<PosReal>* val2;
	Var<PosReal>* val3;

	bool inactivate;

};

class BranchProductProcess : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	BranchProductProcess(BranchVarTree<PosReal>* inprocess1, BranchVarTree<PosReal>* inprocess2, Var<PosReal>* inoffset = 0, bool ininactivate = false)	{
		inactivate = ininactivate;
		process1 = inprocess1;
		process2 = inprocess2;
		offset = inoffset;
		SetWithRoot(false);
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree()	{
		return process1->GetTree();
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	BranchProduct* GetBranchProduct(const Branch* branch)	{
		return dynamic_cast<BranchProduct*>(GetBranchVal(branch));
	}

	protected:

	void RecursiveSpecialUpdate(const Link* from)	{
		if (! from->isRoot())	{
			GetBranchProduct(from->GetBranch())->specialUpdate();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	BranchProduct* CreateBranchVal(const Link* link)	{
		return new BranchProduct(process1->GetBranchVal(link->GetBranch()),process2->GetBranchVal(link->GetBranch()),offset,inactivate);
	}

	BranchVarTree<PosReal>* process1;
	BranchVarTree<PosReal>* process2;
	Var<PosReal>* offset;
	bool inactivate;

};

#endif
