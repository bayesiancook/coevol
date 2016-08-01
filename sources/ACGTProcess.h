
#ifndef ACGT_H
#define ACGT_H


class InstantStat : public Dvar<Profile>	{

	public:

	InstantStat(Var<RealVector>* inlogval, int indim, int inoffset)	{
		setval(Profile(indim));
		bkvalue = Profile(indim);
		logval = inlogval;
		offset = inoffset;
		Register(logval);
		specialUpdate();
	}

	void specialUpdate()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = exp((*logval)[i + offset]);
			total += (*this)[i];
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}

	private:

	Var<RealVector>* logval;
	int offset;

};

class BranchStat : public Dvar<Profile>	{

	public:

	BranchStat(Var<Profile>* inup, Var<Profile>* indown)	{
		setval(Profile(inup->GetDim()));
		bkvalue = Profile(inup->GetDim());
		up = inup;
		down = indown;
		Register(up);
		Register(down);
		specialUpdate();
	}

	void specialUpdate()	{
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = 0.5 * ( (*up)[i] + (*down)[i] );
		}
	}

	private:

	Var<Profile>* up;
	Var<Profile>* down;
};

class StatTree : public NodeValPtrTree<Dvar<Profile> >	{

	public:

	StatTree(NodeVarTree<RealVector>* inprocess, int indim, int inoffset)	{
		dim = indim;
		offset = inoffset;
		process = inprocess;
		RecursiveCreate(GetRoot());
	}

	~StatTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	double GetMeanEntropy()	{
		int n = 0;
		double total = GetTotalEntropy(GetRoot(),n);
		return total / n;
	}

	double GetVarEntropy()	{
		int n = 0;
		double total1 = GetTotalEntropy(GetRoot(),n);
		n = 0;
		double total2 = GetTotalSquareEntropy(GetRoot(),n);
		total1 /= n;
		total2 /= n;
		total2 -= total1 * total1;
		return total2;
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	void RecursiveSpecialUpdate(const Link* from)	{
		GetNodeVal(from->GetNode())->specialUpdate();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	protected:

	double GetTotalEntropy(const Link* from, int& n)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalEntropy(link->Out(),n);
		}
		total += GetNodeVal(from->GetNode())->GetEntropy();
		n++;
		return total;
	}

	double GetTotalSquareEntropy(const Link* from, int& n)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareEntropy(link->Out(),n);
		}
		double tmp = GetNodeVal(from->GetNode())->GetEntropy();
		total += tmp * tmp;
		n++;
		return total;
	}

	Dvar<Profile>* CreateNodeVal(const Link* link)	{
		return new InstantStat(process->GetNodeVal(link->GetNode()),dim,offset);
	}

	private:
	NodeVarTree<RealVector>* process;
	int dim;
	int offset;

};

class BranchStatTree : public BranchValPtrTree<Dvar<Profile> >	{

	public:

	BranchStatTree(NodeVarTree<Profile>* inprocess)	{
		process = inprocess;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~BranchStatTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	protected:

	void RecursiveSpecialUpdate(const Link* from)	{
		if (from->isRoot() && WithRoot())	{
			GetBranchVal(0)->specialUpdate();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->specialUpdate();
			RecursiveSpecialUpdate(link->Out());
		}
	}

	Dvar<Profile>* CreateBranchVal(const Link* link)	{
		return new BranchStat(process->GetNodeVal(link->GetNode()),process->GetNodeVal(link->Out()->GetNode()));
	}

	private:
	NodeVarTree<Profile>* process;

};

class NucMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	NucMatrixTree(Var<Profile>* inrelrate, BranchVarTree<Profile>* instattree) {
		SetWithRoot(true);
		stattree = instattree;
		relrate = inrelrate;
		RecursiveCreate(GetRoot());
	}

	~NucMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),true);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;
	Var<Profile>* relrate;

};

#endif
