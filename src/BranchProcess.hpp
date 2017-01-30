
#ifndef BRANCHPROCESS_H
#define BRANCHPROCESS_H

#include <sstream>

#include "ValTree.h"
#include "RandomTypes.h"

// iid random variables of type Rvar<V>
// associated to each branch of the tree

template<class V> class BranchProcess : public MCMC , public virtual BranchValPtrTree< Rvar<V> > {

	public:

	BranchProcess(Tree* intree, bool inwithroot = false) : tree(intree) {
		this->SetWithRoot(inwithroot);
	}

	Tree* GetTree() {return tree;}

	V val(const Branch* branch)	{
		return this->GetBranchVal(branch)->val();
	}

	void setval(const Branch* branch, V inval)	{
		this->GetBranchVal(branch)->setval(inval);
	}

	virtual void drawSample()	{
		drawSample(this->GetRoot());
	}

	virtual double Move(double tuning)	{
		int n = 0;
		double tot = Move(this->GetRoot(),tuning,n);
		return tot / n;
	}

	virtual double GetLogProb()	{
		return GetLogProb(this->GetRoot());
	}

	void RecursiveRegister(DAGnode* node, const Link* from)	{
		if (this->WithRoot() || (! from->isRoot()))	{
			this->GetBranchVal(from->GetBranch())->Register(node);
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(node,link->Out());
		}
	}

	protected:

	virtual void drawSample(const Link* from)	{
		if (this->WithRoot() && (from->isRoot()))	{
			this->GetBranchVal(0)->Sample();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			this->GetBranchVal(link->GetBranch())->Sample();
			drawSample(link->Out());
		}
	}


	virtual double Move(const Link* from, double tuning, int& count)	{
		double total = 0;
		if (this->WithRoot() && (from->isRoot()))	{
			total += this->GetBranchVal(0)->Move(tuning);
			count++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->Move(tuning);
			count++;
			total += Move(link->Out(),tuning,count);
		}
		return total;
	}

	virtual double GetLogProb(const Link* from)	{
		double total = 0;
		if (this->WithRoot() && (from->isRoot()))	{
			total += this->GetBranchVal(0)->GetLogProb();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->GetLogProb();
			total += GetLogProb(link->Out());
		}
		return total;
	}

	Tree* tree;
};


class GammaTree : public BranchProcess<PosReal> {

	public:

	GammaTree(Tree* intree)	: BranchProcess<PosReal>(intree) {}

	GammaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool inwithroot = false, bool inmeanvar = false)  : BranchProcess<PosReal>(intree) {
		SetWithRoot(inwithroot);
		meanvar = inmeanvar;
		alpha = inalpha;
		beta = inbeta;
		RecursiveCreate(GetRoot());
	}

	~GammaTree()	{
		RecursiveDelete(GetRoot());
	}

	Var<PosReal>* GetAlpha()	{
		return alpha;
	}

	Var<PosReal>* GetBeta()	{
		return beta;
	}

	double GetTotal()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total;

	}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	protected:

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new Gamma(alpha,beta,meanvar);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	bool meanvar;

};

class BetaTree : public BranchProcess<UnitReal> {

	public:

	BetaTree(Tree* intree)	: BranchProcess<UnitReal>(intree) {}

	BetaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool inwithroot = false)  : BranchProcess<UnitReal>(intree) {
		SetWithRoot(inwithroot);
		alpha = inalpha;
		beta = inbeta;
		RecursiveCreate(GetRoot());
	}

	~BetaTree()	{
		RecursiveDelete(GetRoot());
	}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	protected:

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<UnitReal>* CreateBranchVal(const Link* link)	{
		return new Beta(alpha,beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;

};

class ExponentialTree : public BranchProcess<PosReal> {

	public:

	ExponentialTree(Tree* intree, Var<PosReal>* inmeanlength) : BranchProcess<PosReal>(intree) {
		meanlength = inmeanlength;
		RecursiveCreate(GetRoot());
	}

	~ExponentialTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		Exponential* expo = new Exponential(meanlength,Exponential::MEAN);
		return expo;
	}

	private:

	Var<PosReal>* meanlength;

};

/*
template <class V> class IIDBranchProcess : public BranchProcess<V>, public IIDBranchVarTree<Rvar<V> >	{

	IIDBranchProcess(Tree* intree, Rvar<V>& from) : BranchProcess<V>(intree), IIDBranchVarTree<Rvar<V> >(from) {}
};

class GammaTree : public IIDBranchProcess<PosReal> {

	public:

	GammaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDBranchProcess<PosReal>(intree, Gamma(inalpha,inbeta)) {}
};

class ExponentialTree : public BranchProcess<PosReal> {

	public:

	ExponentialTree(Tree* intree, Var<PosReal>* inmeanlength) : IIDBranchProcess<PosReal>(intree, Exponential(inmeanlength,Exponential::Mean)) {}
};
*/


#endif // BRANCHPROCESS_H
