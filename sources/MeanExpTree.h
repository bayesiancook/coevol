

#ifndef MEANEXPTREE_H
#define MEANEXPTREE_H

#include "MultiVariateTreeProcess.h"


class MeanExp: public Dvar<PosReal>{

	private:

		Var<PosReal>* time;
		Var<Real>* up;
		Var<Real>* down;
		Var<Real>* offset;
		BranchValType bval;

	public:

	MeanExp(Var<Real>* inup, Var<Real>* indown, Var<PosReal>* intime, BranchValType inbval = INTEGRAL, Var<Real>* inoffset = 0){
		up =inup;
		down = indown;
		time = intime;
		bval = inbval;
		offset = inoffset;
		Register(up);
		Register(down);
		Register(time);
		Register(offset);
		specialUpdate();
	}

	void specialUpdate(){
		double tmpoffset = 1;
		if (offset)	{
			tmpoffset = exp(offset->val());
		}
		if (time)	{
			if (bval == INTEGRAL)	{
				setval( (exp(up->val()) + exp(down->val()))/2 * time->val() * tmpoffset);
			}
			else	{
				setval( (exp(up->val()) + exp(down->val()))/2 * tmpoffset);
			}
		}
		else	{
			if (bval == INTEGRAL)	{
				cerr << "error in meanexp\n";
				exit(1);
			}
			else	{
				setval(tmpoffset);
			}
		}
	}
};

class MeanExpTree : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	MeanExpTree(NodeVarTree<Real>* inprocess, LengthTree* inlengthtree, BranchValType inbval = INTEGRAL, Var<Real>* inoffset = 0)	{
		SetWithRoot(true);
		process = inprocess;
		lengthtree = inlengthtree;
		bval = inbval;
		offset = inoffset;
		RecursiveCreate(GetRoot());
	}

	~MeanExpTree()	{
		RecursiveDelete(GetRoot());
	}

	Var<PosReal>* GetRootRate() {return GetBranchVal(0);}

	void specialUpdate()	{
		specialUpdate(GetRoot());
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
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += GetBranchVal(link->GetBranch())->val();
			tmp += RecursiveGetMean(link->Out(),n);
			n++;
		}
		return tmp;
	}

	void specialUpdate(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->specialUpdate();
			specialUpdate(link->Out());
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			if (bval == INTEGRAL)	{
				return 0;
			}
			else	{
				return new MeanExp(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), lengthtree->GetBranchVal(link->GetBranch()), bval, offset);
			}
		}
		return new MeanExp(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), lengthtree->GetBranchVal(link->GetBranch()), bval, offset);
	}

	Tree* GetTree() {return process->GetTree();}

	private:

	NodeVarTree<Real>* process;
	Var<Real>* offset;
	LengthTree* lengthtree;
	BranchValType bval;
};


class MeanExpFromMultiVariate: public Dvar<PosReal>{

	protected:

		Var<PosReal>* time;
		Var<RealVector>* up;
		Var<RealVector>* down;
		Var<PosReal>* jitter;
		int index;
		BranchValType bval;
		bool meanexp;

	public:

	MeanExpFromMultiVariate() : jitter(0) {
	}

	MeanExpFromMultiVariate(Var<RealVector>* inup, Var<RealVector>* indown, Var<PosReal>* intime, int inindex, BranchValType inbval, bool inmeanexp, Var<PosReal>* injitter = 0){

		jitter = injitter;
		up =inup;
		down = indown;
		time = intime;
		index = inindex;
		bval = inbval;
		meanexp = inmeanexp;
		Register(up);
		Register(down);
		if (time)	{
			Register(time);
		}
		else	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanExp : null time\n";
				exit(1);
			}
		}
		Register(jitter);
		specialUpdate();
	}

	void specialUpdate(){
		if (meanexp)	{
			double mean = 0;
			double& x1 = (*up)[index];
			double& x2 = (*down)[index];
			if (fabs(x2 - x1) > 1e-6)	{
				mean = (exp(x2) - exp(x1)) / (x2 - x1);
			}
			else	{
				mean = exp((x1 + x2) / 2);
			}

			if (jitter)	{
				mean *= jitter->val();
			}

			if (bval == INTEGRAL)	{
				setval(mean * time->val());
			}
			else	{
				setval(mean);
			}
		}
		else	{
			double f = 1;
			if (jitter)	{
				f = jitter->val();
			}
			if (bval == INTEGRAL)	{
				setval( (exp(up->val()[index]) + exp(down->val()[index]))/2 * time->val() * f);
			}
			else	{
				setval( (exp(up->val()[index]) + exp(down->val()[index]))/2 * f);
			}
		}
	}

};

class MeanExpTreeFromMultiVariate : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	MeanExpTreeFromMultiVariate()	{
	}

	MeanExpTreeFromMultiVariate(LengthTree* inlengthtree, NodeVarTree<RealVector>* inprocess, int inindex, BranchValType inbval, bool withroot, bool inmeanexp, BranchVarTree<PosReal>* injittertree = 0)	{
		process = inprocess;
		jittertree = injittertree;
		lengthtree = inlengthtree;
		index = inindex;
		bval = inbval;
		meanexp = inmeanexp;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	MeanExpTreeFromMultiVariate(MultiVariateTreeProcess* inprocess, int inindex, BranchValType inbval, bool withroot, bool inmeanexp, BranchVarTree<PosReal>* injittertree = 0)	{
		process = inprocess;
		jittertree = injittertree;
		lengthtree = inprocess->GetLengthTree();
		index = inindex;
		bval = inbval;
		meanexp = inmeanexp;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~MeanExpTreeFromMultiVariate()	{
		RecursiveDelete(GetRoot());
	}

	LengthTree* GetLengthTree()	{
		return lengthtree;
	}

	Var<PosReal>* GetRootRate() {return GetBranchVal(0);}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	double GetMax()	{
		double tmp = RecursiveGetMax(GetRoot());
		return tmp;
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

	double RecursiveGetMax(Link* from)	{
		double ret = 0;
		if ((! from->isRoot()) || WithRoot())	{
			ret = GetBranchVal(from->GetBranch())->val();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMax(link->Out());
			if (ret < tmp)	{
				ret = tmp;
			}
		}
		return ret;
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
		if (link->isRoot())	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanExpFromMult::CreateBranchVal\n";
				exit(1);
			}
			return new MeanExpFromMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->GetNode()), 0, index, bval, meanexp, 0);
		}
		Var<PosReal>* tmp = 0;
		if (jittertree)	{
			tmp = jittertree->GetBranchVal(link->GetBranch());
		}
		return new MeanExpFromMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), GetLengthTree()->GetBranchVal(link->GetBranch()), index, bval, meanexp, tmp);
	}

	Tree* GetTree() {return process->GetTree();}

	NodeVarTree<RealVector>* process;
	LengthTree* lengthtree;
	int index;
	BranchValType bval;
	bool meanexp;
	BranchVarTree<PosReal>* jittertree;
};

class MeanExpFromTwoMultiVariate: public MeanExpFromMultiVariate {

	protected:

	Var<Real>* alpha1;
	Var<Real>* alpha2;
	Var<PosReal>* beta1;
	Var<PosReal>* beta2;
	Var<PosReal>* beta3;
	int index2;

	public:

	MeanExpFromTwoMultiVariate(Var<RealVector>* inup, Var<RealVector>* indown, Var<PosReal>* intime, Var<Real>* inalpha1, Var<Real>* inalpha2, Var<PosReal>* inbeta1, Var<PosReal>* inbeta2, Var<PosReal>* inbeta3, int inindex, int inindex2, BranchValType inbval, bool inmeanexp){
		up =inup;
		down = indown;
		time = intime;
		index = inindex;
		index2 = inindex2;
		alpha1 = inalpha1;
		alpha2 = inalpha2;
		beta1 = inbeta1;
		beta2 = inbeta2;
		beta3 = inbeta3;
		bval = inbval;
		meanexp = inmeanexp;
		Register(up);
		Register(down);
		Register(alpha1);
		Register(alpha2);
		Register(beta1);
		Register(beta2);
		Register(beta3);
		if (time)	{
			Register(time);
		}
		else	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanExp : null time\n";
				exit(1);
			}
		}
		specialUpdate();
	}

	void specialUpdate(){
		double mean = 0;
		double x1 = (*up)[index];
		double y1 = (*down)[index];
		double x2 = (*up)[index2];
		double y2 = (*down)[index2];

		if (beta3)	{
			double ex1 = exp(x1);
			double ex2 = exp(x2);
			double ey1 = exp(y1);
			double ey2 = exp(y2);
			if (meanexp && (fabs(y2 - x2) > 1e-6))	{
				double q1 = (ey1 - ex1)/(y1 - x1);
				double q2 = (ey1*ey2 - ex1*ex2)/(y1+y2-x1-x2);
				double q3 = (ey1/ey2 - ex1/ex2)/(y1-y2-x1+x2);
				mean = beta1->val() * q1 + beta2->val() * q2 + beta3->val() * q3;
			}
			else	{
				double zx = ex1 * (beta1->val() + beta2->val() * ex2 + beta3->val() / ex2);
				double zy = ey1 * (beta1->val() + beta2->val() * ey2 + beta3->val() / ey2);
				mean = 0.5 * (zx + zy);
			}
		}
		else	{
			double x = alpha1->val() * x1 + alpha2->val() * x2;
			double y = alpha1->val() * y1 + alpha2->val() * y2;

			if (meanexp && (fabs(y - x) > 1e-6))	{
				mean = (exp(y) - exp(x)) / (y - x);
			}
			else	{
				mean = exp((y + x) / 2);
			}

		}

		if (bval == INTEGRAL)	{
			setval(mean * time->val());
		}
		else	{
			setval(mean);
		}
	}

};

class MeanExpTreeFromTwoMultiVariate : public MeanExpTreeFromMultiVariate	{

	public:

	MeanExpTreeFromTwoMultiVariate(LengthTree* inlengthtree, NodeVarTree<RealVector>* inprocess, int inindex, int inindex2, Var<Real>* inalpha1, Var<Real>* inalpha2, Var<PosReal>* inbeta1, Var<PosReal>* inbeta2, Var<PosReal>* inbeta3, BranchValType inbval, bool withroot, bool inmeanexp)	{
		process = inprocess;
		lengthtree = inlengthtree;
		index = inindex;
		index2 = inindex2;
		alpha1 = inalpha1;
		alpha2 = inalpha2;
		beta1 = inbeta1;
		beta2 = inbeta2;
		beta3 = inbeta3;

		bval = inbval;
		meanexp = inmeanexp;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	MeanExpTreeFromTwoMultiVariate(MultiVariateTreeProcess* inprocess, int inindex, int inindex2, Var<Real>* inalpha1, Var<Real>* inalpha2, Var<PosReal>* inbeta1, Var<PosReal>* inbeta2, Var<PosReal>* inbeta3, BranchValType inbval, bool withroot, bool inmeanexp)	{
		process = inprocess;
		lengthtree = inprocess->GetLengthTree();
		index = inindex;
		index2 = inindex2;
		alpha1 = inalpha1;
		alpha2 = inalpha2;
		beta1 = inbeta1;
		beta2 = inbeta2;
		beta3 = inbeta3;
		bval = inbval;
		meanexp = inmeanexp;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~MeanExpTreeFromTwoMultiVariate()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanExpFromMult::CreateBranchVal\n";
				exit(1);
			}
			return new MeanExpFromTwoMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->GetNode()), 0, alpha1, alpha2, beta1, beta2, beta3, index, index2, bval, meanexp);
		}
		return new MeanExpFromTwoMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), GetLengthTree()->GetBranchVal(link->GetBranch()), alpha1, alpha2, beta1, beta2, beta3, index, index2, bval, meanexp);
	}

	int index2;
	Var<Real>* alpha1;
	Var<Real>* alpha2;
	Var<PosReal>* beta1;
	Var<PosReal>* beta2;
	Var<PosReal>* beta3;
};

class MeanLogitFromMultiVariate: public Dvar<UnitReal>{

	private:

		Var<PosReal>* time;
		Var<RealVector>* up;
		Var<RealVector>* down;
		int index;
		BranchValType bval;

	public:

	MeanLogitFromMultiVariate(Var<RealVector>* inup, Var<RealVector>* indown, Var<PosReal>* intime, int inindex, BranchValType inbval = INTEGRAL){
		up =inup;
		down = indown;
		time = intime;
		index = inindex;
		bval = inbval;
		Register(up);
		Register(down);
		if (time)	{
			Register(time);
		}
		else	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanLogit : null time\n";
				exit(1);
			}
		}
		specialUpdate();
	}

	void specialUpdate(){
		double xup = exp((*up)[index]) / (1 + exp((*up)[index]));
		double xdown = exp((*down)[index]) / (1 + exp((*down)[index]));
		if (bval == INTEGRAL)	{
			setval((xup + xdown)/2 * time->val());
		}
		else	{
			setval((xup + xdown)/2);
		}
	}
};

class MeanLogitTreeFromMultiVariate : public BranchValPtrTree<Dvar<UnitReal> >	{

	public:

	MeanLogitTreeFromMultiVariate(LengthTree* inlengthtree, NodeVarTree<RealVector>* inprocess, int inindex, BranchValType inbval, bool withroot)	{
		process = inprocess;
		lengthtree = inlengthtree;
		index = inindex;
		bval = inbval;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	MeanLogitTreeFromMultiVariate(MultiVariateTreeProcess* inprocess, int inindex, BranchValType inbval, bool withroot)	{
		process = inprocess;
		lengthtree = inprocess->GetLengthTree();
		index = inindex;
		bval = inbval;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~MeanLogitTreeFromMultiVariate()	{
		RecursiveDelete(GetRoot());
	}

	Var<UnitReal>* GetRootRate() {return GetBranchVal(0);}

	LengthTree* GetLengthTree()	{
		return lengthtree;
	}

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
			if (bval == INTEGRAL)	{
				cerr << "error in MeanLogitFromMult::CreateBranchVal\n";
				exit(1);
			}
			return new MeanLogitFromMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->GetNode()), 0, index, bval);
		}
		return new MeanLogitFromMultiVariate(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), GetLengthTree()->GetBranchVal(link->GetBranch()), index, bval);
	}

	Tree* GetTree() {return process->GetTree();}

	private:

	NodeVarTree<RealVector>* process;
	LengthTree* lengthtree;
	int index;
	BranchValType bval;
};

class Ratio: public Dvar<PosReal>{

	private:

		Var<PosReal>* numerator;
		Var<PosReal>* denominator;

	public:

	Ratio(Var<PosReal>* innumerator, Var<PosReal>* indenominator){
		numerator = innumerator;
		denominator = indenominator;
		Register(numerator);
		Register(denominator);
		specialUpdate();
	}

	void specialUpdate(){
		setval(numerator->val() / denominator->val());
	}
};

class RatioTree : public BranchValPtrTree<Dvar<PosReal> >	{

	private:

	BranchVarTree<PosReal>* numtree;
	BranchVarTree<PosReal>* dentree;
	public:

	RatioTree(BranchVarTree<PosReal>* innumtree, BranchVarTree<PosReal>* indentree)	{
		SetWithRoot(false);
		numtree = innumtree;
		dentree = indentree;
		RecursiveCreate(GetRoot());
	}

	~RatioTree()	{
		RecursiveDelete(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
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
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += GetBranchVal(link->GetBranch())->val();
			tmp += RecursiveGetMean(link->Out(),n);
			n++;
		}
		return tmp;
	}

	void specialUpdate(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			GetBranchVal(link->GetBranch())->specialUpdate();
			specialUpdate(link->Out());
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new Ratio(numtree->GetBranchVal(link->GetBranch()), dentree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return numtree->GetTree();}

};

#endif

