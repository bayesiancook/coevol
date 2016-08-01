
#ifndef LOGN_H
#define LOGN_H

#include <sstream>

#include "BaseType.h"
#include "Var.h"
#include "ValTree.h"
#include "ContinuousData.h"

class LogNormalProcessInstantValue : public Rvar<Real>	{

	public:

	// LogNormalProcessInstantValue() {}

	LogNormalProcessInstantValue(Var<PosReal>* intau, Var<Real>* inup = 0, Var<PosReal>* intime = 0, bool inprec = true)	{

		tau = intau;
		prec = inprec;
		up = inup;
		time = intime;

		if (up)	{
			Register(up);
		}
		if (time)	{
			Register(time);
		}
		Register(tau);
		if (up) {
			SetName("normal process instant value\n");
		}
		else	{
			SetName("root normal process instant value\n");
		}
		Sample();
	}

	virtual void drawSample()	{
		if (up)	{
			DrawFiniteTime();
		}
		else	{
			DrawStationary();
		}
	}

	virtual double logProb()	{
		if (up)	{
			return GetFiniteTimeTransitionLogProb();
		}
		return GetStationaryLogProb();
	}

	virtual double GetFiniteTimeTransitionLogProb()	{
		double tmp = val() - up->val();
		if (prec)	{
			return 0.5 * log(tau->val()) - 0.5 * tau->val() * tmp * tmp / time->val();
		}
		return -0.5 * log(tau->val()) - 0.5 / tau->val() * tmp * tmp / time->val();
	}

	// normal process has no stationary distribution!
	// let us define a default: standard normal distribution
	virtual double GetStationaryLogProb()	{
		return 0;
		// return - val() * val();
	}

	virtual void DrawFiniteTime()	{
		if (prec)	{
			setval(Random::sNormal() * sqrt(time->val() / tau->val()) + up->val());
		}
		else	{
			setval(Random::sNormal() * sqrt(time->val() * tau->val()) + up->val());
		}
	}

	virtual void DrawStationary()	{
		setval(Random::sNormal());
	}

	protected:

	Var<PosReal>* tau;
	Var<PosReal>* time;
	Var<Real>* up;
	bool prec;
};


class LogNormalProcessFiniteTimeAverage : public Dvar<PosReal>	{

	public:

	LogNormalProcessFiniteTimeAverage(Var<Real>* inup, Var<Real>* indown, Var<PosReal>* intime, BranchValType inbval)	{
		bval = inbval;
		up = inup;
		down = indown;
		time = intime;

		Register(up);
		Register(down);
		Register(time);
		SetName("log normal process finite average");
		specialUpdate();
	}

	void specialUpdate()	{
		if (bval == INTEGRAL)	{
			value = (0.5 * (exp(up->val()) + exp(down->val()))) * time->val();
		}
		else	{
			value = (0.5 * (exp(up->val()) + exp(down->val())));
		}
	}

	public :

	Var<Real>* up;
	Var<Real>* down;
	Var<PosReal>* time;
	BranchValType bval;
};



class LogNormalTreeProcess : public NodeBranchValPtrTree<Rvar<Real>, Dvar<PosReal> >, public MCMC, public Additive	{

	public:

	LogNormalTreeProcess() {}

	LogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma, BranchValType inbval, bool inprec = true)  {
		SetWithRoot(false);
		bval = inbval;
		tree = intree;
		sigma = insigma;
		prec = inprec;
		cerr << "recursive create\n";
		RecursiveCreate(GetRoot());
		cerr << "ok\n";
	}

	~LogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = 0;
		if (! link->isRoot())	{
			time = GetLengthTree()->GetBranchLength(link->GetBranch());
		}

		// make the new instant value, and return the pointer
		return new LogNormalProcessInstantValue(sigma, vup, time, prec);
	}

	virtual Dvar<PosReal>* CreateBranchVal(const Link* link)	{

		if (link->isRoot())	{
			cerr << "error : create branch val called on root\n";
			exit(1);
		}
		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the instant value associated to this node
		Var<Real>* vdown = GetNodeVal(link->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new LogNormalProcessFiniteTimeAverage(vup, vdown, time, bval);
	}

	// various accessors
	virtual Tree* GetTree() {
		return GetLengthTree()->GetTree();
	}

	LengthTree* GetLengthTree()	{
		return tree;
	}

	Rvar<Real>* GetRootRate() {return GetNodeVal(GetRoot()->GetNode());}

	// various functions that basically iterate a simple action or computation over the tree
	virtual double Move(double tuning)	{
		int n = 0;
		double tot = Move(this->GetRoot(),tuning,n);
		return tot / n;
	}

	virtual double GetLogProb()	{
		return GetLogProb(this->GetRoot());
	}

	virtual void Clamp()	{
		Clamp(this->GetRoot());
	}

	virtual void drawSample()	{
		SampleNode(this->GetRoot());
	}

	double GetExpVal(const Link* link)	{
		return exp(GetNodeVal(link->GetNode())->val());
	}

	double GetTotalLength()	{

		int n = 0;
		double total = GetTotalRate(GetRoot(),n);
		return total;
	}

	double GetMeanRate()	{

		int n = 0;
		double total = GetTotalRate(GetRoot(),n);
		return total / n;
	}

	double GetVarRate()	{

		int n = 0;
		double total1 = GetTotalRate(GetRoot(),n);
		n = 0;
		double total2 = GetTotalSquareRate(GetRoot(),n);
		total1 /= n;
		total2 /= n;
		total2 -= total1 * total1;
		return total2;
	}

	double GetMeanLogRate()	{

		int n = 0;
		double total = GetTotalLogRate(GetRoot(),n);
		return total / n;
	}

	double GetVarLogRate()	{

		int n = 0;
		double total1 = GetTotalLogRate(GetRoot(),n);
		n = 0;
		double total2 = GetTotalSquareLogRate(GetRoot(),n);
		total1 /= n;
		total2 /= n;
		total2 -= total1 * total1;
		return total2;
	}

	void ClampRootAt(double d)	{
		GetNodeVal(GetRoot()->GetNode())->ClampAt(d);
	}

	void ClampAt(ContinuousData* contdata, int index = 0, bool logit=false)	{
		if (index >= contdata->GetNsite())	{
			cerr << "error in LogNormalTreepProcess::clampAt: index overflow\n";
			exit(1);
		}
		int tot = RecursiveClampAt(GetRoot(),contdata,index,logit);
		cerr << "total number of pop size data at leaves : " << tot << '\n';
	}

	int ScalarAddition(double d)	{
		RecursiveScalarAddition(GetRoot(),d);
		return 0;
	}

	void Register(DAGnode* in)  {
		RegisterNodeTree(GetRoot(),in);
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	void Translation(double u){
		RecursiveTranslation(this->GetRoot(),u);
	}

	protected:
	// recursive functions
	virtual void Clamp(Link* from)	{
		this->GetNodeVal(from->GetNode())->Clamp();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			Clamp(link->Out());
		}
	}

	virtual void SampleNode(Link* from)	{
		this->GetNodeVal(from->GetNode())->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			SampleNode(link->Out());
		}
	}

	virtual double Move(Link* from, double tuning, int& count)	{
		double total = this->GetNodeVal(from->GetNode())->Move(tuning);
		count++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += Move(link->Out(),tuning,count);
		}
		return total;
	}

	virtual double GetLogProb(Link* from)	{
		double total = this->GetNodeVal(from->GetNode())->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetLogProb(link->Out());
		}
		return total;
	}

	void RecursiveTranslation(const Link* from, double u)	{
		GetNodeVal(from->GetNode())->setval(GetNodeVal(from->GetNode())->val() + u);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTranslation(link->Out(),u);
		}
	}

	void RecursiveRegister(DAGnode* node, const Link* from)	{
		GetNodeVal(from->GetNode())->Register(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(node,link->Out());
		}
	}

	void RecursiveSpecialUpdate(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
			GetBranchVal(link->GetBranch())->specialUpdate();
		}
	}

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
		if (! GetNodeVal(from->GetNode())->isClamped())	{
			GetNodeVal(from->GetNode())->setval(0);
		}
	}

	void RecursiveScalarAddition(Link* from, double d)	{
		if (! GetNodeVal(from->GetNode())->isClamped())	{
			GetNodeVal(from->GetNode())->ScalarAddition(d);
		}
		else	{
			GetNodeVal(from->GetNode())->ScalarAddition(0);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveScalarAddition(link->Out(),d);
		}
	}

	int RecursiveClampAt(Link* from, ContinuousData* contdata, int index, bool logit)	{
		int tot = 0;
		if (from->isLeaf())	{
			int tax = contdata->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax != -1)	{
				double tmp = contdata->GetState(tax, index);
				if (tmp != -1)	{
					if (tmp <= 0)	{
						cerr << "error in recursive clamp at: negative cont data\n";
						exit(1);
					}
					if (logit)	{
						GetNodeVal(from->GetNode())->ClampAt(log(tmp / (1 - tmp)));
					}
					else	{
						GetNodeVal(from->GetNode())->ClampAt(log(tmp));
					}
				}
				tot++;
			}
		}
		else	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				tot += RecursiveClampAt(link->Out(),contdata,index,logit);
			}
		}
		return tot;
	}

	double GetTotalRate(const Link* from, int& n)	{
		double total = exp(GetNodeVal(from->GetNode())->val());
		n++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalRate(link->Out(),n);
		}
		return total;
	}

	double GetTotalSquareRate(const Link* from, int& n)	{
		double temp =  exp(GetNodeVal(from->GetNode())->val());
		double total = temp * temp;
		n++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareRate(link->Out(),n);
		}
		return total;
	}

	double GetTotalLogRate(const Link* from, int& n)	{
		double total = GetNodeVal(from->GetNode())->val();
		n++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalLogRate(link->Out(),n);
		}
		return total;
	}

	double GetTotalSquareLogRate(const Link* from, int& n)	{
		double temp =  GetNodeVal(from->GetNode())->val();
		double total = temp * temp;
		n++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareLogRate(link->Out(),n);
		}
		return total;
	}

	LengthTree* tree;
	Var<PosReal>* sigma;
	BranchValType bval;
	bool prec;
};

#endif

