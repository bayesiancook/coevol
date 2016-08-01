#ifndef SCALABLENORMALTREEPROCESS_H
#define SCALABLENORMALTREEPROCESS_H

#include "AutocorrelatedProcess.h"
#include "ContinuousData.h"

// implements the random variable corresponding to the instant value of the process
// (say, at a node of a tree)
//
// implements the interface MarkovProcessInstantValue<V>
// with V=Real (Normal process takes real values)
//

class NormalProcessInstantValue : public MarkovProcessInstantValue<Real> {

	public:

	// up is the instant value at time 0
	// "this" is the instant value at time "time"
	// both are Var<Real> (normal process)
	//
	NormalProcessInstantValue(Var<PosReal>* insigma, Var<Real>* inup = 0, Var<PosReal>* intime = 0)	: sigma(insigma)	{

		up = inup;
		time = intime;

		if (up)	{
			Register(up);
		}
		if (time)	{
			Register(time);
		}
		Register(sigma);
		if (up) {
			SetName("normal process instant value\n");
		}
		else	{
			SetName("root normal process instant value\n");
		}
		Sample();
	}

	virtual double GetFiniteTimeTransitionLogProb()	{
		double tmp = val() - up->val();
		// double tmp = val() - up->val() + 0.5 * sigma->val() * time->val();
		return -0.5 * log(sigma->val() * time->val()) - 0.5 * tmp * tmp / time->val() / sigma->val();
	}

	// normal process has no stationary distribution!
	// let us define a default: standard normal distribution
	virtual double GetStationaryLogProb()	{
		return 0;
		// return - val() * val();
	}

	virtual void DrawFiniteTime()	{
		setval(Random::sNormal() * sqrt(time->val() * sigma->val()) + up->val());
	}

	virtual void DrawStationary()	{
		setval(Random::sNormal());
	}

	protected:

	Var<PosReal>* sigma;
};

class LogNormalProcessFiniteTimeAverage : public MarkovProcessBranchValue<Real,PosReal>	{

	public:


	LogNormalProcessFiniteTimeAverage(Var<Real>* inup, Var<Real>* indown, Var<PosReal>* intime, Var<PosReal>* innu)	{
		up = inup;
		down = indown;
		time = intime;
		nu = innu;

		Register(up);
		Register(down);
		Register(time);
		Register(nu);
		SetName("log normal process finite average");
		specialUpdate();
	}

	virtual void SetFiniteTimeAverage()	{
		value = (0.5 * (exp(up->val()) + exp(down->val()))) * time->val() * nu->val();
	}

	public :

	Var<PosReal>* nu;
};


class LogNormalTreeProcess : public RigidAutocorrelatedProcess<Real,PosReal>, public Additive	{

	public:

	LogNormalTreeProcess() {}

	LogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma, Var<PosReal>* innu)  {
		tree = intree;
		sigma = insigma;
		nu = innu;
		RecursiveCreate(GetRoot());
		GetNodeVal(GetRoot()->GetNode())->ClampAt(0);
	}

	~LogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	LengthTree* GetLengthTree()	{
		return tree;
	}

	Rvar<Real>* GetRootRate() {return GetNodeVal(GetRoot()->GetNode());}

	double GetExpVal(const Link* link)	{
		return exp(GetNodeVal(link->GetNode())->val());
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

	protected:

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

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = 0;
		if (! link->isRoot())	{
			time = GetLengthTree()->GetBranchLength(link->GetBranch());
		}

		// make the new instant value, and return the pointer
		return new NormalProcessInstantValue(sigma, vup, time);
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
		return new LogNormalProcessFiniteTimeAverage(vup, vdown, time, nu);
	}

	LengthTree* tree;
	Var<PosReal>* sigma;
	Var<PosReal>* nu;
};



#endif
