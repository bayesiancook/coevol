#ifndef NORMALTREEPROCESS_H
#define NORMALTREEPROCESS_H

#include "AutocorrelatedProcess.h"

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

		Register(sigma);
		if (up)	{
			Register(up);
		}
		if (time)	{
			Register(time);
		}
		Sample();
	}

	virtual double GetFiniteTimeTransitionLogProb()	{
		double tmp = val() - up->val();
		return -log(sigma->val()) - 0.5 * tmp * tmp / sigma->val() / sigma->val() / time->val();
	}

	// normal process has no stationary distribution!
	// let us define a default: standard normal distribution
	virtual double GetStationaryLogProb()	{
		return 0;
		// return - val() * val();
	}

	virtual void DrawFiniteTime()	{
		setval(Random::sNormal() * sigma->val() * sqrt(time->val()) + up->val());
	}

	virtual void DrawStationary()	{
		setval(Random::sNormal());
	}

	protected:

	Var<PosReal>* sigma;
};


// implements an approximation of the cumulated process along a given period of time
// simply the average of the two extreme values
// multiplied by the time
//
class NormalProcessFiniteTimeAverage : public MarkovProcessFiniteTimeAverage<Real>	{

	public:

	NormalProcessFiniteTimeAverage(Var<Real>* inup, Var<Real>* indown, Var<PosReal>* intime)	{
		up = inup;
		down = indown;
		time = intime;
		Register(up);
		Register(down);
		Register(time);
		specialUpdate();
	}

	virtual void SetFiniteTimeAverage()	{
		// SETVAL problem
		value = (0.5 * (up->val() + down->val())) * time->val();
	}
};

// implements an approximation of the exponential of the cumulated process along a given period of time
// simply the average of the exponential of the two extreme values
// multiplied by the time
//
// note that this is a MarkovProcessBranchValue<Real,PosReal>
// because the instant values of the process are real
// but the result of the exponential is a positive real
//
class LogNormalProcessFiniteTimeAverage : public MarkovProcessBranchValue<Real,PosReal>	{

	public:

	LogNormalProcessFiniteTimeAverage(Var<Real>* inup, Var<Real>* indown, Var<PosReal>* intime)	{
		up = inup;
		down = indown;
		time = intime;

		Register(up);
		Register(down);
		Register(time);
		specialUpdate();
	}

	virtual void SetFiniteTimeAverage()	{
		// SETVAL problem
		value = (0.5 * (exp(up->val()) + exp(down->val()))) * time->val();
	}
};


class LogNormalTreeProcess : public RigidAutocorrelatedProcess<Real,PosReal>	{

	public:

	LogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma)  {
		tree = intree;
		sigma = insigma;
		RecursiveCreate(GetRoot());
		GetNodeVal(GetRoot()->GetNode())->ClampAt(0);
	}

	~LogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	LengthTree* GetLengthTree()	{
		return tree;
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

	private:

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

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new NormalProcessInstantValue(sigma, vup, time);
	}

	virtual Dvar<PosReal>* CreateBranchVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the instant value associated to this node
		Var<Real>* vdown = GetNodeVal(link->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new LogNormalProcessFiniteTimeAverage(vup, vdown, time);
	}

	LengthTree* tree;
	Var<PosReal>* sigma;
};


class NormalTreeProcess : public RigidAutocorrelatedProcess<Real,Real>	{

	public:

	NormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma)  {
		tree = intree;
		sigma = insigma;
		RecursiveCreate(GetRoot());
	}

	~NormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	LengthTree* GetLengthTree()	{
		return tree;
	}

	private:

	// create the instant value that will be associated to the node accessible via the link (link->GetNode())
	// this function will be called by RecursiveCreate, such as defined in NodeValPtrTree
	//
	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new NormalProcessInstantValue(sigma, vup, time);
	}

	// create the branch length (= time * meanrate) that will be associated to the branch accessible via the link (link->GetBranch())
	// this function will be called by RecursiveCreate, such as defined in BranchValPtrTree
	//
	virtual Dvar<Real>* CreateBranchVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the instant value associated to this node
		Var<Real>* vdown = GetNodeVal(link->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new NormalProcessFiniteTimeAverage(vup, vdown, time);
	}

	LengthTree* tree;
	Var<PosReal>* sigma;
};



#endif
