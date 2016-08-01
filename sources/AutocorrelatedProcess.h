#ifndef AUTOCORRELATEDPROCESS_H
#define AUTOCORRELATEDPROCESS_H

#include <sstream>

#include "BaseType.h"
#include "Var.h"
#include "ValTree.h"


// a Markov process taking values of type T
// abstract class. interface specifying the general form transition / stationary probabilities
// and implementing Sample() and logProb()
//
// this is a Rvar<T>: a random node of the model's graph

template<class T> class MarkovProcessInstantValue : public virtual Rvar<T>	{

	public:

	// these 4 abstract methods will be model-specific
	virtual double GetFiniteTimeTransitionLogProb() = 0;
	virtual double GetStationaryLogProb() = 0;
	virtual void DrawFiniteTime() = 0;
	virtual void DrawStationary() = 0;

	protected:
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

	Var<T>* up;
	Var<PosReal>* time;
};

// the integral of a Markov process of type T along a given time interval
// and conditional on the values at both ends
//
// this is a Rvar<T>: the integral is random

template<class T> class MarkovProcessFiniteTimeIntegral : public virtual Rvar<T>	{

	// these 2 abstract methods will be model-specific
	virtual double GetFiniteTimeIntegralLogProb() = 0;
	virtual void  DrawFiniteIntegral() = 0;

	protected:
	virtual void drawSample()	{
		if (up)	{ // finite time transition
			this->DrawFiniteTimeIntegral();
		}
		else	{
			this->setval(0);
		}
	}

	virtual double logProb()	{
		if (up)	{ // finite time transition prob
			return GetFiniteTimeIntegralLogProb();
		}
		if (this->val() != 0)	{
			return Random::INFPROB;
		}
		return 0;
	}

	Var<T>* up;
	Var<T>* down;
	Var<PosReal>* time;
};


template<class S, class T> class MarkovProcessBranchValue : public virtual Dvar<T>	{

	// this abstract method will be model-specific
	virtual void SetFiniteTimeAverage() = 0;

	protected:
	void specialUpdate()	{
		SetFiniteTimeAverage();
	}

	Var<S>* up;
	Var<S>* down;
	Var<PosReal>* time;
};

// the expectation of a Markov process of type T along a given time interval
// and conditional on the values at both ends
//
// this is a Dvar<T>: deterministic node in the graph

template<class T> class MarkovProcessFiniteTimeAverage : public MarkovProcessBranchValue<T,T>	{};


// defines an autocorrelated Markov process on a valuated tree (LengthTree)
// taking instant values of type V
//
// defines values (Rvar<V>) at each node (these will be the instant values of the Markov process)
// and on each branch (these will be the integrals along the time)
//
// implements all global Monte Carlo sampling methods
// using recursive calls to the corresponding methods for all node- and branch- associated values
//


template<class U, class V> class FlexibleAutocorrelatedProcess : public MCMC , public NodeBranchValPtrTree< Rvar<U>, Rvar<V> > {

	public:

	virtual LengthTree* GetLengthTree() = 0;

	virtual Tree* GetTree() {
		return GetLengthTree()->GetTree();
	}

	virtual double Move(double tuning)	{
		int n = 0;
		double tot = Move(this->GetRoot(),tuning,n);
		return tot / n;
	}

	virtual double GetLogProb()	{
		return GetLogProb(this->GetRoot());
	}

	protected:

	virtual void drawSample()	{
		SampleNode(this->GetRoot());
		SampleBranch(this->GetRoot());
	}

	virtual void SampleNode(Link* from)	{
		this->GetNodeVal(from->GetNode())->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			SampleNode(link->Out());
		}
	}

	virtual void SampleBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			this->GetBranchVal(link->GetBranch())->Sample();
			SampleBranch(link->Out());
		}
	}

	virtual double Move(Link* from, double tuning, int& count)	{
		double total = this->GetNodeVal(from->GetNode())->Move(tuning);
		count++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->Move(tuning);
			count++;
			total += Move(link->Out(),tuning,count);
		}
		return total;
	}

	virtual double GetLogProb(Link* from)	{
		double total = this->GetNodeVal(from->GetNode())->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->GetLogProb();
			total += GetLogProb(link->Out());
		}
		return total;
	}
};

template<class U, class V> class RigidAutocorrelatedProcess : public MCMC , public NodeBranchValPtrTree< Rvar<U>, Dvar<V> > {

	public:

	virtual LengthTree* GetLengthTree() = 0;

	virtual Tree* GetTree() {
		return GetLengthTree()->GetTree();
	}

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

	protected:

	virtual void drawSample()	{
		SampleNode(this->GetRoot());
	}

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

	/*
	virtual double Move(Link* from, double tuning, int& count)	{
		double total = this->GetNodeVal(from->GetNode())->Move(tuning);
		count++;
		int nlink = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next()){
			nlink++;
		}
		Link** linkarray = new Link*[nlink];
		int k = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next()){
			linkarray[k] = link;
			k++;
		}
		for (int i=0; i<nlink; i++)	{
			total += Move(linkarray[i]->Out(),tuning,count);
		}
		for (int i=nlink-1; i>=0; i--)	{
			total += Move(linkarray[i]->Out(),tuning,count);
		}
		delete[] linkarray;
		total += this->GetNodeVal(from->GetNode())->Move(tuning);
		count++;
		return total;
	}
	*/

	virtual double GetLogProb(Link* from)	{
		double total = this->GetNodeVal(from->GetNode())->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetLogProb(link->Out());
		}
		return total;
	}
};

#endif

