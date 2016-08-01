
#ifndef NODEPROCESS_H
#define NODEPROCESS_H

#include <sstream>

#include "ValTree.h"
#include "RandomTypes.h"

// iid random variables of type Rvar<V>
// associated to each branch of the tree

template<class V> class NodeProcess : public MCMC , public virtual NodeValPtrTree< Rvar<V> > {

	public:

	NodeProcess(Tree* intree) : tree(intree) {
	}

	Tree* GetTree() {return tree;}

	V val(const Node* node)	{
		return this->GetNodeVal(node)->val();
	}

	void setval(const Node* node, V inval)	{
		this->GetNodeVal(node)->setval(inval);
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

	protected:

	virtual void drawSample(const Link* from)	{
		this->GetNodeVal(from->GetNode())->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			this->GetNodeVal(link->GetNode())->Sample();
			drawSample(link->Out());
		}
	}


	virtual double Move(const Link* from, double tuning, int& count)	{
		double total = 0;
		total += this->GetNodeVal(from->GetNode())->Move(tuning);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetNodeVal(link->GetNode())->Move(tuning);
			count++;
			total += Move(link->Out(),tuning,count);
		}
		return total;
	}

	virtual double GetLogProb(const Link* from)	{
		double total = 0;
		total += this->GetNodeVal(from->GetNode())->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetNodeVal(link->GetNode())->GetLogProb();
			total += GetLogProb(link->Out());
		}
		return total;
	}

	Tree* tree;
};

#endif // NODEPROCESS_H
