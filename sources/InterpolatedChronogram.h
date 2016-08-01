#ifndef INTERPOLATEDCHRONOGRAM_H
#define INTERPOLATEDCHRONOGRAM_H

#include "Chronogram.h"


class InterpolatedNodeDate : public Dvar<PosReal>	{

	public:

	InterpolatedNodeDate(Var<PosReal>* inRate, NodeDate* inup, NodeDate* indown, double inp)	{
		Register(inRate);
		up = inup;
		down = indown;
		p = inp;
	}

	void specialUpdate()	{
		setval(up->val() - (up->val() - down->val()) * p);
	}

	private:

	NodeDate* up;
	NodeDate* down;
	double p;

};


class InterpolatedChronogram : public Chronogram	{

	InterpolatedChronogram() {}

	InterpolatedChronogram(Tree* intree, Var<PosReal>* inRate)	{
		SetWithRoot(false);
		tree = intree;
		rate = inRate;
		// Create objects
		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());
		// Set values and make the tree ultra-metric
		double maxage = RecursiveSetNodesValues(GetRoot());
		RecursiveEqualizeLeafNodes(GetRoot(),maxage);
		RecursiveNormalizeTree(GetRoot(),maxage,true);
		RecursiveUpdateBranches(GetRoot());
	}

	~InterpolatedChronogram()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	Rvar<PosReal>* CreateNodeVal (const Link* link){
		if (link->GetDegree() == 2)	{
			int nup = 0;
			const Link* linkup = link->GetUp(nup);
			int ndown = 0;
			const Link* linkdown = link->Next()->GetUp(ndown);
			int N = nup + ndown;
			double p = ((double) nup) / ndown;
			return new InterpolatedNodeDate(rate,GetNodeDate(linkup->GetNode()),GetNodeDate(linkdown->GetNode()),p);
		}
		return new NodeDate(rate);
	}


};


#endif
