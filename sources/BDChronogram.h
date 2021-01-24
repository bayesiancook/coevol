
#ifndef BDCHRONO
#define BDCHRONO

#include "Chronogram.h"
#include <algorithm>
#include <utility>
#include <list>


class BDNodeDate : public NodeDate	{

	protected:

	double min, max;

	public:

	BDNodeDate(Var<PosReal>* inRate, Var<PosReal>* indelta, Var<PosReal>* inkappa, bool inleaf, bool inroot){
		delta = indelta;
		kappa = inkappa;
		leaf = inleaf;
		root = inroot;
		Register(inRate);
		Register(delta);
		Register(kappa);
		SetName("date");
	}

	virtual double logProb(){
		if (leaf)	{
			if (val())	{
				cerr << "error : leaf node : " << val() << '\n';
				exit(1);
			}
			return 0;
		}
		if (root)	{
			if (val() != 1)	{
				cerr << "error : root node : " << val () << '\n';
				exit(1);
			}
			return 0;
		}
		double t1 = 1.0;
		double expo = exp(-delta->val() * val());
		double expo1 = exp(-delta->val() * t1);
		double nu = 1 - delta->val() * expo1 / (kappa->val() + (delta->val() - kappa->val()) * expo1);
		double Prho = delta->val() / (kappa->val() + (delta->val() - kappa->val()) * expo);
		double g = kappa->val() * Prho * Prho * expo / nu;
		double h = log(g);
		if (std::isnan(h))	{
			cerr << "error in bd node date : nan\n";
			exit(1);
		}
		if (std::isinf(h))	{
			cerr << "error in bd node data : inf\n";
			exit(1);
		}
		return h;
	}

	protected:

	Var<PosReal>* delta;
	Var<PosReal>* kappa;
	bool leaf;
	bool root;

};

class BDChronogram : public Chronogram	{

	public:

	BDChronogram(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* indelta, Var<PosReal>* inkappa)	{

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		delta = indelta;
		kappa = inkappa;

		// Create objects
		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		// Set values and make the tree ultra-metric
		double maxage = RecursiveSetNodesValues(GetRoot());
		RecursiveEqualizeLeafNodes(GetRoot(),maxage);
		RecursiveNormalizeTree(GetRoot(),maxage,true);
		RecursiveUpdateBranches(GetRoot());
	}

	Rvar<PosReal>* CreateNodeVal (const Link* link){
		return new BDNodeDate(rate,delta,kappa,link->isLeaf(), link->isRoot());
	}

	Var<PosReal>* delta;
	Var<PosReal>* kappa;
};

#endif

