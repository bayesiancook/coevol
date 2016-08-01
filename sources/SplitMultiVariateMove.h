

#ifndef SPLITMVMOVE_H
#define SPLITMVMOVE_H

#include "SplitTree.h"
#include "MultiVariateTreeProcess.h"
#include "Move.h"

class SplitMultiVariateBranchMove : public BranchValPtrTree<Mnode>, public MCUpdate {


	public:

	SplitMultiVariateBranchMove(Tree* intree, MultiVariateTreeProcess* inprocess, double intuning)	{
		SetWithRoot(false);
		process = inprocess;
		tree = intree;
		tuning = intuning;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree()	{
		return tree;
	}

	int GetDim()	{
		return process->GetDim();
	}

	double Move(double tuning_modulator = 1.0)	{
		int n = 0;
		double tot =  RecursiveMove(GetRoot(),tuning * tuning_modulator ,n);
		return tot / n;
	}

	private:

	SplitTree* GetSplitTree()	{
		SplitTree* tmp = dynamic_cast<SplitTree*>(process->GetTree());
		if (!tmp)	{
			cerr << "error in split mv move: null split tree\n";
			exit(1);
		}
		return tmp;
	}

	Mnode* CreateBranchVal(const Link* from)	{

		const Link* fromimage = GetSplitTree()->GetImage(from);

		if (!fromimage->GetBranch())	{
			cerr << "error in create: null branch\n";
			exit(1);
		}
		Mnode* mnode = new Mnode;
	 	process->GetMultiNormal(fromimage->GetNode())->Register(mnode);
		const Link* extend = fromimage->Out();
		int k = 0;
		while (extend->GetDegree() == 2)	{
			process->GetMultiNormal(extend->GetNode())->Register(mnode);
			extend = extend->Next()->Out();
			k++;
		}
	 	process->GetMultiNormal(extend->GetNode())->Register(mnode);
		return mnode;
	}

	double RecursiveMove(const Link* from, double tuning, int& n)	{
		double tot = 0;
		if (! from->isRoot())	{
			tot += Move(from,tuning);
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveMove(link->Out(),tuning,n);
		}
		return tot;
	}


	double Move(const Link* from, double tuning)	{

		GetBranchVal(from->GetBranch())->Corrupt(true);

		const Link* fromimage = GetSplitTree()->GetImage(from);
		if (!fromimage->GetBranch())	{
			cerr << "error in move: null branch\n";
			exit(1);
		}
		double* delta = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			delta[i] = tuning * (Random::Uniform() - 0.5);
		}

		double alpha =  Random::sGamma(1);
		double beta = Random::sGamma(1);

		const Link* extend = fromimage->Out();
		int n = GetSplitTree()->GetSplitBranch(extend->GetBranch())->GetSplitOrder();
		if (!extend->GetBranch())	{
			cerr << "error in move: null branch\n";
			exit(1);
		}
		int k = 1;

		while (extend->GetDegree() == 2)	{
			double f = ((double) k / n);
			double x = exp(alpha * log(f) + beta * log(1-f));
			process->GetMultiNormal(extend->GetNode())->Shift(delta,x);
			k++;
			extend = extend->Next()->Out();
		}

		delete[] delta;

		double logRatio = GetBranchVal(from->GetBranch())->Update();
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			GetBranchVal(from->GetBranch())->Corrupt(false);
			GetBranchVal(from->GetBranch())->Restore();
		}
		return (double) accepted;
	}

	MultiVariateTreeProcess* process;
	Tree* tree;
	double tuning;

};

class SplitMultiVariateNodeMove : public NodeValPtrTree<Mnode>, public MCUpdate	{

	public:

	SplitMultiVariateNodeMove(MultiVariateTreeProcess* inprocess, double intuning)	{
		process = inprocess;
		tuning = intuning;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree()	{
		return process->GetTree();
	}

	int GetDim()	{
		return process->GetDim();
	}

	double Move(double tuning_modulator = 1.0)	{
		int n = 0;
		double tot =  RecursiveMove(GetRoot(),tuning * tuning_modulator ,n);
		return tot / n;
	}

	private:

	SplitTree* GetSplitTree()	{
		SplitTree* tmp = dynamic_cast<SplitTree*>(process->GetTree());
		if (!tmp)	{
			cerr << "error in split mv move: null split tree\n";
			exit(1);
		}
		return tmp;
	}

	Mnode* CreateNodeVal(const Link* from)	{
		if (from->GetDegree() == 2)	{
			return 0;
		}
		Mnode* mnode = new Mnode;
		process->GetMultiNormal(from->GetNode())->Register(mnode);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			const Link* extend = link->Out();
			while (extend->GetDegree() == 2)	{
				process->GetMultiNormal(extend->GetNode())->Register(mnode);
				extend = extend->Next()->Out();
			}
			process->GetMultiNormal(extend->GetNode())->Register(mnode);
		}
		return mnode;
	}

	double RecursiveMove(const Link* from, double tuning, int& n)	{
		if (from->GetDegree() == 2)	{
			return RecursiveMove(from->Next()->Out(),tuning,n);
		}
		else	{
			double tot = 0;
			tot += Move(from,tuning);
			n++;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				tot += RecursiveMove(link->Out(),tuning,n);
			}
			return tot;
		}
	}


	double Move(const Link* from, double tuning)	{

		GetNodeVal(from->GetNode())->Corrupt(true);
		double* delta = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			delta[i] = tuning * (Random::Uniform() - 0.5);
		}

		process->GetMultiNormal(from->GetNode())->Shift(delta);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			int n = GetSplitTree()->GetSplitBranch(link->GetBranch())->GetSplitOrder();
			int k = n-1;
			const Link* extend = link->Out();
			while (extend->GetDegree() == 2)	{
				process->GetMultiNormal(extend->GetNode())->Shift(delta,((double) k) / n);
				k--;
				extend = extend->Next()->Out();
			}
		}

		delete[] delta;

		double logRatio = GetNodeVal(from->GetNode())->Update();
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			GetNodeVal(from->GetNode())->Corrupt(false);
			GetNodeVal(from->GetNode())->Restore();
		}
		return (double) accepted;
	}

	MultiVariateTreeProcess* process;
	double tuning;

};

#endif

