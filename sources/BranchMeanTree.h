
#ifndef BRANCHMEANTREE_H
#define BRANCHMEANTREE_H

template<class U> class BranchMeanTree  {

	public:

	BranchMeanTree(Tree* intree, int indim) {
		tree = intree;
		dim = indim;
		Reset();
	}

	Tree* GetTree() {return tree;}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(BranchVarTree<U>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->Print(os);
	}

	private:

	void RecursiveAdd(BranchVarTree<U>* sample, LengthTree* chronogram, Link* from)	{
		/*
		if (from->isRoot() && (WithRoot()))	{
			GetBranchVal(from->GetBranch())->Add(sample->GetBranchVal(from->GetBranch()));
		}
		*/
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out());
			branchval[link->GetBranch()].Add(*sample->GetBranchVal(link->GetBranch()));
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveReset(Link* from)	{
		/*
		if (from->isRoot() && (WithRoot()))	{
			branchval[from->GetBranch()] = U(dim);
			branchval[from->GetBranch()].SetAtZero();
		}
		*/
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			branchval[link->GetBranch()] = U(dim);
			branchval[link->GetBranch()].SetAtZero();
			meantime[link->GetBranch()] = 0 ;
			vartime[link->GetBranch()] = 0;
		}
	}

	void RecursiveNormalise(Link* from)	{
		/*
		if (from->isRoot() && (WithRoot()))	{
			branchval[from->GetBranch()].ScalarMultiplication(1.0 / size);
		}
		*/
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			branchval[link->GetBranch()].ScalarMultiplication(1.0 / size);
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			ostringstream s;
			if (link->Out()->isLeaf())	{
				string t = link->Out()->GetNode()->GetName();
				unsigned int i = 0;
				char c = ' ';
				while ((i<t.size()) && (c != '@'))	{
					c = t[i];
					if (c != '@')	{
						s << c;
					}
					i++;
				}
				s << '@';
			}
			s << Translate(branchval[link->GetBranch()]);
			link->Out()->GetNode()->SetName(s.str());
			RecursiveSetName(os, link->Out());
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	string Translate(U from)	{
		ostringstream s;
		for (int i=0; i<dim; i++)	{
			s << from[i];
			if (i<dim-1)	{
				s << '&';
			}
		}
		return s.str();
	}

	Tree* tree;
	int size;
	int dim;


	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;
	map<const Branch*, U> branchval;

};

#endif

