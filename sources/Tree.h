
#ifndef TREE_H
#define TREE_H

#include <cstdlib>
#include <string>
#include <iostream>
#include "TaxonSet.h"
#include "Random.h"

class Node {

	private:
	int index;
	string name;

	public:

	Node() : index(0) , name("") {}
	Node(string s) : index(0) , name(s) {}
	Node(const Node* from) : index(from->index) , name(from->name) {}

	virtual ~Node() {}

	virtual string GetName() const {return name;}
	virtual void SetName(string inname) {name = inname;}
	int GetIndex() const {return index;}
	void SetIndex(int i) {index = i;}
};

class Branch {

	private:
	int index;
	string name;

	public:

	Branch() : index(0) , name("") {}
	Branch(string s) : index(0) , name(s) {}
	Branch(const Branch* from) : index(from->index), name(from->name) {}

	virtual ~Branch() {}

	virtual string GetName() const {return name;}
	virtual void SetName(string inname) {name = inname;}
	int GetIndex() const {return index;}
	void SetIndex(int i) {index = i;}
};

class Link	{

	private:
	Link* next;
	Link* out;
	Branch* branch;
	Node* node;

	public:

	// double* tbl;

	Link()	{
		// tbl = 0;
		next = out = this;
		branch = 0;
		node = 0;
	}

	Link(const Link* from)	{
		// tbl = 0;
		next = out = this;
		node = 0;
		branch = 0;
	}

	Link* Next() const {return next;}
	Link* Out() const {return out;}
	Branch* GetBranch() const {return branch;}
	Node* GetNode() const {return node;}

	void SetBranch(Branch* inbranch)	{
		/*
		if (branch == inbranch)	{
			cout << "error in Link::SetBranch: branch and new branch are the same\n";
			exit(1);
		}
		delete branch;
		*/
		branch = inbranch;
	}
	void SetNode(Node* innode)	{
		/*
		if (node == innode)	{
			cout << "error in Link::SetNode: node and new node are the same\n";
			exit(1);
		}
		delete node;
		*/
		node = innode;
	}

	void SetOut(Link* inout)	{
		out = inout;
	}

	void SetNext(Link* innext)	{
		next = innext;
	}

	void AppendTo(Link* link)	{
		if (link)	{
			link->next = this;
		}
	}

	void Insert(Link* link)	{ // insert link after this
		link->next = next;
		next = link;
	}

	void InsertOut(Link* link)	{ // insert link as out
		link->out = this;
		out = link;
	}

	bool isLeaf() const	{
		return (next == this);
	}

	bool isUnary() const	{
		return (next->Next() == this && !isLeaf());
	}

	bool isRoot() const {
		return (out == this);
	}

	// degree : number of branches connecting to the node associated to this link
	int GetDegree()	const {
		int d = 1;
		const Link* link = next;
		while (link!=this)	{
			d++;
			link=link->next;
		}
		return d;
	}

	const Link* GetUp(int& d) const	{
		cerr << "in getup\n";
		exit(1);
		d = 1;
		const Link* link = Out();
		// const Link* link = link->Out();
		while (link->GetDegree() == 2)	{
			link = link->Next()->Out();
			d++;
		}
		return link;
		return 0;
	}
};


class NewickTree {

	public:

	virtual ~NewickTree() {}
	virtual const Link* GetRoot() const = 0;

	void ToStream(ostream& os) const;
	void ToStream(ostream& os, const Link* from) const;
	double ToStreamSimplified(ostream& os, const Link* from) const;

	const Link* GetLeftMostLink(const Link* from) const 	{
		if (from->isLeaf())	{
			return from;
		}
		return GetLeftMostLink(from->Next()->Out());
	}

	const Link* GetRightMostLink(const Link* from) const {
		if (from->isLeaf())	{
			return from;
		}
		const Link* link = from->Next();
		while (link->Next() != from)	{
			link = link->Next();
		}
		return GetRightMostLink(link->Out());
	}

	string GetLeftMost(const Link* from) const 	{
		if (from->isLeaf())	{
			return GetNodeName(from);
		}
		return GetLeftMost(from->Next()->Out());
	}

	string GetRightMost(const Link* from) const {
		if (from->isLeaf())	{
			return GetNodeName(from);
		}
		const Link* link = from->Next();
		while (link->Next() != from)	{
			link = link->Next();
		}
		return GetRightMost(link->Out());
	}

	static void Simplify()	{
		simplify = true;
	}

	void PrintTab(ostream& os)	{
		RecursivePrintTab(os,GetRoot());
	}

	void RecursivePrintTab(ostream& os, const Link* from)	{
		os << GetLeftMost(from) << '\t' << GetRightMost(from) << '\t' << GetNodeName(from) << '\n';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivePrintTab(os,link->Out());
		}
	}

	int GetNnode()	{
		return RecursiveGetNnode(GetRoot());
	}

	int GetNinternalNode()	{
		return RecursiveGetNinternalNode(GetRoot());
	}

	int RecursiveGetNinternalNode(const Link* from)	{
		int n = 0;
		if (! from->isLeaf())	{
			n++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveGetNinternalNode(link->Out());
		}
		return n;
	}

	int RecursiveGetNnode(const Link* from)	{
		int n = 1;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveGetNnode(link->Out());
		}
		return n;
	}

	protected:

	virtual string GetNodeName(const Link* link) const = 0;
	virtual string GetBranchName(const Link* link) const = 0;

	virtual string GetLeafNodeName(const Link* link) const {
		return GetNodeName(link);
	}

	static bool simplify;

};


class Tree : public NewickTree {

	public:

	Tree();
	// default constructor: set member pointers to 0

	Tree(const Tree* from);
	// clones the entire Link structure
	// but does NOT clone the Nodes and Branches
	// calls RecursiveClone

	Tree(string filename);
	// create a tree by reading into a file (netwick format expected)
	// calls ReadFromStream

	virtual ~Tree();
	// calls RecursiveDelete
	// but does NOT delete the Nodes and Branches

	// Delete the leaf pointing by the next link and set everithing right.
	void DeleteNextLeaf(Link* previous);

	// Delete the unary Node wich from is paart of and set everithing right.
	void DeleteUnaryNode(Link* from);

	Link* GetRoot() const {return root;}
	// const Link* GetRoot() const {return root;}
	const TaxonSet* GetTaxonSet() const {return taxset;}

	void RootAt(Link* newroot);

	void RegisterWith(const TaxonSet* taxset);
	// Registers all leaves of the tree with an external TaxonSet
	// the taxon set defines a map between taxon names and indices (between 0 and P-1)
	// the tree is recursively traversed
	// each leaf's name is looked for in the map of the taxon set
	// if not found : an error is emitted
	// otherwise, the leaf's index is set equal to the index of the corresponding taxon

	bool RegisterWith(const TaxonSet* taxset, Link* from, int& tot);
	// recursive function called by RegisterWith

	double GetBranchLength(const Link* link) const	{
		return atof(GetBranchName(link).c_str());
	}

	double GetMaxHeight(const Link* from) const	{

		double max = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = GetMaxHeight(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		if (! from->isRoot())	{
			max += GetBranchLength(from);
		}
		return max;
	}

	double GetMinHeight(const Link* from) const	{

		double min = -1;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = GetMinHeight(link->Out());
			if ((min == -1) || (min > tmp))	{
				min = tmp;
			}
		}
		if (! from->isRoot())	{
			min += GetBranchLength(from);
		}
		return min;
	}

	void ToStreamRenorm(const Link* from, ostream& os, double normfactor)	{
		
		if (from->isLeaf())	{
			os << GetNodeName(from);
		}
		else	{
			os << "(";
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				ToStreamRenorm(link->Out(),os,normfactor);
				if (link->Next() != from)	{
					os << ",";
				}
			}
			os << ")";
			os << GetNodeName(from);
		}
		if (from->isRoot())	{
			os << ";\n";
		}
		else	{
			os << ":";
			os << GetBranchLength(from) * normfactor;
		}
	}

	string GetBranchName(const Link* link) const {return link->GetBranch()->GetName();}

	string GetNodeName(const Link* link) const {
		return link->GetNode()->GetName();
		/*
		if (! link->isLeaf())	{
			return link->GetNode()->GetName();
		}
		string s = link->GetNode()->GetName();
		unsigned int l = s.length();
		unsigned int i = 0;
		while ((i < l) && (s[i] != '_')) i++;
		if (i == l)	{
			cerr << "error in get name\n";
			exit(1);
		}
		i++;
		return s.substr(i,l-i);
		*/
	}
	// trivial accessors
	// they can be useful to override, so as to bypass Branch::GetName() and Node::GetName()

	void EraseInternalNodeName();
	void EraseInternalNodeName(Link* from);

	// void Print(ostream& os,const Link* from) const ;
	// void Print(ostream& os) const;
	// printing int netwick format

	/*
	void RecursiveCreateTBL(Link* from, int Nstate);
	void RecursiveDeleteTBL(Link* from);
	*/

	// for creating and deleting the arrays of conditional likelihoods associated with each link
	// called by PhyloProcess

	unsigned int GetSize()	const {
		return GetSize(GetRoot());
	}

	int GetSize(const Link* from)	const {
		if (from->isLeaf())	{
			return 1;
		}
		else	{
			int total = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				total += GetSize(link->Out());
			}
			return total;
		}
		return 0;
	}

	int GetFullSize(const Link* from)	const {
		if (from->isLeaf())	{
			return 1;
		}
		else	{
			int total = 1;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				total += GetFullSize(link->Out());
			}
			return total;
		}
		return 0;
	}

	virtual const Link* GetLCA(string tax1, string tax2)	{
		bool found1 = false;
		bool found2 = false;
		const Link* link= RecursiveGetLCA(GetRoot(),tax1,tax2,found1,found2);
		// Print(cerr);
		// cerr << tax1 << '\t' << tax2 << '\n';
		// Print(cerr,link);
		// cerr << '\n' << '\n';
		/*
		cerr << link << '\t' << tax1 << '\t' << tax2 << '\t';
		if (link)	{
			cerr << GetLeftMost(link) << '\t' << GetRightMost(link) << '\n';
		}
		else	{
			cerr << '\n';
		}
		*/
		return link;
	}

	virtual const Link* GetLCA(const Link* from1, const Link* from2)	{
		bool found1 = false;
		bool found2 = false;
		const Link* link= RecursiveGetLCA(GetRoot(),from1,from2,found1,found2);
		return link;
	}

	void Subdivide(Link* from, int Ninterpol);

	string Reduce(Link* from = 0)	{

		if (! from)	{
			from = GetRoot();
		}
		if (from->isLeaf())	{
			cerr << from->GetNode()->GetName() << '\n';;
			return from->GetNode()->GetName();
		}
		else	{
			string name = "None";
			for (Link* link = from->Next(); link!=from; link=link->Next())	{
				string tmp = Reduce(link->Out());
				if (tmp == "diff")	{
					name = "diff";
				}
				else if (name == "None")	{
					name = tmp;
				}
				else if (name != tmp)	{
					name = "diff";
				}
			}
			cerr << '\t' << name << '\n';
			from->GetNode()->SetName(name);
			return name;
		}
		return "";
	}

	void PrintReduced(ostream& os, const Link* from = 0)	{
		if (!from)	{
			from = GetRoot();
		}
		if (from->GetNode()->GetName() != "diff")	{
			os << from->GetNode()->GetName();
		}
		else	{
			os << '(';
			for (const Link* link = from->Next(); link!=from; link=link->Next())	{
				PrintReduced(os,link->Out());
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
		if (from->isRoot())	{
			os << ";\n";
		}
	}

	const Link* ChooseInternalNode()	{
		int n = CountInternalNodes(GetRoot());
		int m = (int) (n * Random::Uniform());
		const Link* tmp;
		const Link* chosen = ChooseInternalNode(GetRoot(),tmp,m);
		if (! chosen)	{
			cerr << "error in choose internal node: null pointer\n";
			exit(1);
		}
		return chosen;
	}

	int CountInternalNodes(const Link* from);
	const Link* ChooseInternalNode(const Link* from, const Link*& chosenup, int& n);
	int CountNodes(const Link* from);
	const Link* ChooseNode(const Link* from, const Link*& chosenup, int& n);

	protected:

	// returns 0 if not found
	// returns link if found (then found1 and found2 must
	const Link* RecursiveGetLCA(const Link* from, string tax1, string tax2, bool& found1, bool& found2)	{
		const Link* ret= 0;
		if (from->isLeaf())	{
			// found1 |= (from->GetNode()->GetName() == tax1);
			// found2 |= (from->GetNode()->GetName() == tax2);
			string name1 = GetLeafNodeName(from).substr(0,tax1.size());
			string name2 = GetLeafNodeName(from).substr(0,tax2.size());
			found1 |= (name1 == tax1);
			found2 |= (name2 == tax2);
			/*
			found1 |= (GetLeafNodeName(from) == tax1);
			found2 |= (GetLeafNodeName(from) == tax2);
			*/
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		else	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				bool tmp1 = false;
				bool tmp2 = false;
				const Link* ret2 = RecursiveGetLCA(link->Out(),tax1,tax2,tmp1,tmp2);
				found1 |= tmp1;
				found2 |= tmp2;
				if (ret2)	{
					if (ret)	{
						cerr << "error : found node twice\n";
						cerr << tax1 << '\t' << tax2 << '\n';
						ToStream(cerr,ret2->Out());
						cerr << '\n';
						ToStream(cerr,ret->Out());
						cerr << '\n';
						exit(1);
					}
					ret = ret2;
				}
			}
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		return ret;
	}

	const Link* RecursiveGetLCA(const Link* from, const Link* from1, const Link* from2, bool& found1, bool& found2)	{
		const Link* ret= 0;
		found1 |= (from == from1);
		found2 |= (from == from2);
		if (! ret)	{
			if (found1 && found2)	{
				ret = from;
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			bool tmp1 = false;
			bool tmp2 = false;
			const Link* ret2 = RecursiveGetLCA(link->Out(),from1,from2,tmp1,tmp2);
			found1 |= tmp1;
			found2 |= tmp2;
			if (ret2)	{
				if (ret)	{
					cerr << "error : found node twice\n";
					cerr << from1 << '\t' << from2 << '\n';
					ToStream(cerr,ret2->Out());
					cerr << '\n';
					ToStream(cerr,ret->Out());
					cerr << '\n';
					exit(1);
				}
				ret = ret2;
			}
		}
		if (! ret)	{
			if (found1 && found2)	{
				ret = from;
			}
		}
		return ret;
	}

	void ReadFromStream(istream& is);
	// reading a tree from a stream:
	// recursively invokes the two following functions

	Link* ParseGroup(string input, Link* from);
	// a group is an expression of one of the two following forms:
	//
	// 	(Body)Node_name
	// 	(Body)Node_name:Branch_name
	//
	// where Body is a list, and Node_name and Branch_name are 2 strings
	// Node_name may be an empty string
	//
	// Node_name cannot contain the ':' character, but Branch_name can
	// thus, if the group reads "(BODY)A:B:C"
	// then Node_name = "A" and Branch_name = "B:C"

	Link* ParseList(string input, Node* node);
	// a list is an expression of the form X1,X2,...Xn
	// where Xi is a group

	void RecursiveClone(const Link* from, Link* to);
	// Used by Tree(const Tree* from)
	// only clone the Links, and their mutual relations
	// does not copy the Node or Branch objects

	// deletes the link structure
	// does not delete the Node or Branch objects
	void RecursiveDelete(Link* from);

	void SetRoot(Link* link) {root = link;}

	// data fields
	// just 2 pointers, to the root and to a list of taxa
	Link* root;
	const TaxonSet* taxset;
};



#endif // TREE_H
