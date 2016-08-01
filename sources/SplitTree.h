
#ifndef SPLITTREE_H
#define SPLITTREE_H

#include "Tree.h"

class SplitBranch : public Branch	{

	public:

	SplitBranch(const Branch* inmother, int inn) : Branch(inmother), mother(inmother) , n(inn) {}
	virtual ~SplitBranch() {}

	const Branch* GetMother() {return mother;}
	double GetSplitFraction() {return 1.0 / n;}
	int GetSplitOrder() {return n;}

	private:

	const Branch* mother;
	int n;

};


class SplitNode : public Node	{

	public:

	SplitNode(const Node* inup, const Node* indown, int ink, int inn) : Node(indown), up(inup), down(indown), k(ink), n(inn) {
		if (!inup || !indown)	{
			cerr << "error : null node\n";
			exit(1);
		}
	}
	virtual ~SplitNode() {}

	const Node* GetDown() {return down;}
	const Node* GetUp() {return up;}
	double GetSplitFraction() {return ((double) k) / n;}

	bool Interpolates() {return (up != down);}

	private:

	const Node* up;
	const Node* down;
	int k;
	int n;
};

class SplitTree : public Tree	{

	public:

	SplitTree(const Tree* infrom, int insplitn);

	SplitBranch* GetSplitBranch(Branch* branch)	{
		SplitBranch* splitbranch = dynamic_cast<SplitBranch*>(branch);
		if (!splitbranch)	{
			cerr << "error in SplitBranch::GetSplitBranch: null pointer\n";
			cerr << branch << '\t' << splitbranch << '\n';
			exit(1);
		}
		return splitbranch;
	}

	SplitNode* GetSplitNode(Node* node)	{
		SplitNode* tmp = dynamic_cast<SplitNode*>(node);
		if (! tmp)	{
			cerr << "error in GetSplitNode\n";
			exit(1);
		}
		return tmp;
	}
	// image of a node of the initial tree in the final split tree
	const Link* GetImage(const Link* link)	{
		if (linkmap.find(link) == linkmap.end())	{
			cerr << "error in SplitTree::GetImage\n";
			exit(1);
		}
		return linkmap[link];
	}

	private:

	const Tree* from;
	int splitn;
	void RecursiveCloneAndSubdivide(const Link* from, Link* to);
	void MakeNewBranch(const Branch* from, Link* link1, Link* link2, const Node* node1, const Node* node2, int n);

	map<const Link*, const Link*> linkmap;

};

SplitTree::SplitTree(const Tree* infrom, int insplitn) : Tree() , from(infrom), splitn(insplitn)  {
	taxset = from->GetTaxonSet();
	root = new Link(from->GetRoot());
	linkmap[from->GetRoot()] = root;
	root->InsertOut(root);
	RecursiveCloneAndSubdivide(from->GetRoot(),root);
}

void SplitTree::MakeNewBranch(const Branch* from, Link* link1, Link* link2, const Node* node1, const Node* node2, int n)	{

	int k = 0;
	Link* current = link1;
	while (k<n)	{
		SplitBranch* branch = new SplitBranch(from, n);
		current->SetBranch(branch);
		if (k != n-1)	{
			Node* node = new SplitNode(node1,node2,(n-k-1),n);
			Link* l1 = new Link();
			l1->SetBranch(branch);
			l1->SetNode(node);
			current->SetOut(l1);
			l1->SetOut(current);
			Link* l2 = new Link();
			l2->SetNode(node);
			l1->SetNext(l2);
			l2->SetNext(l1);
			current = l2;
		}
		else	{
			link2->SetBranch(branch);
			link2->SetOut(current);
			current->SetOut(link2);
		}
		k++;
	}
}

void SplitTree::RecursiveCloneAndSubdivide(const Link* from, Link* to)	{

	Node* node = new SplitNode(from->GetNode(),from->GetNode(),1,1);
	to->SetNode(node);
	const Link* linkfrom = from->Next();
	Link* linkto = to;
	while (linkfrom != from)	{
		Link* newnext = new Link(linkfrom); // newnext points to same node and branch as linkfrom
		linkmap[linkfrom] = newnext;
		newnext->SetNode(node);
		linkto->Insert(newnext);
		Link* newout = new Link(linkfrom->Out()); // idem, same node and branch as linkfrom->Out()
		newout->SetNode(linkfrom->Out()->GetNode());
		linkmap[linkfrom->Out()] = newout;
		int n = splitn;
		if (!n)	{
			cerr << "error in split tree: cannot split a branch into " << n << "component\n";
			exit(1);
		}
		if (n < 0)	{
			int tmp = ((int) (-n * atof(linkfrom->GetBranch()->GetName().c_str()))) + 1;
			n = tmp;
		}

		MakeNewBranch(linkfrom->GetBranch(),newnext,newout,linkfrom->Out()->GetNode(),linkfrom->GetNode(),n);
		RecursiveCloneAndSubdivide(linkfrom->Out(),newout);
		linkfrom = linkfrom->Next();
		linkto = linkto->Next();
	}
}

#endif
