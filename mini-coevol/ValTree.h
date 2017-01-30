
#ifndef VALTREE_H
#define VALTREE_H

#include <sstream>

#include "BaseType.h"
#include "Var.h"
#include "AbstractTree.h"

// implements recursive allocations and basic accessors

enum BranchValType {MEAN,INTEGRAL};

template<class V> class BranchValPtrTree;

template<class U> class NodeValPtTree;

template<class U, class V> class NodeBranchValPtrTree; // : public NodeValPtrTree<U>, public BranchValPtrTree<V>;

template<class V> class _BranchValPtrTree : public virtual AbstractTree {

	public:

	bool WithRoot()	{
		return withRoot;
	}

	void SetWithRoot(bool in)	{
		withRoot = in;
	}

	V* GetBranchVal(const Branch* branch)	{
		return branchval[branch];
	}

	void SetBranchVal(const Branch* branch, V* in)	{
		if ((! WithRoot()) && (! branch))	{
			cerr << "error in _BranchValPtrTree: null branch pointer\n";
			throw;
		}
		branchval[branch] = in;
	}

	virtual void RecursiveCreate(const Link* from)	{
		if (WithRoot() && from->isRoot())	{
			V* v= CreateBranchVal(from);
			branchval[0] = v;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			V* v= CreateBranchVal(link);
			branchval[link->GetBranch()] = v;
			RecursiveCreate(link->Out());
		}
	}

	virtual ~_BranchValPtrTree()	{}

	virtual void RecursiveDelete(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDelete(link->Out());
			delete this->GetBranchVal(link->GetBranch());
			this->SetBranchVal(link->GetBranch(),0);
		}
		if (WithRoot() && from->isRoot())	{
			delete this->GetBranchVal(0);
		}
	}

	friend ostream& operator<<(ostream& os, _BranchValPtrTree<V>& b)	{
		b.ToStream(os,b.GetRoot());
		os << '\n';
		return os;
	}

	friend istream& operator>>(istream& is, _BranchValPtrTree<V>& b)  {
		b.FromStream(is,b.GetRoot());
		return is;
	}

	void ToStream(ostream& os, const Link* from)	{
		if (WithRoot() && from->isRoot())	{
			os <<  *(this->GetBranchVal(0)) << '\t';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			os <<  *(this->GetBranchVal(link->GetBranch())) << '\t';
			ToStream(os, link->Out());
		}
	}

	void FromStream(istream& is, Link* from)	{
		if (WithRoot() && from->isRoot())	{
			is >>  *(this->GetBranchVal(0));
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			is >>  *(this->GetBranchVal(link->GetBranch()));
			FromStream(is, link->Out());
		}
	}

	protected:

	virtual V* CreateBranchVal(const Link* link) = 0;

	map<const Branch*,V*> branchval;

	bool withRoot;
};

template<class T> class BranchValPtrTree : public virtual _BranchValPtrTree<T> {

};

template<class V> class BranchValPtrTree<Rvar<V> > : public virtual _BranchValPtrTree<Rvar<V> >, public virtual BranchVarTree<V>	{

	public:

	virtual Rvar<V>* GetBranchVal(const Branch* branch)	{
		return _BranchValPtrTree<Rvar<V> >::GetBranchVal(branch);
	}
};

template<class V> class BranchValPtrTree<Dvar<V> > : public virtual _BranchValPtrTree<Dvar<V> >, public virtual BranchVarTree<V>	{

	public:

	virtual Dvar<V>* GetBranchVal(const Branch* branch)	{
		return _BranchValPtrTree<Dvar<V> >::GetBranchVal(branch);
	}
};

template<> class BranchValPtrTree< Rvar<PosReal> > : public virtual _BranchValPtrTree<Rvar<PosReal> >, public virtual RandomLengthTree {

	public:

	virtual Rvar<PosReal>* GetBranchVal(const Branch* branch)	{
		return _BranchValPtrTree<Rvar<PosReal> >::GetBranchVal(branch);
	}
};

template<> class BranchValPtrTree< Dvar<PosReal> > : public virtual _BranchValPtrTree<Dvar<PosReal> >, public virtual LengthTree {

	public:

	virtual Dvar<PosReal>* GetBranchVal(const Branch* branch)	{
		return _BranchValPtrTree<Dvar<PosReal> >::GetBranchVal(branch);
	}
};

template<class U> class _NodeValPtrTree   : public virtual AbstractTree {

	public:

	U* GetNodeVal(const Node* node)	{
		return nodeval[node];
	}

	void SetNodeVal(const Node* node, U* in)	{
		nodeval[node] = in;
	}

	virtual ~_NodeValPtrTree()	{}

	virtual void RecursiveCreate(const Link* from)	{
		U* u = CreateNodeVal(from);
		nodeval[from->GetNode()] = u;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out());
		}
	}

	virtual void RecursiveDelete(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDelete(link->Out());
		}
		delete this->GetNodeVal(from->GetNode());
		this->SetNodeVal(from->GetNode(),0);
	}

	friend ostream& operator<<(ostream& os, _NodeValPtrTree<U>& b)	{
		b.ToStream(os,b.GetRoot());
		return os;
	}

	friend istream& operator>>(istream& is, _NodeValPtrTree<U>& b)  {
		b.FromStream(is,b.GetRoot());
		return is;
	}

	void ToStream(ostream& os, const Link* from)	{
		os <<  *(this->GetNodeVal(from->GetNode())) << '\n';
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStream(os, link->Out());
		}
	}

	void FromStream(istream& is, Link* from)	{
		is >>  *(this->GetNodeVal(from->GetNode()));
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			FromStream(is, link->Out());
		}
	}


	protected:

	virtual U* CreateNodeVal(const Link* link) = 0;

	map<const Node*,U*> nodeval;

};

template<class U> class NodeValPtrTree : public virtual _NodeValPtrTree<U> {

};

template<class V> class NodeValPtrTree<Rvar<V> > : public virtual _NodeValPtrTree<Rvar<V> >, public virtual NodeVarTree<V>	{

	public:

	virtual Rvar<V>* GetNodeVal(const Node* node)	{
		return _NodeValPtrTree<Rvar<V> >::GetNodeVal(node);
	}
};

template<class V> class NodeValPtrTree<Dvar<V> > : public virtual _NodeValPtrTree<Dvar<V> >, public virtual NodeVarTree<V>	{

	public:

	virtual Dvar<V>* GetNodeVal(const Node* node)	{
		return _NodeValPtrTree<Dvar<V> >::GetNodeVal(node);
	}
};

template<class U, class V> class NodeBranchValPtrTree : public virtual NodeValPtrTree<U>, public virtual BranchValPtrTree<V> {

	public:

	virtual ~NodeBranchValPtrTree()	{}

	// by default: makes 2 sweeps, for creating nodes, and then branches
	virtual void RecursiveCreate(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveCreate: should not be called on root\n";
			exit(1);
		}
		RecursiveCreateNode(from);
		RecursiveCreateBranch(from);
	}

	virtual void RecursiveCreateNode(Link* from)	{
		U* u = this->CreateNodeVal(from);
		this->nodeval[from->GetNode()] = u;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreateNode(link->Out());
		}
	}

	virtual void RecursiveCreateBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			V* v= this->CreateBranchVal(link);
			this->branchval[link->GetBranch()] = v;
			RecursiveCreateBranch(link->Out());
		}
	}

	virtual void RecursiveDeleteNode(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteNode(link->Out());
		}
		delete this->nodeval[from->GetNode()];
	}

	virtual void RecursiveDeleteBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteBranch(link->Out());
			delete this->branchval[link->GetBranch()];
		}
	}

	virtual void RecursiveDelete(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveDelete: should not be called on root\n";
			exit(1);
		}
		RecursiveDeleteBranch(from);
		RecursiveDeleteNode(from);
	}


	void ToStream(ostream& os)	{
		NodeValPtrTree<U>::ToStream(os,this->GetRoot());
		os << '\n';
		BranchValPtrTree<V>::ToStream(os,this->GetRoot());
	}

	friend ostream& operator<<(ostream& os, NodeBranchValPtrTree<U,V>& b)	{
		b.ToStream(os);
		return os;
	}

	void FromStream(istream& is)	{
		NodeValPtrTree<U>::FromStream(is,this->GetRoot());
		BranchValPtrTree<V>::FromStream(is,this->GetRoot());
	}

	friend istream& operator>>(istream& is, NodeBranchValPtrTree<U,V>& b)  {
		b.FromStream(is);
		return is;
	}
};

template<class U, class V> class NodeBranchValPtrTree<Rvar<U>, Dvar<V> > : public virtual NodeValPtrTree<Rvar<U> >, public virtual BranchValPtrTree<Dvar<V> >, public NodeBranchVarTree<U,V>  {

	public:

	virtual ~NodeBranchValPtrTree()	{}

	// by default: makes 2 sweeps, for creating nodes, and then branches
	virtual void RecursiveCreate(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveCreate: should not be called on root\n";
			exit(1);
		}
		RecursiveCreateNode(from);
		RecursiveCreateBranch(from);
	}

	virtual void RecursiveCreateNode(Link* from)	{
		Rvar<U>* u = this->CreateNodeVal(from);
		this->nodeval[from->GetNode()] = u;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreateNode(link->Out());
		}
	}

	virtual void RecursiveCreateBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			Dvar<V>* v= this->CreateBranchVal(link);
			this->branchval[link->GetBranch()] = v;
			RecursiveCreateBranch(link->Out());
		}
	}

	virtual void RecursiveDeleteNode(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteNode(link->Out());
		}
		delete this->nodeval[from->GetNode()];
	}

	virtual void RecursiveDeleteBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteBranch(link->Out());
			delete this->branchval[link->GetBranch()];
		}
	}

	virtual void RecursiveDelete(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveDelete: should not be called on root\n";
			exit(1);
		}
		RecursiveDeleteBranch(from);
		RecursiveDeleteNode(from);
	}


	void ToStream(ostream& os)	{
		NodeValPtrTree<Rvar<U> >::ToStream(os,this->GetRoot());
		os << '\n';
		BranchValPtrTree<Dvar<V> >::ToStream(os,this->GetRoot());
	}

	friend ostream& operator<<(ostream& os, NodeBranchValPtrTree<Rvar<U> ,Dvar<V> >& b)	{
		b.ToStream(os);
		return os;
	}

	void FromStream(istream& is)	{
		NodeValPtrTree<Rvar<U> >::FromStream(is,this->GetRoot());
		BranchValPtrTree<Dvar<V> >::FromStream(is,this->GetRoot());
	}

	friend istream& operator>>(istream& is, NodeBranchValPtrTree<Rvar<U> ,Dvar<V> >& b)  {
		b.FromStream(is);
		return is;
	}
};

template<class U, class V> class NodeBranchValPtrTree<Dvar<U>, Dvar<V> > : public virtual NodeValPtrTree<Dvar<U> >, public virtual BranchValPtrTree<Dvar<V> >, public NodeBranchVarTree<U,V>  {

	public:

	virtual ~NodeBranchValPtrTree()	{}

	// by default: makes 2 sweeps, for creating nodes, and then branches
	virtual void RecursiveCreate(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveCreate: should not be called on root\n";
			exit(1);
		}
		RecursiveCreateNode(from);
		RecursiveCreateBranch(from);
	}

	virtual void RecursiveCreateNode(Link* from)	{
		Dvar<U>* u = this->CreateNodeVal(from);
		this->nodeval[from->GetNode()] = u;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreateNode(link->Out());
		}
	}

	virtual void RecursiveCreateBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			Dvar<V>* v= this->CreateBranchVal(link);
			this->branchval[link->GetBranch()] = v;
			RecursiveCreateBranch(link->Out());
		}
	}

	virtual void RecursiveDeleteNode(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteNode(link->Out());
		}
		delete this->nodeval[from->GetNode()];
	}

	virtual void RecursiveDeleteBranch(Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeleteBranch(link->Out());
			delete this->branchval[link->GetBranch()];
		}
	}

	virtual void RecursiveDelete(Link* from)	{
		if (! from->isRoot())	{
			cerr << "error in NodeBranchValPtrTree::RecursiveDelete: should not be called on root\n";
			exit(1);
		}
		RecursiveDeleteBranch(from);
		RecursiveDeleteNode(from);
	}


	void ToStream(ostream& os)	{
		NodeValPtrTree<Dvar<U> >::ToStream(os,this->GetRoot());
		os << '\n';
		BranchValPtrTree<Dvar<V> >::ToStream(os,this->GetRoot());
	}

	friend ostream& operator<<(ostream& os, NodeBranchValPtrTree<Dvar<U> ,Dvar<V> >& b)	{
		b.ToStream(os);
		return os;
	}

	void FromStream(istream& is)	{
		NodeValPtrTree<Dvar<U> >::FromStream(is,this->GetRoot());
		BranchValPtrTree<Dvar<V> >::FromStream(is,this->GetRoot());
	}

	friend istream& operator>>(istream& is, NodeBranchValPtrTree<Dvar<U> ,Dvar<V> >& b)  {
		b.FromStream(is);
		return is;
	}
};


#endif
