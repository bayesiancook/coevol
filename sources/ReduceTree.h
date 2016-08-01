#ifndef TREECUTER_H
#define TREECUTER_H

#include "Chronogram.h"



/* Create copy of Tree where and link corresponding branches */
template<class T> class SubBranchTree : public virtual BranchValPtrTree<T>{

	Tree* tree;
	BranchValPtrTree<T>* originaltree;
	map<Branch*,Branch*> tofrom;

	public :

	SubBranchTree(Tree* intree, BranchValPtrTree<T>* inoriginaltree, map<Branch*,Branch*>& intofrom){
		tree=intree;
		originaltree=inoriginaltree;
		tofrom=intofrom;
		SetWithRoot(originaltree->WithRoot());
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree(){return tree;}
	Link* GetRoot() {return GetTree()->GetRoot();}


	private :

	T* CreateBranchVal(const Link* link){
		return originaltree->GetBranchVal(tofrom[link->GetBranch()]);
	}

};

/* Create copy of Tree where and link corresponding branches */
template<class T> class SubNodeTree : public virtual NodeValPtrTree<T>{

	Tree* tree;
	NodeValPtrTree<T>* originaltree;
	map<Node*,Node*> tofrom;

	public :

	SubNodeTree(Tree* intree, NodeValPtrTree<T>* inoriginaltree, map<Node*,Node*>& intofrom){
		tree=intree;
		originaltree=inoriginaltree;
		tofrom=intofrom;
		RecursiveCreate(GetRoot());
	}

	Tree* GetTree(){return tree;}
	Link* GetRoot() {return GetTree()->GetRoot();}


	private :

	T* CreateNodeVal(const Link* link){
		return originaltree->GetNodeVal(tofrom[link->GetNode()]);
	}

};


class ReduceTree {

	bool withunarynode;

	public :

	ReduceTree(Tree* intreefrom, TaxonSet* taxonset, bool inwithunarynode = true){
		treefrom=intreefrom;
		withunarynode = inwithunarynode;
		tree = new Tree(treefrom);
		RecursiveSetToFrom(tree->GetRoot(),treefrom->GetRoot());
		int tot =0;
		tree->RegisterWith(taxonset, tree->GetRoot(), tot);
	}

	~ReduceTree(){}


	template<class T>
	BranchValPtrTree<T>* GetSubBranchTree(BranchValPtrTree<T>* intree){
		if(!withunarynode){
			cerr << "Cannot call ReduceTree::GetSubBranchTree with no unary node\n";
			exit(0);
		}
		if(treefrom != intree->GetTree()){
			cerr << "bad call of GetSubBranchTree";
			exit(0);
		}
		return new SubBranchTree<T>(tree, intree, branchtofrom);
	}


	template<class T>
	NodeValPtrTree<T>* GetSubNodeTree( NodeValPtrTree<T>* intree ){
		if(treefrom != intree->GetTree()){
			cerr << "bad call of GetSubNodeTree";
			exit(0);
		}
		return new SubNodeTree<T>(tree, intree, nodetofrom);
	}




	private :
	Tree* treefrom;
	Tree* tree;

	map<Branch*,Branch*> branchtofrom;
	map<Node*,Node*> nodetofrom;

	//Fill a map king each branch of the actual tree to a branch in the original length tree
	void RecursiveSetToFrom( const Link* fromto , const Link* fromfrom ){
		branchtofrom[fromto->GetBranch()]=fromfrom->GetBranch();
		nodetofrom[fromto->GetNode()]=fromfrom->GetNode();
		const Link* linkto = fromto->Next();
		for( const Link* linkfrom = fromfrom->Next() ; linkfrom != fromfrom ; linkfrom = linkfrom->Next() ){
			RecursiveSetToFrom( linkto->Out() , linkfrom->Out() );
			linkto = linkto->Next();
		}
	}
};





#endif //TREECUTER_H
