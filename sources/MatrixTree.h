#ifndef MATRIXTREE_H
#define MATRIXTREE_H

#include "GTRSubMatrix.h"

#include "CodonSequenceAlignment.h"
#include "CodonSubMatrix.h"
#include "MG3OmegaCodonSubMatrix.h"
#include "MGOmega3CodonSubMatrix.h"
#include "GCProcess.h"

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree, Var<PosReal>* inrootomega, map<const Link*, int>* inmissingmap = 0) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		missingmap = inmissingmap;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if ((missingmap) && ((*missingmap)[link]))	{
			return 0;
		}
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;
	map<const Link*, int>* missingmap;

};




class GCMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GCMatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	NucMatrixTree* nucmatrixtree;
	Var<PosReal>* rootomega;

};

class GC3MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GC3MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree1, NucMatrixTree* innucmatrixtree2, NucMatrixTree* innucmatrixtree3, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrixtree1 = innucmatrixtree1;
		nucmatrixtree2 = innucmatrixtree2;
		nucmatrixtree3 = innucmatrixtree3;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GC3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMG3OmegaCodonSubMatrix(statespace,nucmatrixtree1->GetBranchVal(link->GetBranch()),nucmatrixtree2->GetBranchVal(link->GetBranch()),nucmatrixtree3->GetBranchVal(link->GetBranch()),rootomega);
		}
		return new RandomMG3OmegaCodonSubMatrix(statespace,nucmatrixtree1->GetBranchVal(link->GetBranch()),nucmatrixtree2->GetBranchVal(link->GetBranch()),nucmatrixtree3->GetBranchVal(link->GetBranch()),omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	NucMatrixTree* nucmatrixtree1;
	NucMatrixTree* nucmatrixtree2;
	NucMatrixTree* nucmatrixtree3;
	Var<PosReal>* rootomega;

};

class Omega3MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	Omega3MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatstree, LengthTree* inomegatv0tree, LengthTree* inomegatvgctree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		omegatvgctree = inomegatvgctree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~Omega3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega3CodonSubMatrix(statespace,nucmatrix,rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega3CodonSubMatrix(statespace,nucmatrix,omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()),omegatvgctree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	LengthTree* omegatvgctree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};

class Omega3X3MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	Omega3X3MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* insynratetstree, LengthTree* insynratetvgctree, LengthTree* inomegatstree, LengthTree* inomegatv0tree, LengthTree* inomegatvgctree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		synratetvgctree = insynratetvgctree;
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		omegatvgctree = inomegatvgctree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~Omega3X3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega3X3CodonSubMatrix(statespace,nucmatrix,rootomega,rootomega,rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega3X3CodonSubMatrix(statespace,nucmatrix,synratetstree->GetBranchVal(link->GetBranch()), synratetvgctree->GetBranchVal(link->GetBranch()), omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()),omegatvgctree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* synratetstree;
	LengthTree* synratetvgctree;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	LengthTree* omegatvgctree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};

class Omega2X2MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	Omega2X2MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* insynratetstree, LengthTree* inomegatstree, LengthTree* inomegatv0tree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~Omega2X2MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega2X2CodonSubMatrix(statespace,nucmatrix,rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega2X2CodonSubMatrix(statespace,nucmatrix,synratetstree->GetBranchVal(link->GetBranch()), omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* synratetstree;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};

class GCOmega3MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GCOmega3MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* inomegatstree, LengthTree* inomegatv0tree, LengthTree* inomegatvgctree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		omegatvgctree = inomegatvgctree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GCOmega3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega3CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega3CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()),omegatvgctree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	LengthTree* omegatvgctree;
	NucMatrixTree* nucmatrixtree;
	Var<PosReal>* rootomega;

};

class GCOmega3X3MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GCOmega3X3MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* insynratetstree, LengthTree* insynratetvgctree, LengthTree* inomegatstree, LengthTree* inomegatv0tree, LengthTree* inomegatvgctree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		synratetvgctree = insynratetvgctree;
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		omegatvgctree = inomegatvgctree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GCOmega3X3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega3X3CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega,rootomega,rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega3X3CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),synratetstree->GetBranchVal(link->GetBranch()),synratetvgctree->GetBranchVal(link->GetBranch()), omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()),omegatvgctree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* synratetstree;
	LengthTree* synratetvgctree;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	LengthTree* omegatvgctree;
	NucMatrixTree* nucmatrixtree;
	Var<PosReal>* rootomega;

};

class GCOmega2X2MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	GCOmega2X2MatrixTree(CodonStateSpace* instatespace, NucMatrixTree* innucmatrixtree, LengthTree* insynratetstree, LengthTree* inomegatstree, LengthTree* inomegatv0tree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		omegatstree = inomegatstree;
		omegatv0tree = inomegatv0tree;
		nucmatrixtree = innucmatrixtree;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~GCOmega2X2MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmega2X2CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),rootomega,rootomega,rootomega);
		}
		return new RandomMGOmega2X2CodonSubMatrix(statespace,nucmatrixtree->GetBranchVal(link->GetBranch()),synratetstree->GetBranchVal(link->GetBranch()), omegatstree->GetBranchVal(link->GetBranch()),omegatv0tree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatstree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* synratetstree;
	LengthTree* omegatstree;
	LengthTree* omegatv0tree;
	NucMatrixTree* nucmatrixtree;
	Var<PosReal>* rootomega;

};


#endif
