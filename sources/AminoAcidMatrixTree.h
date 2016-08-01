
#ifndef AMINOACIDMATRIXTREE_H
#define AMINOACIDMATRIXTREE_H

#include "AminoAcidOmegaSubMatrix.h"

class AminoAcidOmegaMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	AminoAcidOmegaMatrixTree(SimilarityMatrix* insimilaritymatrix, Var<Profile>* inrelrate, Var<Profile>* instationary, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		relrate = inrelrate;
		stationary = instationary;
		rootomega = inrootomega;
		similaritymatrix = insimilaritymatrix;
		RecursiveCreate(GetRoot());
	}

	~AminoAcidOmegaMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot()){
			return new RandomAAOmegaSubMatrix(similaritymatrix, relrate, stationary,rootomega);
		}
		return new RandomAAOmegaSubMatrix(similaritymatrix, relrate, stationary,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	LengthTree* omegatree;
	Var<Profile>* relrate;
	Var<Profile>* stationary;
	Var<PosReal>* rootomega;
	SimilarityMatrix * similaritymatrix;
};

class SplitAminoAcidOmegaMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	SplitAminoAcidOmegaMatrixTree(SplitAAMatrix* insplitaamat, int insplittype, SimilarityMatrix* insimilaritymatrix, Var<Profile>* inrelrate, Var<Profile>* instationary, LengthTree* intstvtree, LengthTree* inomegatstree, LengthTree* inomegatvtree, Var<PosReal>* inroottstv, Var<PosReal>* inrootomegats, Var<PosReal>* inrootomegatv) {
		SetWithRoot(true);
		tstvtree = intstvtree;
		omegatstree = inomegatstree;
		omegatvtree = inomegatvtree;
		relrate = inrelrate;
		stationary = instationary;
		roottstv = inroottstv;
		rootomegats = inrootomegats;
		rootomegatv = inrootomegatv;
		similaritymatrix = insimilaritymatrix;
		splitaaMatrix = insplitaamat;
		splittype = insplittype;
		cerr << "recursive create\n";
		RecursiveCreate(GetRoot());
		cerr << "rec ok\n";
	}

	~SplitAminoAcidOmegaMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot()){
			RandomSubMatrix* tmp = new RandomSplitAAOmegaSubMatrix(similaritymatrix, relrate, stationary, roottstv, rootomegats, rootomegatv, splitaaMatrix, splittype);
			return tmp;
		}
		if (omegatstree)	{
			return new RandomSplitAAOmegaSubMatrix(similaritymatrix, relrate, stationary, tstvtree->GetBranchVal(link->GetBranch()), omegatstree->GetBranchVal(link->GetBranch()), omegatvtree->GetBranchVal(link->GetBranch()), splitaaMatrix, splittype);
		}
		RandomSubMatrix* tmp = new RandomSplitAAOmegaSubMatrix(similaritymatrix, relrate, stationary, tstvtree->GetBranchVal(link->GetBranch()), 0, omegatvtree->GetBranchVal(link->GetBranch()), splitaaMatrix, splittype);
		return tmp;
	}

	Tree* GetTree() {return omegatvtree->GetTree();}

	private:

	LengthTree* tstvtree;
	LengthTree* omegatstree;
	LengthTree* omegatvtree;
	Var<Profile>* relrate;
	Var<Profile>* stationary;
	Var<PosReal>* roottstv;
	Var<PosReal>* rootomegats;
	Var<PosReal>* rootomegatv;
	SimilarityMatrix * similaritymatrix;
	SplitAAMatrix* splitaaMatrix;
	int splittype;
};

class AminoAcidDoubleOmegaMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	AminoAcidDoubleOmegaMatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, SimilarityMatrix* insimilaritymatrix, LengthTree* inomegadNdStree, LengthTree* inomegaKrKctree, Var<PosReal>* inrootomegadNdS, Var<PosReal>* inrootomegaKrKc) {
		SetWithRoot(true);
		omegadNdStree = inomegadNdStree;
		omegaKrKctree = inomegaKrKctree;
		rootomegadNdS = inrootomegadNdS;
		rootomegaKrKc = inrootomegaKrKc;
		statespace = instatespace;
		similaritymatrix = insimilaritymatrix;
		nucmatrix = innucmatrix;
		RecursiveCreate(GetRoot());
	}

	~AminoAcidDoubleOmegaMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot()){
			return new RandomAADoubleOmegaSubMatrix(statespace, nucmatrix, rootomegadNdS, rootomegaKrKc, similaritymatrix );
		}
		return new RandomAADoubleOmegaSubMatrix(statespace, nucmatrix, omegadNdStree->GetBranchVal(link->GetBranch()), omegaKrKctree->GetBranchVal(link->GetBranch()), similaritymatrix);
	}

	Tree* GetTree() {return omegaKrKctree->GetTree();}

	private:

	LengthTree* omegadNdStree;
	LengthTree* omegaKrKctree;
	Var<PosReal>* rootomegadNdS;
	Var<PosReal>* rootomegaKrKc;
	SimilarityMatrix * similaritymatrix;
	CodonStateSpace* statespace;
	RandomSubMatrix* nucmatrix;
};


class AAregMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	AAregMatrixTree(Var<Profile>* instationary, Var<RealVector>* inslope, Var<RealVector>* inintercept, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		stationary = instationary;
		rootomega = inrootomega;
		slope = inslope;
		intercept = inintercept;
		RecursiveCreate(GetRoot());
	}

	~AAregMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot()){
			return new RandomAAregSubMatrix(stationary, rootomega, slope, intercept);
		}
		return new RandomAAregSubMatrix(stationary, omegatree->GetBranchVal(link->GetBranch()), slope, intercept);

	}

	Tree* GetTree() {return omegatree->GetTree();}


	private:

	LengthTree* omegatree;
	Var<Profile>* stationary;
	Var<PosReal>* rootomega;
	Var<RealVector>* slope;
	Var<RealVector>* intercept;
};


#endif

