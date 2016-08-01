#ifndef GCBRANCHPROCESS_H
#define GCBRANCHPROCESS_H

#include "RandomTypes.h"
#include "Normal.h"
#include "ValTree.h"
#include "GTRSubMatrix.h"
#include "HKYSubMatrix.h"
#include "NucSubMatrix.h"
#include "BranchProcess.h"

class InstantStat : public Dvar<Profile>	{

	public:

	InstantStat(Var<UnitReal>* ingc)	{
		setval(Profile(Nnuc));
		bkvalue = Profile(Nnuc);
		gc = ingc;
		Register(gc);
		specialUpdate();
	}

	void specialUpdate()	{
		(*this)[0] = (*this)[3] = 0.5 * (1 - gc->val());
		(*this)[1] = (*this)[2] = 0.5 * gc->val();
	}

	double GetGCContent()	{
		return gc->val();
	}

	private:

	Var<UnitReal>* gc;

};

class GCStatTree : public BranchValPtrTree<InstantStat>	{

	public:

	// for the root
	// creates a special
	GCStatTree(BranchVarTree<UnitReal>* inprocess, Var<UnitReal>* inrootval)	{
		rootval = inrootval;
		if (!rootval)	{
			cerr << "error: should specify a root value\n";
			exit(1);
		}
		SetWithRoot(true);
		process = inprocess;
		RecursiveCreate(GetRoot());
	}

	~GCStatTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	double GetMeanGCContent()	{
		int n = 0;
		double total = GetTotalGCContent(GetRoot(),n);
		return total / n;
	}

	double GetVarGCContent()	{
		int n = 0;
		double total1 = GetTotalGCContent(GetRoot(),n);
		n = 0;
		double total2 = GetTotalSquareGCContent(GetRoot(),n);
		total1 /= n;
		total2 /= n;
		total2 -= total1 * total1;
		return total2;
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	protected:

	void RecursiveSpecialUpdate(const Link* from)	{
		GetBranchVal(from->GetBranch())->specialUpdate();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	double GetTotalGCContent(const Link* from, int& n)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalGCContent(link->Out(),n);
		}
		total += GetBranchVal(from->GetBranch())->GetGCContent();
		n++;
		return total;
	}

	double GetTotalSquareGCContent(const Link* from, int& n)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareGCContent(link->Out(),n);
		}
		double tmp = GetBranchVal(from->GetBranch())->GetGCContent();
		total += tmp * tmp;
		n++;
		return total;
	}

	InstantStat* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new InstantStat(rootval);
		}
		return new InstantStat(process->GetBranchVal(link->GetBranch()));
	}

	private:
	BranchVarTree<UnitReal>* process;
	Var<UnitReal>* rootval;

};

class NucMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:
	NucMatrixTree() {}

};

class GTRGCNucMatrixTree : public NucMatrixTree	{


	public:

	GTRGCNucMatrixTree(Var<Profile>* inrelrate, GCStatTree* instattree, bool innormalise) {
		SetWithRoot(instattree->WithRoot());
		stattree = instattree;
		relrate = inrelrate;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~GTRGCNucMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return stattree->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),normalise);
	}

	private:

	GCStatTree* stattree;
	Var<Profile>* relrate;
	bool normalise;

};

class HKYNucMatrixTree : public NucMatrixTree	{


	public:

	HKYNucMatrixTree(BranchVarTree<PosReal>* inkappatree, Var<Profile>* instat, Var<Profile>* intsrr, Var<Profile>* intvrr, bool innormalise) {
		SetWithRoot(true);
		stat = instat;
		kappatree = inkappatree;
		tsrr = intsrr;
		tvrr = intvrr;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~HKYNucMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return kappatree->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new HKYRandomSubMatrix(kappatree->GetBranchVal(link->GetBranch()),stat,tsrr,tvrr,normalise);
	}

	private:

	Var<Profile>* stat;
	BranchVarTree<PosReal>* kappatree;
	Var<Profile>* tsrr;
	Var<Profile>* tvrr;
	bool normalise;

};

class HKYGCNucMatrixTree : public NucMatrixTree	{


	public:

	HKYGCNucMatrixTree(BranchVarTree<PosReal>* inkappatree, GCStatTree* instattree, Var<Profile>* intsrr, Var<Profile>* intvrr, bool innormalise) {
		SetWithRoot(instattree->WithRoot());
		stattree = instattree;
		kappatree = inkappatree;
		tsrr = intsrr;
		tvrr = intvrr;

		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~HKYGCNucMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return stattree->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		return new HKYRandomSubMatrix(kappatree->GetBranchVal(link->GetBranch()),stattree->GetBranchVal(link->GetBranch()),tsrr,tvrr,normalise);
	}

	private:

	GCStatTree* stattree;
	BranchVarTree<PosReal>* kappatree;
	Var<Profile>* tsrr;
	Var<Profile>* tvrr;
	bool normalise;

};


class Nuc3X3MatrixTree : public NucMatrixTree	{


	public:

	Nuc3X3MatrixTree(Var<Profile>* inrr, Var<Profile>* instat, LengthTree* insynratetstree, LengthTree* insynratetvgctree, bool innormalise) {
	// Nuc3X3MatrixTree(Var<Profile>* inrr, Var<Profile>* instat, LengthTree* insynratetstree, LengthTree* insynratetvgctree, Var<PosReal>* inrootval) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		synratetvgctree = insynratetvgctree;
		rr = inrr;
		stat = instat;
		// rootval = inrootval;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~Nuc3X3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		/*
		if (link->isRoot())	{
			return new RandomNuc3X3SubMatrix(rootval,rootval,stat,rr,normalise);
		}
		*/
		return new RandomNuc3X3SubMatrix(synratetstree->GetBranchVal(link->GetBranch()), synratetvgctree->GetBranchVal(link->GetBranch()),stat,rr,normalise);
	}

	Tree* GetTree() {return synratetstree->GetTree();}

	private:

	LengthTree* synratetstree;
	LengthTree* synratetvgctree;
	// Var<PosReal>* rootval;
	Var<Profile>* rr;
	Var<Profile>* stat;
	bool normalise;

};

class Nuc2X2MatrixTree : public NucMatrixTree	{


	public:

	Nuc2X2MatrixTree(Var<Profile>* inrr, Var<Profile>* instat, LengthTree* insynratetstree, bool innormalise)	{ // , Var<PosReal>* inrootval) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		rr = inrr;
		stat = instat;
		// rootval = inrootval;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~Nuc2X2MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		/*
		if (link->isRoot())	{
			return new RandomNuc2X2SubMatrix(rootval,stat,rr);
		}
		*/
		return new RandomNuc2X2SubMatrix(synratetstree->GetBranchVal(link->GetBranch()),stat,rr,normalise);
	}

	Tree* GetTree() {return synratetstree->GetTree();}

	private:

	LengthTree* synratetstree;
	// Var<PosReal>* rootval;
	Var<Profile>* rr;
	Var<Profile>* stat;
	bool normalise;

};

class GCNuc3X3MatrixTree : public NucMatrixTree	{


	public:

	GCNuc3X3MatrixTree(Var<Profile>* inrr, GCStatTree* instattree, LengthTree* insynratetstree, LengthTree* insynratetvgctree, bool innormalise)	{ // , Var<PosReal>* inrootval) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		synratetvgctree = insynratetvgctree;
		rr = inrr;
		stattree = instattree;
		normalise = innormalise;
		// rootval = inrootval;
		RecursiveCreate(GetRoot());
	}

	~GCNuc3X3MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		/*
		if (link->isRoot())	{
			return new RandomNuc3X3SubMatrix(rootval,rootval,stattree->GetBranchVal(link->GetBranch()),rr);
		}
		*/
		return new RandomNuc3X3SubMatrix(synratetstree->GetBranchVal(link->GetBranch()), synratetvgctree->GetBranchVal(link->GetBranch()),stattree->GetBranchVal(link->GetBranch()),rr,normalise);
	}

	Tree* GetTree() {return synratetstree->GetTree();}

	private:

	LengthTree* synratetstree;
	LengthTree* synratetvgctree;
	// Var<PosReal>* rootval;
	Var<Profile>* rr;
	GCStatTree* stattree;
	bool normalise;

};

class GCNuc2X2MatrixTree : public NucMatrixTree	{


	public:

	GCNuc2X2MatrixTree(Var<Profile>* inrr, GCStatTree* instattree, LengthTree* insynratetstree, bool innormalise)	{ // Var<PosReal>* inrootval) {
		SetWithRoot(true);
		synratetstree = insynratetstree;
		rr = inrr;
		stattree = instattree;
		// rootval = inrootval;
		normalise = innormalise;
		RecursiveCreate(GetRoot());
	}

	~GCNuc2X2MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		/*
		if (link->isRoot())	{
			return new RandomNuc2X2SubMatrix(rootval,stattree->GetBranchVal(link->GetBranch()),rr);
		}
		*/
		return new RandomNuc2X2SubMatrix(synratetstree->GetBranchVal(link->GetBranch()), stattree->GetBranchVal(link->GetBranch()),rr,normalise);
	}

	Tree* GetTree() {return synratetstree->GetTree();}

	private:

	LengthTree* synratetstree;
	// Var<PosReal>* rootval;
	Var<Profile>* rr;
	GCStatTree* stattree;
	bool normalise;

};


#endif

