#ifndef BGC_H
#define BGC_H

#include "BranchProductProcess.h"

#include "BGCSubMatrix.h"
#include "BGCCpGSubMatrix.h"

int BGCMutSelSubMatrix::bgccount = 0;
int BGCCpGMutSelSubMatrix::bgccount = 0;

class BGCInstantVal : public Dvar<PosReal>	{

	public:

	BGCInstantVal(Var<RealVector>* inlogbgcup, Var<RealVector>* inlogbgcdown, Var<Real>* inlogbgcoffset, int inindex = 0)	{
		logbgcup = inlogbgcup;
		logbgcdown = inlogbgcdown;
		logbgcoffset = inlogbgcoffset;
		index = index;
		Register(logbgcup);
		Register(logbgcdown);
		Register(logbgcoffset);
		specialUpdate();
	}

	~BGCInstantVal() {}

	void specialUpdate()	{
		double bgcup = exp(logbgcoffset->val() + (*logbgcup)[index]);
		double bgcdown = exp(logbgcoffset->val() + (*logbgcdown)[index]);
		setval(0.5 * (bgcup + bgcdown));
	}

	protected:

	Var<RealVector>* logbgcup;
	Var<RealVector>* logbgcdown;
	Var<Real>* logbgcoffset;
	int index;

};

class BGCProcess : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	BGCProcess(NodeVarTree<RealVector>* inprocess, Var<Real>* inlogbgcoffset, int inindex = 0)	{
		index = inindex;
		process = inprocess;
		logbgcoffset = inlogbgcoffset;
		SetWithRoot(false);
		RecursiveCreate(GetRoot());
	}

	~BGCProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	BGCInstantVal* GetBGCInstantVal(const Branch* branch)	{
		return dynamic_cast<BGCInstantVal*>(GetBranchVal(branch));
	}

	protected:

	BGCInstantVal* CreateBranchVal(const Link* link)	{
		return new BGCInstantVal(process->GetNodeVal(link->GetNode()),process->GetNodeVal(link->Out()->GetNode()),logbgcoffset,index);
	}

	void RecursiveSpecialUpdate(const Link* from)	{
		if (! from->isRoot())	{
			GetBGCInstantVal(from->GetBranch())->specialUpdate();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	private	:

	NodeVarTree<RealVector>* process;
	int index;
	Var<Real>* logbgcoffset;
};

class MutMatrixTree : public NucMatrixTree	{


	public:

	MutMatrixTree() {}

	MutMatrixTree(BranchVarTree<PosReal>* inlambdaprocess, Var<PosReal>* inlambda, Var<PosReal>* inat2cg, Var<PosReal>* inat2gc, Var<PosReal>* inat2ta, Var<PosReal>* ingc2cg, bool innormalise)	{
		lambdaprocess = inlambdaprocess;
		lambda = inlambda;
		at2cg = inat2cg;
		at2gc = inat2gc;
		at2ta = inat2ta;
		gc2cg = ingc2cg;
		normalise = innormalise;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~MutMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return lambdaprocess->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		/*
		if (link->isRoot())	{
		}
		*/
		return new RandomMutSubMatrix(lambda,at2cg,at2gc,at2ta,gc2cg,true,lambdaprocess->GetBranchVal(link->GetBranch()));
	}

	protected :

	BranchVarTree<PosReal>* lambdaprocess;
	Var<PosReal>* lambda;
	Var<PosReal>* at2cg;
	Var<PosReal>* at2gc;
	Var<PosReal>* at2ta;
	Var<PosReal>* gc2cg;
	bool normalise;
};

class BGCMatrixTree: public NucMatrixTree	{

	public:

	BGCMatrixTree() {}

	BGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise, int indiscgam, NucMatrixTree* inmutmatrixtree = 0)	{
	// BGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, Var<PosReal>* inrootbgc, bool innormalise)	{
		process = inprocess;
		mutmatrix = inmutmatrix;
		mutmatrixtree = inmutmatrixtree;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		discgam = indiscgam;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~BGCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	int GetBGCOverflowCount()	{
		return RecursiveGetBGCOverflowCount(GetRoot());
	}

	BGCMutSelSubMatrix* GetBGCMutSelSubMatrix(const Branch* branch)	{
		BGCMutSelSubMatrix*  tmp = dynamic_cast<BGCMutSelSubMatrix*>(GetBranchVal(branch));
		if (!tmp)	{
			cerr << "error in BGCMatrixTree: " << GetBranchVal(branch) << '\n';
			exit(1);
		}
		return tmp;
	}

	Tree* GetTree() {return process->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (mutmatrix)	{
			if (link->isRoot())	{
				// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
				return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,0,normalise,discgam);
			}
			// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
			return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),alpha,normalise,discgam);
		}
		else	{
			if (link->isRoot())	{
				return new RandomBGCMutSelSubMatrix(mutmatrixtree->GetBranchVal(link->GetBranch()),rootbgc,0,normalise,discgam);
			}
			return new RandomBGCMutSelSubMatrix(mutmatrixtree->GetBranchVal(link->GetBranch()),process->GetBranchVal(link->GetBranch()),alpha,normalise,discgam);
		}
	}

	int RecursiveGetBGCOverflowCount(const Link* from)	{
		int tot = GetBGCMutSelSubMatrix(from->GetBranch())->GetBGCOverflowCount();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveGetBGCOverflowCount(link->Out());
		}
		return tot;
	}

	protected :

	RandomSubMatrix* mutmatrix;
	NucMatrixTree* mutmatrixtree;
	BranchVarTree<PosReal>* process;
	Var<PosReal>* rootbgc;
	Var<PosReal>* alpha;
	bool normalise;
	int discgam;
};

class NonRevBGCMatrixTree: public BGCMatrixTree	{

	public:

	NonRevBGCMatrixTree() {}

	NonRevBGCMatrixTree(LengthTree* inlengthtree, BranchVarTree<PosReal>* inprocess, RandomNonRevSubMatrix* innonrevmutmatrix, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise, int indiscgam, int indiscn)	{
		lengthtree = inlengthtree;
		process = inprocess;
		nonrevmutmatrix = innonrevmutmatrix;
		mutmatrix = innonrevmutmatrix;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		discn = indiscn;
		discgam = indiscgam;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~NonRevBGCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
			return new RandomBGCNonRevMutSelSubMatrix(nonrevmutmatrix, lengthtree->GetBranchVal(link->GetBranch()),rootbgc,0,normalise,discgam,discn);
		}
		// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
		return new RandomBGCNonRevMutSelSubMatrix(nonrevmutmatrix, lengthtree->GetBranchVal(link->GetBranch()),process->GetBranchVal(link->GetBranch()),alpha,normalise,discgam,discn);
	}

	int RecursiveGetBGCOverflowCount(const Link* from)	{
		int tot = GetBGCMutSelSubMatrix(from->GetBranch())->GetBGCOverflowCount();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveGetBGCOverflowCount(link->Out());
		}
		return tot;
	}

	protected :

	RandomNonRevSubMatrix* nonrevmutmatrix;
	int discn;
	int discgam;
	LengthTree* lengthtree;
};

class AlphaBGCMatrixTree: public BGCMatrixTree	{

	public:

	AlphaBGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, BranchVarTree<PosReal>* inalphaprocess, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise)	{
	// BGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, Var<PosReal>* inrootbgc, bool innormalise)	{
		process = inprocess;
		mutmatrix = inmutmatrix;
		alphaprocess = inalphaprocess;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~AlphaBGCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
			return new RandomAlphaBGCMutSelSubMatrix(mutmatrix,rootbgc,0,0,normalise);
		}
		// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
		return new RandomAlphaBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),alphaprocess->GetBranchVal(link->GetBranch()),alpha,normalise);
	}

	private	:

	BranchVarTree<PosReal>* alphaprocess;
};

class HotSpotBGCMatrixTree: public BGCMatrixTree	{

	public:

	HotSpotBGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, BranchVarTree<PosReal>* inalphaprocess, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise)	{
		process = inprocess;
		mutmatrix = inmutmatrix;
		alphaprocess = inalphaprocess;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~HotSpotBGCMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
			return new RandomHotSpotBGCMutSelSubMatrix(mutmatrix,rootbgc,0,0,normalise);
		}
		// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
		return new RandomHotSpotBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),alphaprocess->GetBranchVal(link->GetBranch()),alpha,normalise);
	}

	private	:

	BranchVarTree<PosReal>* alphaprocess;
};


class BGCCpGMatrixTree: public NucMatrixTree	{

	public:

	BGCCpGMatrixTree() {}

	BGCCpGMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, RandomSubMatrix* inctmutmatrix, Var<PosReal>* incpgrate, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise, int indiscgam)	{
	// BGCMatrixTree(BranchVarTree<PosReal>* inprocess, RandomSubMatrix* inmutmatrix, Var<PosReal>* inrootbgc, bool innormalise)	{
		process = inprocess;
		mutmatrix = inmutmatrix;
		ctmutmatrix = inctmutmatrix;
		cpgrate = incpgrate;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		discgam = indiscgam;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~BGCCpGMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	int GetBGCOverflowCount()	{
		return RecursiveGetBGCOverflowCount(GetRoot());
	}

	BGCCpGMutSelSubMatrix* GetBGCCpGMutSelSubMatrix(const Branch* branch)	{
		BGCCpGMutSelSubMatrix*  tmp = dynamic_cast<BGCCpGMutSelSubMatrix*>(GetBranchVal(branch));
		if (!tmp)	{
			cerr << "error in BGCCpGMatrixTree: " << GetBranchVal(branch) << '\n';
			exit(1);
		}
		return tmp;
	}

	Tree* GetTree() {return process->GetTree();}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
			return new RandomBGCCpGMutSelSubMatrix(mutmatrix,ctmutmatrix,cpgrate,rootbgc,0,normalise,discgam);
		}
		// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
		return new RandomBGCCpGMutSelSubMatrix(mutmatrix,ctmutmatrix,cpgrate,process->GetBranchVal(link->GetBranch()),alpha,normalise,discgam);
	}

	int RecursiveGetBGCOverflowCount(const Link* from)	{
		int tot = GetBGCCpGMutSelSubMatrix(from->GetBranch())->GetBGCOverflowCount();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveGetBGCOverflowCount(link->Out());
		}
		return tot;
	}

	protected :

	RandomSubMatrix* mutmatrix;
	RandomSubMatrix* ctmutmatrix;
	Var<PosReal>* cpgrate;
	BranchVarTree<PosReal>* process;
	Var<PosReal>* rootbgc;
	Var<PosReal>* alpha;
	bool normalise;
	int discgam;
};

class NonRevBGCCpGMatrixTree: public BGCCpGMatrixTree	{

	public:

	NonRevBGCCpGMatrixTree() {}

	NonRevBGCCpGMatrixTree(LengthTree* inlengthtree, BranchVarTree<PosReal>* inprocess, RandomNonRevSubMatrix* innonrevmutmatrix, RandomNonRevSubMatrix* inctnonrevmutmatrix, Var<PosReal>* incpgrate, Var<PosReal>* inalpha, Var<PosReal>* inrootbgc, bool innormalise, int indiscgam, int indiscn)	{
		lengthtree = inlengthtree;
		process = inprocess;
		nonrevmutmatrix = innonrevmutmatrix;
		ctnonrevmutmatrix = inctnonrevmutmatrix;
		cpgrate = incpgrate;
		mutmatrix = innonrevmutmatrix;
		ctmutmatrix = inctnonrevmutmatrix;
		alpha = inalpha;
		rootbgc = inrootbgc;
		normalise = innormalise;
		discn = indiscn;
		discgam = indiscgam;
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	~NonRevBGCCpGMatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			// return new RandomBGCMutSelSubMatrix(mutmatrix,rootbgc,normalise);
			return new RandomBGCCpGNonRevMutSelSubMatrix(nonrevmutmatrix, ctnonrevmutmatrix, cpgrate, lengthtree->GetBranchVal(link->GetBranch()),rootbgc,0,normalise,discgam,discn);
		}
		// return new RandomBGCMutSelSubMatrix(mutmatrix,process->GetBranchVal(link->GetBranch()),normalise);
		return new RandomBGCCpGNonRevMutSelSubMatrix(nonrevmutmatrix, ctnonrevmutmatrix, cpgrate, lengthtree->GetBranchVal(link->GetBranch()),process->GetBranchVal(link->GetBranch()),alpha,normalise,discgam,discn);
	}

	/*
	int RecursiveGetBGCOverflowCount(const Link* from)	{
		int tot = GetBGCMutSelSubMatrix(from->GetBranch())->GetBGCOverflowCount();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveGetBGCOverflowCount(link->Out());
		}
		return tot;
	}
	*/

	protected :

	RandomNonRevSubMatrix* nonrevmutmatrix;
	RandomNonRevSubMatrix* ctnonrevmutmatrix;
	int discn;
	int discgam;
	LengthTree* lengthtree;
};

#endif

