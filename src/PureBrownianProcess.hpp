
#ifndef PUREBROWNIANPROCESS_H
#define	PUREBROWNIANPROCESS_H

#include "ValTree.h"
#include "RandomBrownianPath.h"
#include "Chronogram.h"

class PureBrownianProcess : public BranchValPtrTree<Rvar<BrownianBridge> >, public MCMC
{
	public:
		PureBrownianProcess() {}
		PureBrownianProcess(Chronogram *intree, Var<CovMatrix> *insigma, Var<PosReal>* inagescale = 0, GlobalScalingFunction* inscalefunction = 0);
		virtual ~PureBrownianProcess();

		virtual Rvar<BrownianBridge> * CreateBranchVal(const Link* link);

		virtual double Move(double tuning);
		virtual double Move(Link* from, double tuning, int& count);

		virtual double GetLogProb();
		virtual double GetLogProb(Link* from);

		virtual void Clamp();
		virtual void Clamp(Link* from);

		virtual void drawSample();
		virtual void SampleBranch(Link* from);

		void RebuildSegments();
		void RebuildSegments(const Link* from);

		// various accessors
		virtual Tree* GetTree() {
			return GetLengthTree()->GetTree();
		}

		Chronogram* GetLengthTree()	{
			return tree;
		}

		RandomBrownianPath* GetRandomBrownianPath(const Link *from) {return dynamic_cast<RandomBrownianPath*>(GetBranchVal(from->GetBranch()));}

		void CheckLength() {
			CheckLength(GetRoot());
		}

		void RecursiveProposeMove(const Link* from, double tuning, CovMatrix& cov)	{
			if (! from->isRoot())	{
				GetRandomBrownianPath(from)->ProposeMoveCov(cov,tuning);
			}
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				RecursiveProposeMove(link->Out(),tuning,cov);
			}
		}

		void CheckLength(const Link* from) {
			for(const Link* link = from->Next(); link!=from; link = link->Next()) {
				CheckLength(link->Out());
				GetBranchVal(link->GetBranch())->checkLength(tree->GetNodeVal(link->GetNode())->val(), tree->GetNodeVal(link->Out()->GetNode())->val(), tree->GetBranchVal(link->GetBranch())->val());
			}
		}

		protected:

		Chronogram *tree;	   //The chronogram of the tree
		Var<CovMatrix> *sigma;	//covariance
		Var<PosReal>* agescale;
		GlobalScalingFunction* scalefunction;

};

class ConjugatePureBrownianProcess : public virtual PureBrownianProcess	{

	public:

	ConjugatePureBrownianProcess(Chronogram *intree, ConjugateInverseWishart *insigma, Var<PosReal>* inagescale = 0, GlobalScalingFunction* inscalefunction = 0);
	
	ConjugateRandomBrownianPath* GetConjugateRandomBrownianPath(const Link* link)	{
		ConjugateRandomBrownianPath* m = dynamic_cast<ConjugateRandomBrownianPath*> (GetBranchVal(link->GetBranch()));
		return m;
	}


	protected:

	virtual Rvar<BrownianBridge> * CreateBranchVal(const Link* link);

	ConjugateInverseWishart* sigma;

	Var<PosReal>* agescale;
	GlobalScalingFunction* scalefunction;
};



#endif	/* PUREBROWNIANPROCESS_H */

