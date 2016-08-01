
#ifndef COMMUTATIVEBROWNIANPROCESS_H
#define	COMMUTATIVEBROWNIANPROCESS_H

#include "BrownianProcess.h"
#include "ValTree.h"

class CommutativeIntegralPath : public Dvar<PosReal> {
	public:

		CommutativeIntegralPath(RandomBrownianPath *inbrownianPath, Var<RealVector> *inup, Var<RealVector> *indown, int inindex = 0) {

			brownianPath = inbrownianPath;
			up = inup;
			down = indown;
			index = inindex;
			SetName("Commutative Brownian Path Average");

			Register(up);
			Register(down);
			Register(brownianPath);
			Register(brownianPath->getUp());
			Register(brownianPath->getDown());


			specialUpdate();
		}
		void specialUpdate() {
			 setval(brownianPath->getIntegral(index, up->val()[index], down->val()[index]));
			 // setval(brownianPath->getIntegral(0, up->val()[0], down->val()[0]));
		}
	private:

		RandomBrownianPath *brownianPath;
		Var<RealVector> *up, *down;
		int index;
};

class CommutativeBrownianProcess : public BranchValPtrTree<Dvar<PosReal> >
{
	public:

		CommutativeBrownianProcess(BrownianProcess* inbrownianProcess, int inindex = 0) {
			SetWithRoot(false);
			brownianProcess = inbrownianProcess;
			index = inindex;
			RecursiveCreate(GetRoot());
		}
		virtual ~CommutativeBrownianProcess() {
			RecursiveDelete(GetRoot());
		}

		virtual Dvar<PosReal>* CreateBranchVal(const Link* link) {
			// grep the instant value associated to the node immediately upstream
			MultiNormal* vdown = brownianProcess->GetInstantProcess()->GetMultiNormal(link->Out());

			// grep the instant value associated to this node
			MultiNormal* vup = brownianProcess->GetInstantProcess()->GetMultiNormal(link);

			// grep the brownian path of the branch
			RandomBrownianPath* brownianPath = brownianProcess->GetPureBrownianProcess()->GetRandomBrownianPath(link);

			// make the new transition matrix, and return the pointer
			return new CommutativeIntegralPath(brownianPath, vup, vdown, index);
		}

		void specialUpdate() {
			RecursiveSpecialUpdate(GetRoot());
		}

		void RecursiveSpecialUpdate(const Link* from) {
			for(Link* link=from->Next(); link!=from; link=link->Next()) {
				RecursiveSpecialUpdate(link->Out());
				GetBranchVal(link->GetBranch())->specialUpdate();
			}
		}

		// various accessors
		virtual Tree* GetTree() {
			return brownianProcess->GetTree();
		}

	protected:
	private:

		BrownianProcess* brownianProcess;
		int index;

};



#endif	/* COMMUTATIVEBROWNIANPROCESS_H */

