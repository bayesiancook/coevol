
#ifndef CONJUGATELOGNORMALTREEPROCESS_H
#define CONJUGATELOGNORMALTREEPROCESS_H

#include "PrecisionNormalTreeProcess.h"
#include "ConjugatePoisson.h"

class ConjugateNormalProcessInstantValue : public virtual ConjugateSampling<Real>, public virtual NormalProcessInstantValue   {

	public:

	// up is the instant value at time 0
	// "this" is the instant value at time "time"
	// both are Var<Real> (normal process)
	//
	ConjugateNormalProcessInstantValue(Var<PosReal>* intau, Var<Real>* inup = 0, Var<PosReal>* intime = 0) : NormalProcessInstantValue(intau,inup,intime) {
		conj_tau = dynamic_cast<PoissonSemiConjugate*>(intau);
		if (up && conj_tau)	{
			conjugate_up.insert(conj_tau);
		}
	}

	void AddSufficientStatistic(SemiConjPrior* parent) {
		if (up && conj_tau && (parent == conj_tau))	{
			conj_tau->AddToShape(0.5);
			conj_tau->AddToScale(0.5 * IndependentContrast());
		}
	}

	double IndependentContrast()	{
		if (! up)	{
			cerr << "error in ConjugateNormalProcess : independent contrast called on root\n";
			throw;
		}
		if (fabs(time->val()) < 1e-7)	{
			cerr << "error in ConjugateNormalProcess : time too small : " << time->val() << '\n';
			throw;
		}
		double tmp = val() - up->val();
		return tmp * tmp / time->val();
	}

	protected:

	PoissonSemiConjugate* conj_tau;
};


class ConjugateLogNormalTreeProcess : public LogNormalTreeProcess	{

	public:

	ConjugateLogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma)  {
		tree = intree;
		sigma = insigma;
		RecursiveCreate(GetRoot());
		GetNodeVal(GetRoot()->GetNode())->ClampAt(0);
	}

	~ConjugateLogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

		// make the new instant value, and return the pointer
		return new ConjugateNormalProcessInstantValue(sigma, vup, time);
	}

};


#endif
