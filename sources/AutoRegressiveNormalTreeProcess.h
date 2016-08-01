#ifndef AUTOREGRESSIVENORMALTREEPROCESS_H
#define AUTOREGRESSIVENORMALTREEPROCESS_H

#include "PrecisionNormalTreeProcess.h"


class  AutoRegressiveNormalProcessInstantValue : public NormalProcessInstantValue	{

	public:

	AutoRegressiveNormalProcessInstantValue(Var<PosReal>* intau, Var<PosReal>* inphi, Var<Real>* inup = 0, Var<PosReal>* intime = 0, bool inprec = true)	{
		tau = intau;
		prec = inprec;
		up = inup;
		time = intime;
		phi = inphi;

		if (up)	{
			Register(up);
		}
		Register(phi);
		if (time)	{
			Register(time);
		}
		Register(tau);
		if (up) {
			SetName("normal process instant value\n");
		}
		else	{
			SetName("root normal process instant value\n");
		}
		Sample();
	}

	virtual double GetFiniteTimeTransitionLogProb()	{

		double mean = 0;
		double sigma = 0;
		if (prec)	{
			sigma = 1.0 / tau->val();
			cerr << "precision?\n";
			exit(1);
		}
		else	{
			sigma = tau->val();
		}
		double expo = exp(-phi->val() * time->val());
		double mu = sqrt ((1 - exp(-2 * phi->val() * time->val())) / 2 / phi->val());
		double d = (val() - up->val() * expo - mean * (1-expo)) / mu;
		double total = -0.5 * (log(sigma) + d*d/sigma) - log(mu);
		return total;

		/*
		double tXSX = 0;
		double* dval = new double[GetDim()];
		for (int i=0; i<GetDim() ; i++) {
			dval[i] =  (val()[i] - up->val()[i] * expo - (*mean)[i] * (1 - expo))/mu;
		}
		for (int i=0; i<GetDim() ; i++) {
			tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
			for (int j=0; j<i ; j++) {
				tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
			}
		}
		double d = -0.5 * (sigma->GetLogDeterminant() + tXSX) -  GetDim() * log(mu);
		delete[] dval;
		return d;
		*/
	}

	// normal process has no stationary distribution!
	// let us define a default: standard normal distribution
	virtual double GetStationaryLogProb()	{
		double mean = 0;
		double mu = 1.0 / sqrt(2 * phi->val());
		double d = (val() - mean) / mu;
		double sigma = 0;
		if (prec)	{
			sigma = 1.0 / tau->val();
		}
		else	{
			sigma = tau->val();
		}
		double tot = -0.5 * (log(sigma) + d * d / sigma) - log(mu);
		return tot;

		/*
		for (int i=0; i<GetDim() ; i++) {
			dval[i] =  (val()[i] - (*mean)[i]) / mu;
		}
		for (int i=0; i<GetDim() ; i++) {
			tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
			for (int j=0; j<i ; j++) {
				tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
			}
		}
		double d = -0.5 * (sigma->GetLogDeterminant() + tXSX) -  GetDim() * log(mu);
		return d;
		*/

	}

	virtual void DrawFiniteTime()	{
		if (prec)	{
			setval(Random::sNormal() * sqrt(time->val() / tau->val()) + up->val());
		}
		else	{
			setval(Random::sNormal() * sqrt(time->val() * tau->val()) + up->val());
		}
	}

	virtual void DrawStationary()	{
		setval(0);
	}

	protected:

	Var<PosReal>* phi;
};


class AutoRegressiveLogNormalTreeProcess : public LogNormalTreeProcess	{

	public:


	AutoRegressiveLogNormalTreeProcess() {}

	AutoRegressiveLogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma, Var<PosReal>* inphi, BranchValType inbval, bool inprec = true)  {
		SetWithRoot(false);
		bval = inbval;
		tree = intree;
		sigma = insigma;
		phi = inphi;
		prec = inprec;
		cerr << "recursive create\n";
		RecursiveCreate(GetRoot());
		cerr << "ok\n";
		// GetNodeVal(GetRoot()->GetNode())->ClampAt(0);
	}

	~AutoRegressiveLogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = 0;
		if (! link->isRoot())	{
			time = GetLengthTree()->GetBranchLength(link->GetBranch());
		}

		// make the new instant value, and return the pointer
		return new AutoRegressiveNormalProcessInstantValue(sigma, phi, vup, time, prec);
	}

	void GetMeanAndVar(double& mean, double& var)	{
		mean = 0;
		var = 0;
		int count = 0;
		RecursiveGetMeanAndVar(GetRoot(),mean,var,count);
		mean/=count;
		var/=count;
		var-=mean*mean;
	}

	protected:

	void RecursiveGetMeanAndVar(const Link* from, double& mean, double& var, int& count)	{
		if (! from->isRoot())	{
			double tmp = GetBranchVal(from->GetBranch())->val();
			mean += tmp;
			var += tmp*tmp;
			count++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetMeanAndVar(link->Out(),mean,var,count);
		}
	}

	Var<PosReal>* phi;
};

class FlexRhoAutoRegressiveLogNormalTreeProcess : public AutoRegressiveLogNormalTreeProcess	{

	public:


	FlexRhoAutoRegressiveLogNormalTreeProcess(LengthTree* intree, Var<PosReal>* insigma, LengthTree* inphiprocess, Var<PosReal>* inphi, BranchValType inbval, bool inprec = true)  {
		SetWithRoot(false);
		bval = inbval;
		tree = intree;
		sigma = insigma;
		phi = inphi;
		phiprocess = inphiprocess;
		prec = inprec;
		cerr << "recursive create\n";
		RecursiveCreate(GetRoot());
		cerr << "ok\n";
		// GetNodeVal(GetRoot()->GetNode())->ClampAt(0);
	}

	~FlexRhoAutoRegressiveLogNormalTreeProcess()	{
		RecursiveDelete(GetRoot());
	}

	virtual Rvar<Real>* CreateNodeVal(const Link* link)	{

		// grep the instant value associated to the node immediately upstream
		Var<Real>* vup = GetNodeVal(link->Out()->GetNode());

		// grep the time associated to the branch
		Var<PosReal>* time = 0;
		if (! link->isRoot())	{
			time = GetLengthTree()->GetBranchLength(link->GetBranch());
			return new AutoRegressiveNormalProcessInstantValue(sigma, phiprocess->GetBranchVal(link->GetBranch()), vup, time, prec);
		}

		// make the new instant value, and return the pointer
		return new AutoRegressiveNormalProcessInstantValue(sigma, phi, vup, time, prec);
	}

	protected:

	Var<PosReal>* phi;
	LengthTree* phiprocess;
};


#endif
