
#include "CalibratedChronogram.h"
#include <list>


class ImproperUniform: public Rvar<PosReal>	{

	public:

	ImproperUniform(DAGnode* inroot = 0)	{
		root = inroot;
		if (root)	{
			Register(root);
		}
		Sample();
	}

	double logProb()	{
		return 0;
	}

	void drawSample()	{
		setval(1);
	}

	protected:

	DAGnode* root;
};


class DivBDCalibratedNodeDate : public CalibratedNodeDate	{

	private:

	Var<PosReal>* lambda;
	Var<PosReal>* mu;
	int isroot;

	public:

	DivBDCalibratedNodeDate(Var<PosReal>* inRate, Var<PosReal>* inScale, Var<PosReal>* inlambda, Var<PosReal>* inmu, int inisroot) : CalibratedNodeDate(inRate,inScale), lambda(inlambda), mu(inmu), isroot(inisroot)	{
		Register(lambda);
		Register(mu);
	}

	protected:


	double logg(double t)	{
		double ret = 0;
		if (mu->val() < 1e-10)	{
			ret = -lambda->val()*t - 2*log(lambda->val());
		}
		else if (lambda->val() > mu->val())	{
			double d = (lambda->val() - mu->val()) * t;
			double e = 0;
			if (d < 250)	{
				e = exp(-d);
			}
			ret = - d - 2 * log(lambda->val() - mu->val()*e);
		}
		else	{
			double d = (mu->val() - lambda->val()) * t;
			double e = 0;
			if (d < 250)	{
				e = exp(-d);
			}
			ret = - d - 2 * log(mu->val() - lambda->val()*e);
		}

		if (std::isnan(ret))	{
			cerr << "internal log prob is nan\n";
			exit(1);
		}

		return ret;
	}

	double internalLogProb()	{
		double t = val() * scale->val();
		// double t = val() * scale->val() / 100;
		return logg(t) + log(scale->val());
	}

	double rootLogProb()	{
		int ntaxa = isroot;
		double t = val() * scale->val();
		// double t = val() * scale->val() / 100;

		double ret = logg(t);
		if (mu->val() < 1e-10)	{
			ret += (3*ntaxa-4)*log(lambda->val()) - lambda->val() * t;
		}
		else if (lambda->val() > mu->val())	{
			double d = (lambda->val() - mu->val()) * t;
			double e = 0;
			if (d < 250)	{
				e = exp(-d);
			}
			ret += (ntaxa-1) * log(lambda->val()) - log(mu->val()) + (2*ntaxa-2)*log(lambda->val() - mu->val()) + log(log(lambda->val() / (lambda->val() - mu->val()*e)));
		}
		else	{
			double d = (mu->val() - lambda->val()) * t;
			double e = 0;
			if (d < 250)	{
				e = exp(-d);
			}
			ret += (ntaxa-1) * log(lambda->val()) - log(mu->val()) + (2*ntaxa-2)*log(mu->val() - lambda->val()) + log(log(lambda->val() * e / (mu->val() - lambda->val()*e)));
		}

		if (std::isnan(ret))	{
			cerr << "root log prob is nan\n";
			cerr << "lambda : " << lambda->val() << '\n';
			cerr << "mu     : " << mu->val() << '\n';
			exit(1);
		}
		return ret;
	}

	double logProb()	{
		if (isroot == -1)	{
			return 0;
		}
		double total = 0;
		// here add log prob from birth death
		if (isroot)	{
			total += rootLogProb();
		}
		else	{
			total += internalLogProb();
		}
		if ((olderlimit != -1) && (value > olderlimit/scale->val()))	{
			// cerr << "too old\n";
			total -= 10000;
		}
		if ((youngerlimit != -1) && (value < youngerlimit/scale->val()))	{
			// cerr << "too young\n";
			total -= 10000;
		}
		return total;
	}



};

class DivBDCalibratedChronogram : public CalibratedChronogram	{


	public:

	DivBDCalibratedChronogram(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* inlambda, Var<PosReal>* inmu, CalibrationSet* incalibset)	{

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		lambda = inlambda;
		mu = inmu;

		scale = new ImproperUniform(rate);
		scale->setval(100);

		calibset = incalibset;

		Ntaxa = GetTree()->GetSize(GetRoot());
		cerr << "Ntaxa : " << Ntaxa << '\n';

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		RecursiveSetCalibrations(GetRoot());

		RecursiveYoungerLimit(GetRoot());
		RecursiveOlderLimit(GetRoot(),-1);
		double d = RecursiveDrawAges(GetRoot());

		cerr << "age  " << d << '\n';

		RecursiveNormalizeTree(GetRoot(),d,false);
		RecursiveUpdateBranches(GetRoot());
		scale->setval(d);
		map<DAGnode*,int> tmpmap;
		scale->FullCorrupt(tmpmap);
		scale->FullUpdate();

		double tmp = RecursiveGetLogProb(GetRoot());
		cerr << tmp << '\n';
		if (std::isnan(tmp))	{
			cerr << "log prob is nan\n";
			exit(1);
		}
		/*
		if (fabs(tmp) > 1e-6)	{
			cerr << "error : nodes out of calib\n";
			exit(1);
		}
		*/
	}

	~DivBDCalibratedChronogram()	{
		RecursiveDeleteBranch(GetRoot());
		RecursiveDeleteNode(GetRoot());
		delete scale;
	}

	double GetYuleLogProb()	{
		return (Ntaxa - 2) * log(lambda->val()) - lambda->val() * scale->val() * (GetTotalTime()) + (Ntaxa - 2) * log(scale->val());
	}

	void RecursiveGetAgeList(const Link* from, list<double>& agelist)	{
		if (! from->isLeaf())	{
			double temp = GetNodeVal(from->GetNode())->val();
			agelist.push_back(temp);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetAgeList(link->Out(),agelist);
		}
	}


	void GetAgeList(list<double>& agelist)	{

		RecursiveGetAgeList(GetRoot(),agelist);
		agelist.sort();
	}

	protected:

	Rvar<PosReal>* CreateNodeVal (const Link* link){
		if (link->isRoot())	{
			return new DivBDCalibratedNodeDate(rate,scale,lambda,mu,Ntaxa);
		}
		else if (link->isLeaf())	{
			return new DivBDCalibratedNodeDate(rate,scale,lambda,mu,-1);
		}
		return new DivBDCalibratedNodeDate(rate,scale,lambda,mu,0);
	}

	private:

	unsigned int Ntaxa;
	Var<PosReal>* lambda;
	Var<PosReal>* mu;

};

