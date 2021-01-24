
#ifndef MEANJITTER_H
#define MEANJITTER_H

#include "BranchProcess.h"
#include "MultiVariateTreeProcess.h"

class WhiteNoiseBranchMean : public Rvar<Real>	{

	public:

	WhiteNoiseBranchMean(Var<Real>* inmean, Var<PosReal>* invariance, Var<PosReal>* intime)	{

		mean = inmean;
		variance = invariance;
		time = intime;
		Register(mean);
		Register(variance);
		Register(time);
		Sample();
	}

	double logProb()	{
		return -0.5 * log(2 * Pi * variance->val() / time->val()) -0.5 * (this->val() - mean->val()) * (this->val() - mean->val()) / variance->val() * time->val();
	}

	protected:

	void drawSample()	{
		setval((Random::sNormal() + mean->val()) * sqrt(variance->val() / time->val()));
	}

	private:

	Var<Real>* mean;
	Var<PosReal>* variance;
	Var<PosReal>* time;

};

class WhiteNoiseProcess : public BranchProcess<Real>	{

	public:

	WhiteNoiseProcess(LengthTree* intree, Var<Real>* inmean, Var<PosReal>* invar)  : BranchProcess<Real>(intree->GetTree()) {
		SetWithRoot(false);
		mean = inmean;
		var = invar;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~WhiteNoiseProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return tree->GetTree();}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	protected:

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<Real>* CreateBranchVal(const Link* link)	{
		return new WhiteNoiseBranchMean(mean,var,tree->GetBranchVal(link->GetBranch()));
	}

	private:

	Var<Real>* mean;
	Var<PosReal>* var;
	LengthTree* tree;

};

class GammaWhiteNoiseBranchMean : public Rvar<PosReal>	{

	public:

	GammaWhiteNoiseBranchMean(Var<PosReal>* inalpha, Var<PosReal>* intime, Var<PosReal>* inmean = 0, bool inintegral = false, bool invariance = false)	{

		alpha = inalpha;
		mean = inmean;
		time = intime;
		integral = inintegral;
		variance = invariance;
		Register(alpha);
		if (mean)	{
			Register(mean);
		}
		Register(time);
		Sample();
	}

	double logProb()	{
		if (integral)	{
			double m = time->val();
			double v = alpha->val() * time->val();
			if (mean)	{
				m *= mean->val();
				v *= mean->val() * mean->val();
			}
			double a = m*m/v;
			double b = m/v;
			return a * log(b) - Random::logGamma(a) + (a-1) * log(this->val()) - b * this->val();
		}
		double a = 1.0;
		if (variance)	{
			a = 1.0 / alpha->val() * time->val();
		}
		else	{
			a = alpha->val() * time->val();
		}
		return a * log(a) - Random::logGamma(a) + (a-1) * log(this->val()) - a * this->val();
	}

	protected:

	void drawSample()	{
		if (integral)	{
			double m = time->val();
			double v = alpha->val() * time->val();
			/*
			if (mean)	{
				m *= mean->val();
				v *= mean->val() * mean->val();
			}
			*/
			double a = m*m/v;
			double b = m/v;
			double tmp = Random::Gamma(a,b);
			if (tmp < 1e-4)	{
				tmp = 1e-4;
			}
			if (mean)	{
				tmp *= mean->val();
			}
			setval(tmp);
			if (std::isnan(((double) val())))	{
				cerr << "nan in gamma white noise\n";
				cerr << m << '\t' << v << '\t' << a << '\t' << b << '\n';
				cerr << "time : " << time->val() << '\n';
				cerr << "alpha : " << alpha->val() << '\n';
				if (mean)	{
					cerr << "mean : " << mean->val() << '\n';
				}
				exit(1);
			}
		}
		else	{
			double a = 1.0;
			if (variance)	{
				a = 1.0 / alpha->val() * time->val();
			}
			else	{
				a = alpha->val() * time->val();
			}
			setval(Random::Gamma(a,a));
			if (val() < 1e-4)	{
				setval(1e-4);
			}
		}
	}

	private:

	Var<PosReal>* alpha;
	Var<PosReal>* mean;
	Var<PosReal>* time;
	bool integral;
	bool variance;

};

class GammaWhiteNoiseProcess : public BranchProcess<PosReal>	{

	public:

	GammaWhiteNoiseProcess(LengthTree* intree, Var<PosReal>* inalpha, Var<PosReal>* inmean = 0, bool inintegral = false, bool invariance = false)  : BranchProcess<PosReal>(intree->GetTree()) {
		SetWithRoot(false);
		alpha = inalpha;
		mean = inmean;
		tree = intree;
		integral = inintegral;
		variance = invariance;
		RecursiveCreate(GetRoot());
	}

	~GammaWhiteNoiseProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return tree->GetTree();}

	double GetTotal()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total;
	}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	protected:

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			this->GetBranchVal(link->GetBranch())->setval(1.0);
			RecursiveReset(link->Out());
		}
	}

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new GammaWhiteNoiseBranchMean(alpha,tree->GetBranchVal(link->GetBranch()),mean,integral,variance);
	}

	private:

	Var<PosReal>* alpha;
	Var<PosReal>* mean;
	LengthTree* tree;
	bool integral;
	bool variance;
};


class MeanJitteredLogit: public Dvar<UnitReal>{

	private:

		MultiNormal* up;
		MultiNormal* down;
		int index;
		Var<Real>* whitenoisemean;

	public:

	MeanJitteredLogit(MultiNormal* inup, MultiNormal* indown, int inindex, Var<Real>* inwhitenoisemean){
		up =inup;
		down = indown;
		index = inindex;
		Register(up);
		Register(down);
		if (whitenoisemean)	{
			Register(whitenoisemean);
		}
		specialUpdate();
	}

	void specialUpdate(){
		double x = 0.5 * ((*up)[index] + (*down)[index]);
		if (whitenoisemean)	{
			x += whitenoisemean->val();
		}
		double y = exp(x) / (1.0 + exp(x));
		setval(y);
	}
};

class MeanJitteredLogitTree: public BranchValPtrTree<Dvar<UnitReal> >	{

	public:

	MeanJitteredLogitTree(MultiVariateTreeProcess* inprocess, int inindex, BranchProcess<Real>* inwnprocess, bool withroot = true)	{
		process = inprocess;
		wnprocess = inwnprocess;
		index = inindex;
		SetWithRoot(withroot);
		RecursiveCreate(GetRoot());
	}

	~MeanJitteredLogitTree()	{
		RecursiveDelete(GetRoot());
	}

	Var<UnitReal>* GetRootRate() {return GetBranchVal(0);}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	double GetVar()	{
		int n =0;
		double tmp1 = RecursiveGetMeanSquare(GetRoot(),n);
		int m = 0;
		double tmp2 = RecursiveGetMean(GetRoot(),m);
		tmp1 /= n;
		tmp2 /= n;
		tmp1 -= tmp2 * tmp2;
		return tmp1;
	}

	double GetMean()	{
		int n =0;
		double tmp = RecursiveGetMean(GetRoot(),n);
		return tmp / n;
	}

	double GetTotal()	{
		int n =0;
		double tmp = RecursiveGetMean(GetRoot(),n);
		return tmp;
	}

	protected:

	double RecursiveGetMean(Link* from, int& n)	{
		double tmp = 0;
		if ((! from->isRoot()) || WithRoot())	{
			tmp += GetBranchVal(from->GetBranch())->val();
			n++;
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += RecursiveGetMean(link->Out(),n);
		}
		return tmp;
	}

	double RecursiveGetMeanSquare(Link* from, int& n)	{
		double tmp = 0;
		if ((! from->isRoot()) || WithRoot())	{
			double t = GetBranchVal(from->GetBranch())->val();
			tmp += t * t;
			n++;
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			tmp += RecursiveGetMeanSquare(link->Out(),n);
		}
		return tmp;
	}

	void specialUpdate(Link* from)	{
		if ((! from->isRoot()) || WithRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<UnitReal>* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new MeanJitteredLogit(process->GetMultiNormal(link), process->GetMultiNormal(link), index, 0);
		}
		return new MeanJitteredLogit(process->GetMultiNormal(link), process->GetMultiNormal(link->Out()), index, wnprocess->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return process->GetTree();}

	private:

	MultiVariateTreeProcess* process;
	BranchProcess<Real>* wnprocess;
	int index;
};

#endif
