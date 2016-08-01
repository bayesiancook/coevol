
#ifndef ALPHASTABLE_H
#define ALPHASTABLE_H

#include "BranchProcess.h"
#include "stabledistribution.h"

class RandomAlphaStablePath : public Rvar<RealVector>	{

	public:

	RandomAlphaStablePath(Var<PosReal>* inlength, RandomAlphaStablePath* inup, Var<PosReal>* inalpha, Var<PosReal>* inbeta, int innTreeSegments, Var<PosReal>* inrho = 0, Var<PosReal>* insigma = 0) {
		
		nTreeSegments = innTreeSegments;

		length = inlength;
		up = inup;

		alpha = inalpha;
		beta = inbeta;
		rho = inrho;
		sigma = insigma;

		int nseg = 0;
		if (length)	{
			nseg = (int) (nTreeSegments*length->val());
			if(nseg<3) {
				nseg =3;
			}
		}

		setval(RealVector(nseg+1));
		bkvalue = RealVector(nseg+1);

		Register(length);
		Register(up);
		Register(alpha);
		Register(beta);
		Register(rho);
		Register(sigma);

		Sample();
	}

	void drawSample()	{

		if (! GetNSegments())	{
			(*this)[0] = 0;
		}
	
		else	{
			(*this)[0] = up->GetFinalValue();
			double segmentlength = GetSegmentLength();
			for (int i=0; i<GetNSegments(); i++)	{
				(*this)[i+1] = (*this)[i] + DrawFiniteTime(segmentlength);
			}
		}
	}

	double logProb()	{

		double logprob = 0;
		if (alpha->val() < 0.2)	{
			return log(0);
		}
		if (rho && (alpha->val() != 2.0))	{
			cerr << "error in alpha stable process: poisson jumps only with brownian motion\n";
			exit(1);
		}
		double segmentlength = GetSegmentLength();
		for (int i=0; i<GetNSegments(); i++)	{
			logprob += FiniteTimeLogProb((*this)[i+1] - (*this)[i], segmentlength);

		}
		return logprob;
	}
	 
	double GetInitValue()	{
		return (*this)[0];
	}

	double GetFinalValue()	{
		return (*this)[GetNSegments()];
	}

	int GetNSegments()	{
		return GetDim() - 1;
	}

	double GetSegmentLength() {
		if (! length)	{
			return 0;
		}
		return length->val() / GetNSegments();
	}

	void Translate(double m)	{

		int n = GetNSegments();
		for (int i=0; i<=n; i++)	{
			(*this)[i] += m;
		}
	}

	void TranslateFromBase(double m)	{

		int n = GetNSegments();
		if (n)	{
			for (int i=0; i<=n; i++)	{
				(*this)[i] += m * ((double) (n-i)) / n;
			}
		}
		else	{
			(*this)[0] += m;
		}
	}

	void TranslateFromTip(double m)	{

		int n = GetNSegments();
		if (n)	{
			for (int i=0; i<=n; i++)	{
				(*this)[i] += m * ((double) i) / n;
			}
		}
		else	{
			(*this)[0] += m;
		}
	}

	double DrawFiniteTime(double time)	{

		double ret = 0;
			double scale = beta->val() * sqrt(time);
			ret = scale * Random::sNormal();
		/*
		// Cauchy
		if (alpha->val() == 1.0)	{
			double u = Random::Uniform() - 0.5;
			double scale = beta->val();
			ret = scale * tan(Pi * u);
		}

		// Brownian
		else if (alpha->val() == 2.0)	{
			double scale = beta->val() * sqrt(time);
			ret = scale * Random::sNormal();
		}

		else	{
			cerr << "error: only cauchy or gaussian\n";
			exit(1);
		}
		*/

		return ret;
	}

	double FiniteTimeLogProb(double delta, double time)	{

		if (bkalpha != alpha->val())	{
			bkalpha = alpha->val();
			as.setAlpha(alpha->val());
		}

		double logprob = 0;
		if (rho)	{

			double varb = beta->val() * beta->val() * time;
			double s2 = sigma->val() * sigma->val();

			double scale = rho->val() * time;
			double expo = exp(-scale);
			double totweight = 0;
			double prob = 0;
			double fact = 1;
			double num = 1;
			int i = 0;
			double poissoncutoff = 0.99;
			int max = 20;
			while ((i<max) && (totweight < poissoncutoff))	{
				double w = expo * num / fact;
				totweight += w;
				double var = i * s2 + varb;
				prob += w * exp(-0.5 * delta * delta / var) / sqrt(var);
				i++;
				fact *= i;
				num *= scale;
			}
			logprob = log(prob/totweight);
			/*
			if (i == max)	{
				cerr << "max : " << totweight << '\t' << rho->val() << '\t' << time << '\n';
			}
			*/
		}
		else	{
			double scale = beta->val() * exp(1.0/alpha->val() * log(time));
			double y = delta / scale;
			logprob  = -log(scale) + as.logPDF(y,1.0);
		}
		/*
		double logprob2 = as.logPDF(delta,scale);
		double error = fabs(logprob - logprob2);
		if (error > 1e-5)	{
			cerr << "error in finite time log prob\n";
			cerr << logprob << '\t' << logprob2 << '\n';
		}
		*/
		return logprob;

		/*
		double logprob = 0;
		// Cauchy
		if (alpha->val() == 1.0)	{
			double scale = beta->val() * time;
			double y = delta/scale;
			logprob = -log(scale) - log(1 + y*y);
			if (max < fabs(y))	{
				max = fabs(y);
			}
		}

		// Brownian
		else if (alpha->val() == 2.0)	{
			double scale = beta->val() * sqrt(time);
			double y = delta/scale;
			logprob = -log(scale) - 0.5 * y*y;
			if (max < fabs(y))	{
				max = fabs(y);
			}
		}

		// other cases
		else	{
			cerr << "error: only cauchy or gaussian\n";
			exit(1);
			double scale = beta->val() * exp(log(time)/alpha->val());
			logprob = -log(scale) + GetStandardlnPDF(delta/scale,alpha->val());
		}

		return logprob;
		*/
	}

	double logNormalDensity(double scale, double x)	{
		return -0.5 * log(2 * Pi) - log(scale) - 0.5 * x * x / scale / scale;
	}

	double logCauchyDensity(double scale, double x)	{
		return -log (Pi * scale * (1 + x*x/scale/scale));
	}

	double CauchyDensity(double scale, double x)	{
		return 1.0/ (Pi * scale * (1 + x*x/scale/scale));
	}

	double DrawConditionalCauchy(double x0, double x1, double t0, double t1, double tuning, int n, int p, double* grid, double* cumul)	{

		double scale = beta->val() * tuning;
		if ((! grid) || (! cumul))	{
			cerr << "error: grid cumul not allocated\n";
			exit(1);
		}
		
		double total = 0;
		for (int i=0; i<=2*n; i++)	{
			double xgrid = tuning * p * (i-n) / n;
			double tmp = CauchyDensity(scale*t0,xgrid-x0) * CauchyDensity(scale*t1,xgrid-x1);
			total += tmp;
			grid[i] = tmp;
			cumul[i] = total;
		}

		double u = total * Random::Uniform();
		int k = 0;
		while ((k<=2*n) && (u>cumul[k]))	{
			k++;
		}
		if (k > 2*n)	{
			cerr << "error in draw cond cauchy: overflow\n";
			exit(1);
		}

		double v = Random::Uniform() - 0.5;
		double x = tuning * p * (k - n + v) / n;

		return x;
	}

	double DrawConditionalNormal(double x0, double x1, double t0, double t1, double tuning)	{

		double tau0 = 1.0 / beta->val() / beta->val() / t0;
		double tau1 = 1.0 / beta->val() / beta->val() / t1;

		double m = (x0*tau0 + x1*tau1) / (tau0 + tau1);

		double scale = tuning / sqrt(tau0 + tau1);

		double y = scale*Random::sNormal();
		double x = m + y;
		if (isnan(x))	{
			cerr << "draw cond normal: nan \n";
			cerr << x0 << '\t' << x1 << '\t' << t0 << '\t' << t1 << '\t' << tuning << '\t' << tau0 << '\t' << tau1 << '\t' << m << '\t' << scale << '\n';
			exit(1);
		}
		return x;
	}

	// resample conditional on end points
	double ResampleConditional(double tuning, int n=100, int p=10)	{

		int nsegments = GetNSegments();

		double* y = new double[nsegments+1];
		y[0] = 0;
		y[nsegments] = 0;

		// Cauchy
		if (alpha->val() == 1.0)	{

			double* grid = new double[2*n+1];
			double* cumul = new double[2*n+1];


			for (int i=1; i<nsegments; i++)	{
				y[i] = DrawConditionalCauchy(y[i-1], y[nsegments], length->val()/nsegments, length->val()*(nsegments-i)/nsegments, tuning, n, p, grid, cumul);
				(*this)[i] += y[i];
			}
			
			delete[] grid;
			delete[] cumul;
		}

		else {

			for (int i=1; i<nsegments; i++)	{
				y[i] = DrawConditionalNormal(y[i-1], y[nsegments], length->val()/nsegments, length->val()*(nsegments-1)/nsegments, tuning);
				(*this)[i] += y[i];
			}
		}

		delete[] y;
		return 0;
	}

	virtual double ProposeMove(double tuning)	{

		ResampleConditional(tuning);
		return 0;
	}

	double GetExpRelVar()	{
		if (length)	{
			double tmp = GetExpMean();
			return GetExpMeanSquare() / tmp * tmp;
		}
		return 1;
	}

	double GetExpMean()	{

		if (length)	{
			return GetExpIntegral() / length->val();
		}
		return 0;
	}

	double GetExpIntegral()	{

		double segmentlength = GetSegmentLength();
		double integral = 0;
		for (int i=0; i<GetNSegments(); i++)	{
			integral += 0.5 * (exp((*this)[i]) + exp((*this)[i+1])) * segmentlength;
		}
		if (isinf(integral))	{
			cerr << "get exp int: inf\n";
			cerr << "length : " << length->val() << '\n';
			cerr << "nseg : " << GetNSegments() << '\n';
			for (int i=0; i<=GetNSegments(); i++)	{
				cerr << (*this)[i] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
		if (isnan(integral))	{
			cerr << "get exp int: nan\n";
			cerr << "length : " << length->val() << '\n';
			cerr << "nseg : " << GetNSegments() << '\n';
			for (int i=0; i<=GetNSegments(); i++)	{
				cerr << (*this)[i] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
		return integral;
	}

	double GetExpMeanSquare()	{

		if (! length)	{
			double tmp = exp((*this)[0]);
			return tmp*tmp;
		}

		double integral = 0;
		for (int i=0; i<GetNSegments(); i++)	{
			double tmp = 0.5 * (exp((*this)[i]) + exp((*this)[i+1]));
			integral += tmp*tmp;
		}
		integral /= GetNSegments();

		if (isinf(integral))	{
			cerr << "get exp int: inf\n";
			cerr << "length : " << length->val() << '\n';
			cerr << "nseg : " << GetNSegments() << '\n';
			for (int i=0; i<=GetNSegments(); i++)	{
				cerr << (*this)[i] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
		if (isnan(integral))	{
			cerr << "get exp int: nan\n";
			cerr << "length : " << length->val() << '\n';
			cerr << "nseg : " << GetNSegments() << '\n';
			for (int i=0; i<=GetNSegments(); i++)	{
				cerr << (*this)[i] << '\t';
			}
			cerr << '\n';
			exit(1);
		}
		return integral;
	}

	int RescaleBridge(double factor)	{

		int nseg = GetNSegments();
		if (nseg)	{
			for (int i=1; i<nseg; i++)	{
				double x = (*this)[0] + ((*this)[nseg] - (*this)[0]) * i / nseg;
				(*this)[i] = x + factor * ((*this)[i] - x);
			}
			return nseg-1;
		}
		else	{
			return 0;
		}
	}

	double AlphaRescaleBridge(double alpha1, double alpha2)	{

		double logh = 0;
		int nseg = GetNSegments();
		if (nseg)	{
			for (int i=1; i<nseg; i++)	{
				double x = (*this)[0] + ((*this)[nseg] - (*this)[0]) * i / nseg;
				double dev = fabs((*this)[i] - x);
				double sign = 1.0;
				if ((*this)[i] < x)	{
					sign = -1.0;
				}
				(*this)[i] = x + sign * exp( alpha2/alpha1 * log(dev));
				logh += (alpha1 / alpha2 - 1) * log(dev);
			}
			logh += (nseg - 1) * (log(alpha1) - log(alpha2));
		}
		return logh;
	}

	int Rescale(double mean, double factor)	{

		for (int i=0; i<=GetNSegments(); i++)	{
			(*this)[i] = mean + factor * ((*this)[i] - mean);
		}

		// should count 1 for root path
		int n = GetNSegments();
		if (!n)	{
			n = 1;
		}
		return n;
	}

	double AlphaRescale(double mean, double alpha1, double alpha2)	{

		double logh = 0;
		int nseg = GetNSegments();
		if (nseg)	{
			for (int i=1; i<nseg; i++)	{
				double dev = fabs((*this)[i] - mean);
				double sign = 1.0;
				if ((*this)[i] < mean)	{
					sign = -1.0;
				}
				(*this)[i] = mean + sign * exp( alpha2/alpha1 * log(dev));
				logh += (alpha1 / alpha2 - 1) * log(dev);
			}
			logh += (nseg - 1) * (log(alpha1) - log(alpha2));
		}
		return logh;
	}

	void GetTotal(double& total, int& n)	{

		int nseg = GetNSegments();
		if (! nseg)	{
			total += (*this)[0];
			n++;
		}
		else	{
			for (int i=1; i<=nseg; i++)	{
				total += (*this)[i];
			}
			n += nseg;
		}
	}

	protected:

	Var<PosReal>* length;
	RandomAlphaStablePath* up;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	Var<PosReal>* rho;
	Var<PosReal>* sigma;

	int nTreeSegments;

	static StableDistribution as;
	static double bkalpha;
};

StableDistribution RandomAlphaStablePath::as(1.0,1.0);
double RandomAlphaStablePath::bkalpha = 1.0;


class AlphaStableProcess : public BranchProcess<RealVector>	{
// class AlphaStableProcess : public MCMC, public BranchValPtrTree< Rvar<RealVector> >	{

	public:

	AlphaStableProcess(LengthTree* inlengthtree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, int innTreeSegments, Var<PosReal>* inrho = 0, Var<PosReal>* insigma = 0) : BranchProcess<RealVector>(inlengthtree->GetTree(), true) {

		nTreeSegments = innTreeSegments;
		lengthtree = inlengthtree;
		alpha = inalpha;
		beta = inbeta;
		rho = inrho;
		sigma = insigma;

		RecursiveCreate(GetRoot(),GetRoot());
	}

	~AlphaStableProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	RandomAlphaStablePath* GetPath(const Branch* branch)	{
		RandomAlphaStablePath* tmp = dynamic_cast<RandomAlphaStablePath*>(GetBranchVal(branch));
		if (GetBranchVal(branch) && (! tmp))	{
			cerr << "error in alpha stable process: null path\n";
			cerr << tmp << '\t' << GetBranchVal(branch) << '\t' << branch << '\n';
			exit(1);
		}
		return tmp;
	}

	void Translate(double m)	{

		RecursiveTranslate(GetRoot(),m);
	}
		
	int RescaleBridge(double factor)	{

		int n = RecursiveRescaleBridge(GetRoot(),factor);
		return n;
	}
		
	int Rescale(double mean, double factor)	{

		int n = RecursiveRescale(GetRoot(),mean,factor);
		return n;
	}
		
	int AlphaRescaleBridge(double alpha1, double alpha2)	{

		int n = RecursiveAlphaRescaleBridge(GetRoot(),alpha1,alpha2);
		return n;
	}
		
	int AlphaRescale(double mean, double alpha1, double alpha2)	{

		int n = RecursiveAlphaRescale(GetRoot(),mean,alpha1,alpha2);
		return n;
	}
		
	double GetMean()	{
		double total = 0;
		int n = 0;
		RecursiveGetTotal(GetRoot(),total,n);
		return total/n;
	}

	protected:

	void RecursiveGetTotal(const Link* from, double& total, int& n)	{

		GetPath(from->GetBranch())->GetTotal(total,n);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetTotal(link->Out(),total,n);
		}
	}

	void RecursiveTranslate(const Link* from, double m)	{

		GetPath(from->GetBranch())->Translate(m);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTranslate(link->Out(),m);
		}
	}

	int RecursiveRescaleBridge(const Link* from, double factor)	{

		int n = GetPath(from->GetBranch())->RescaleBridge(factor);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveRescaleBridge(link->Out(),factor);
		}
		return n;
	}

	int RecursiveRescale(const Link* from, double mean, double factor)	{

		int n = GetPath(from->GetBranch())->Rescale(mean,factor);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveRescale(link->Out(),mean,factor);
		}
		return n;
	}

	int RecursiveAlphaRescaleBridge(const Link* from, double alpha1, double alpha2)	{

		int n = GetPath(from->GetBranch())->AlphaRescaleBridge(alpha1,alpha2);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveAlphaRescaleBridge(link->Out(),alpha1,alpha2);
		}
		return n;
	}

	int RecursiveAlphaRescale(const Link* from, double mean, double alpha1, double alpha2)	{

		int n = GetPath(from->GetBranch())->AlphaRescale(mean,alpha1,alpha2);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveAlphaRescale(link->Out(),mean,alpha1,alpha2);
		}
		return n;
	}

	void RecursiveCreate(const Link* from, const Link* prev)	{

		branchval[from->GetBranch()] = CreateBranchVal(from,prev);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out(),from);
		}
	}

	Rvar<RealVector>* CreateBranchVal(const Link* link)	{
		cerr << "error in alpha stable process: in CreateBranchVal(const Link* link)\n";
		exit(1);
	}

	Rvar<RealVector>* CreateBranchVal(const Link* link, const Link* prevlink)	{

		Rvar<RealVector>* tmp = 0;
		if (link->isRoot())	{
			tmp = new RandomAlphaStablePath(0,0,alpha,beta,nTreeSegments,rho,sigma);
		}
		else	{
			tmp = new RandomAlphaStablePath(lengthtree->GetBranchVal(link->GetBranch()),GetPath(prevlink->GetBranch()),alpha,beta,nTreeSegments,rho,sigma);
		}
		return tmp;
	}
	
	LengthTree* lengthtree;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	Var<PosReal>* rho;
	Var<PosReal>* sigma;
	int nTreeSegments;
};



class AlphaStableExpIntegral : public Dvar<PosReal>	{

	public:

	AlphaStableExpIntegral(Var<PosReal>* inlength, RandomAlphaStablePath* inpath, RandomAlphaStablePath* inup)	{

		length = inlength;
		path = inpath;
		up = inup;

		Register(length);
		Register(path);
		Register(up);

		specialUpdate();
	}

	void specialUpdate()	{
		setval(path->GetExpIntegral());
	}

	protected:

	Var<PosReal>* length;
	RandomAlphaStablePath* path;
	RandomAlphaStablePath* up;
};

class AlphaStableExpIntegralTree : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	AlphaStableExpIntegralTree(LengthTree* inlengthtree, AlphaStableProcess* inprocess)	{
		lengthtree = inlengthtree;
		process = inprocess;
		SetWithRoot(false);
		RecursiveCreate(GetRoot(),GetRoot());
	}

	~AlphaStableExpIntegralTree()	{
		RecursiveDelete(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	protected:

	void specialUpdate(Link* from)	{
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	void RecursiveCreate(const Link* from, const Link* prev)	{

		if (! from->isRoot())	{
			branchval[from->GetBranch()] = CreateBranchVal(from,prev);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out(),from);
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		cerr << "error in alpha stable exp integral: in CreateBranchVal(const Link* link)\n";
		exit(1);
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link, const Link* prevlink)	{

		return new AlphaStableExpIntegral(lengthtree->GetBranchVal(link->GetBranch()),process->GetPath(link->GetBranch()),process->GetPath(prevlink->GetBranch()));
	}

	LengthTree* lengthtree;
	AlphaStableProcess* process;
};

class AlphaStableExpRelVar : public Dvar<PosReal>	{

	public:

	AlphaStableExpRelVar(Var<PosReal>* inlength, RandomAlphaStablePath* inpath, RandomAlphaStablePath* inup)	{

		length = inlength;
		path = inpath;
		up = inup;

		Register(length);
		Register(path);
		Register(up);

		specialUpdate();
	}

	void specialUpdate()	{
		setval(path->GetExpRelVar());
	}

	protected:

	Var<PosReal>* length;
	RandomAlphaStablePath* path;
	RandomAlphaStablePath* up;
};

class AlphaStableExpRelVarTree : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	AlphaStableExpRelVarTree(LengthTree* inlengthtree, AlphaStableProcess* inprocess)	{
		lengthtree = inlengthtree;
		process = inprocess;
		SetWithRoot(false);
		RecursiveCreate(GetRoot(),GetRoot());
	}

	~AlphaStableExpRelVarTree()	{
		RecursiveDelete(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	protected:

	void specialUpdate(Link* from)	{
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	void RecursiveCreate(const Link* from, const Link* prev)	{

		if (! from->isRoot())	{
			branchval[from->GetBranch()] = CreateBranchVal(from,prev);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out(),from);
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		cerr << "error in alpha stable exp integral: in CreateBranchVal(const Link* link)\n";
		exit(1);
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link, const Link* prevlink)	{

		return new AlphaStableExpRelVar(lengthtree->GetBranchVal(link->GetBranch()),process->GetPath(link->GetBranch()),process->GetPath(prevlink->GetBranch()));
	}

	LengthTree* lengthtree;
	AlphaStableProcess* process;
};

class DoubleAlphaStableExpIntegral : public Dvar<PosReal>	{

	public:

	DoubleAlphaStableExpIntegral(Var<PosReal>* inlength, RandomAlphaStablePath* inpath1, RandomAlphaStablePath* inup1, RandomAlphaStablePath* inpath2, RandomAlphaStablePath* inup2)	{

		length = inlength;
		path1 = inpath1;
		up1 = inup1;
		path2 = inpath2;
		up2 = inup2;

		Register(length);
		Register(path1);
		Register(up1);
		Register(path2);
		Register(up2);

		specialUpdate();
	}

	void specialUpdate()	{
		setval(GetExpIntegral());
	}

	double GetExpIntegral()	{

		double segmentlength = path1->GetSegmentLength();
		double current = (*path1)[0] + (*path2)[0];
		double total = 0;
		for (int i=0; i<path1->GetNSegments(); i++)	{
			double next = (*path1)[i+1] + (*path2)[i+1];
			total += 0.5 * (exp(current) + exp(next)) * segmentlength;
			current = next;
		}
		return total;
	}

	protected:

	Var<PosReal>* length;
	RandomAlphaStablePath* path1;
	RandomAlphaStablePath* up1;
	RandomAlphaStablePath* path2;
	RandomAlphaStablePath* up2;
};

class DoubleAlphaStableExpIntegralTree : public BranchValPtrTree<Dvar<PosReal> >	{

	public:

	DoubleAlphaStableExpIntegralTree(LengthTree* inlengthtree, AlphaStableProcess* inprocess1, AlphaStableProcess* inprocess2)	{
		process1 = inprocess1;
		process2 = inprocess2;
		SetWithRoot(false);
		RecursiveCreate(GetRoot(),GetRoot());
	}

	~DoubleAlphaStableExpIntegralTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return lengthtree->GetTree();}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}

	protected:

	void specialUpdate(Link* from)	{
		if ((! from->isRoot()) || WithRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	void RecursiveCreate(const Link* from, const Link* prev)	{

		if (! from->isRoot())	{
			branchval[from->GetBranch()] = CreateBranchVal(from,prev);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out(),from);
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		cerr << "error in alpha stable exp integral: in CreateBranchVal(const Link* link)\n";
		exit(1);
	}

	Dvar<PosReal>* CreateBranchVal(const Link* link, const Link* prevlink)	{

		return new DoubleAlphaStableExpIntegral(lengthtree->GetBranchVal(link->GetBranch()),process1->GetPath(link->GetBranch()),process1->GetPath(prevlink->GetBranch()),process2->GetPath(link->GetBranch()),process2->GetPath(prevlink->GetBranch()));
	}

	LengthTree* lengthtree;
	AlphaStableProcess* process1;
	AlphaStableProcess* process2;
};

// additional move:
// choose a node at random
// move focal value and then translate all paths around node ?
class AlphaStableNodeMove: public NodeValPtrTree<Mnode>, public MCUpdate {

	public:

	AlphaStableNodeMove(AlphaStableProcess* inprocess, double intuning)	{

		process = inprocess;
		tuning = intuning;
		RecursiveCreate(GetRoot());
	}

	~AlphaStableNodeMove()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	double Move(double tuning_modulator = 1.0)	{
		int n = 0;
		double tot =  RecursiveMove(GetRoot(),tuning * tuning_modulator ,n);
		return tot / n;
	}

	protected:

	double RecursiveMove(const Link* from, double tuning, int& n)	{

		double tot = 0;
		tot += Move(from,tuning);
		n++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveMove(link->Out(),tuning,n);
		}
		return tot;
	}

	double Move(const Link* from, double tuning_mod)	{

		GetNodeVal(from->GetNode())->Corrupt(true);

		double m = tuning * tuning_mod * (Random::Uniform() - 0.5);

		double bk = process->GetPath(from->GetBranch())->GetFinalValue();
		process->GetPath(from->GetBranch())->TranslateFromTip(m);

		if (! from->isLeaf())	{

			double bk1 = process->GetPath(from->Next()->GetBranch())->GetInitValue();
			double bk2 = process->GetPath(from->Next()->Next()->GetBranch())->GetInitValue();

			if ((fabs(bk -bk1) > 1e-6) || (fabs(bk - bk2)>1e-6))	{
				cerr << "error in node move: node values differ\n";
				cerr << bk << '\t' << bk1 << '\t' << bk2 << '\n';
				if (from->isRoot())	{
					cerr << "is root\n";
				}
				exit(1);
			}

			process->GetPath(from->Next()->GetBranch())->TranslateFromBase(m);
			process->GetPath(from->Next()->Next()->GetBranch())->TranslateFromBase(m);
		}

		double logratio = GetNodeVal(from->GetNode())->Update();
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			GetNodeVal(from->GetNode())->Corrupt(false);
			GetNodeVal(from->GetNode())->Restore();
		}
		return (double) accepted;
	}

	Mnode* CreateNodeVal(const Link* link)	{

		Mnode* mnode = new Mnode;
		process->GetBranchVal(link->GetBranch())->Register(mnode);
		if (!link->isLeaf())	{
			process->GetBranchVal(link->Next()->GetBranch())->Register(mnode);
			process->GetBranchVal(link->Next()->Next()->GetBranch())->Register(mnode);
		}
		return mnode;
	}

	AlphaStableProcess* process;
	double tuning;

};

class AlphaStableRescaleBridgeMove : public MCUpdate, public Mnode	{

	public:

	AlphaStableRescaleBridgeMove(AlphaStableProcess* inprocess, Rvar<PosReal>* inbeta, double intuning)	{
		process = inprocess;
		beta = inbeta;
		RecursiveRegister(GetRoot());
		beta->Register(this);
		tuning = intuning;
	}

	double Move(double tuning_mod) {

		Corrupt(true);

		double m = tuning * tuning_mod * (Random::Uniform() - 0.5);
		double e = exp(m);
		int n = process->RescaleBridge(e);
		beta->ScalarMultiplication(1.0/e);
		double logh = m * (n-1);
		double logratio = Update() + logh;

		int accept = (log(Random::Uniform()) < logratio);

		if (! accept)	{
			Corrupt(false);
			Restore();
		}

		return accept;	
	}
		
	protected:

	const Link* GetRoot() {
		return process->GetRoot();
	}

	void RecursiveRegister(const Link* from)	{

		process->GetPath(from->GetBranch())->Register(this);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	AlphaStableProcess* process;
	Rvar<PosReal>* beta;
	double tuning;

};

class AlphaStableRescaleMove : public MCUpdate, public Mnode	{

	public:

	AlphaStableRescaleMove(AlphaStableProcess* inprocess, Rvar<PosReal>* inbeta, double intuning)	{
		process = inprocess;
		beta = inbeta;
		RecursiveRegister(GetRoot());
		beta->Register(this);
		tuning = intuning;
	}

	double Move(double tuning_mod) {

		Corrupt(true);

		double mean = process->GetMean();
		double m = tuning * tuning_mod * (Random::Uniform() - 0.5);
		double e = exp(m);
		int n = process->Rescale(mean,e);
		beta->ScalarMultiplication(1.0/e);
		double logh = m * (n-1);
		double logratio = Update() + logh;

		int accept = (log(Random::Uniform()) < logratio);

		if (! accept)	{
			Corrupt(false);
			Restore();
		}

		return accept;	
	}
		
	protected:

	const Link* GetRoot() {
		return process->GetRoot();
	}

	void RecursiveRegister(const Link* from)	{

		process->GetPath(from->GetBranch())->Register(this);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	AlphaStableProcess* process;
	Rvar<PosReal>* beta;
	double tuning;

};

class AlphaStableAlphaRescaleBridgeMove : public MCUpdate, public Mnode	{

	public:

	AlphaStableAlphaRescaleBridgeMove(AlphaStableProcess* inprocess, Rvar<PosReal>* inalpha, double intuning)	{
		process = inprocess;
		alpha = inalpha;
		RecursiveRegister(GetRoot());
		alpha->Register(this);
		tuning = intuning;
	}

	double Move(double tuning_mod) {

		Corrupt(true);

		double alpha1 = alpha->val();
		double logh = alpha->ProposeMove(tuning * tuning_mod);
		double alpha2 = alpha->val();
		logh += process->AlphaRescaleBridge(alpha1,alpha2);
		double logratio = Update() + logh;

		int accept = (log(Random::Uniform()) < logratio);

		if (! accept)	{
			Corrupt(false);
			Restore();
		}

		return accept;	
	}
		
	protected:

	const Link* GetRoot() {
		return process->GetRoot();
	}

	void RecursiveRegister(const Link* from)	{

		process->GetPath(from->GetBranch())->Register(this);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	AlphaStableProcess* process;
	Rvar<PosReal>* alpha;
	double tuning;

};

class AlphaStableAlphaRescaleMove : public MCUpdate, public Mnode	{

	public:

	AlphaStableAlphaRescaleMove(AlphaStableProcess* inprocess, Rvar<PosReal>* inalpha, double intuning)	{
		process = inprocess;
		alpha = inalpha;
		RecursiveRegister(GetRoot());
		alpha->Register(this);
		tuning = intuning;
	}

	double Move(double tuning_mod) {

		Corrupt(true);

		double mean = process->GetMean();
		double alpha1 = alpha->val();
		double logh = alpha->ProposeMove(tuning * tuning_mod);
		double alpha2 = alpha->val();
		logh += process->AlphaRescale(mean,alpha1,alpha2);
		double logratio = Update() + logh;

		int accept = (log(Random::Uniform()) < logratio);

		if (! accept)	{
			Corrupt(false);
			Restore();
		}

		return accept;	
	}
		
	protected:

	const Link* GetRoot() {
		return process->GetRoot();
	}

	void RecursiveRegister(const Link* from)	{

		process->GetPath(from->GetBranch())->Register(this);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	AlphaStableProcess* process;
	Rvar<PosReal>* alpha;
	double tuning;

};


class ASUgamCompMove : public MCUpdate, public Mnode {

	AlphaStableProcess* tree;
	BranchProcess<PosReal>* ugamtree;
	double tuning;

	public:

	ASUgamCompMove(BranchProcess<PosReal>* inugamtree, AlphaStableProcess* intree, double intuning){
		tree = intree;
		ugamtree = inugamtree;
		tuning = intuning;
		tree->RecursiveRegister(this,tree->GetRoot());
		ugamtree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->Translate(u);
		int n = ugamtree->ScalarMultiplication(m);
		double loghastings = -n * u;
		double logratio = loghastings + Update();
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

#endif
