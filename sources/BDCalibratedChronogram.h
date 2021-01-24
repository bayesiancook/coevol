
#ifndef BDCALCHRONO
#define BDCALCHRONO

#include "CalibratedChronogram.h"
#include <algorithm>
#include <utility>
#include <list>

class CalDateCopy : public Dvar<PosReal>	{

	public:

	CalDateCopy(CalibratedNodeDate* innodedate, Var<PosReal>* inChi, Var<PosReal>* inChi2, double inalpha, double inbeta, Var<PosReal>* inScale)	{
		nodedate = innodedate;
		Chi = inChi;
		Chi2 = inChi2;
		alpha = inalpha;
		beta = inbeta;
		Scale = inScale;

		Register(nodedate);
		Register(Chi);
		Register(Chi2);
		Register(Scale);
		specialUpdate();
	}

	double GetLogGG()	{
		ComputeLogGG();
		return loggg;
	}

	void ComputeLogGG()	{
		if (! nodedate->isValueUpdated())	{
			cerr << "error in compute log gg\n";
		}
		// only if node is not calibrated
		loggg = 0;
		if (! nodedate->isCalibrated())	{
			double t0 = Scale->val() * beta / alpha;
			double t = nodedate->val() * t0;
			loggg = BD_logg(Chi->val(),Chi2->val(),t,t0);
			if (std::isinf(loggg))	{
				cerr << "in BD Cal Chrono: inf\n";
				loggg = 1e4;
			}
		}
	}

	void specialUpdate()	{
		// ComputeLogGG();
		setval(nodedate->val());
	}

	private:

	double BD_logg(double p1, double p2, double t, double t0)	{

		double ret = 0;
		if (fabs(p1) > 1e-6)	{
			double e = exp(-p1*t);

			// the P0 below is in fact P0 / rho, where
			// P0 is the probability that a lineage appearing at time 0 is not extinct at time t
			// rho is the sampling fraction
			double P0 = p1 / (p2 + (p1-p2)*e);
			double e0 = exp(-p1*t0);
			double lognu = log(p2) + log(1 - e0) - log(p2*(1-e0) + p1*e0);
			ret = log(p2) - lognu + 2 * log(P0) - p1 * t;

			//
			// I have tried this below, but this seems a bit dangerous, numerically...
			//
			// double nu = p2 * (1-e0) / (p2*(1 - e0) + p1*e0);
			// double expret = p2 / nu * P0 * P0 * exp(-p1 * t);
			// ret = log(expret);
			if (std::isnan(ret))	{
				cerr << "numerical error in BD_logg\n";
				cerr << p1 << '\t' << p2 << '\n';
				cerr << t << '\t' << t0 << '\n';
				cerr << P0 << '\n';
				exit(1);
			}
			if (std::isinf(ret))	{
				cerr << "numerical error in BD_logg\n";
				cerr << p1 << '\t' << p2 << '\n';
				cerr << t << '\t' << t0 << '\n';
				cerr << P0 << '\n';
				exit(1);
			}
		}
		else	{
			ret = log( (1 + p2*t0) / t0 / (1 + p2*t) / (1 + p2*t));
			if (std::isnan(ret))	{
				cerr << "numerical error in BD_logg\n";
				cerr << p1 << '\t' << p2 << '\n';
				cerr << t << '\t' << t0 << '\n';
				exit(1);
			}
			if (std::isinf(ret))	{
				cerr << "numerical error in BD_logg\n";
				cerr << p1 << '\t' << p2 << '\n';
				cerr << t << '\t' << t0 << '\n';
				exit(1);
			}
		}
		return ret;
	}

	double loggg;

	CalibratedNodeDate* nodedate;
	Var<PosReal>* Chi;
	Var<PosReal>* Chi2;
	Var<PosReal>* Scale;
	double alpha;
	double beta;
};

class CalCopyTree : public NodeValPtrTree<CalDateCopy>	{

	public:

	CalCopyTree(CalibratedChronogram* inchrono, Var<PosReal>* inchi, Var<PosReal>* inchi2, double inalpha, double inbeta, Var<PosReal>* inscale)	{
		chrono = inchrono;
		chi = inchi;
		chi2 = inchi2;
		alpha = inalpha;
		beta = inbeta;
		scale = inscale;
		RecursiveCreate(GetRoot());
	}

	~CalCopyTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return chrono->GetTree();}

	double GetTotalLogGG()	{
		return RecursiveGetLogGG(GetRoot());
	}

	protected:

	double RecursiveGetLogGG(const Link* from)	{
		double total = 0;
		if ((! from->isRoot()) && (! from->isLeaf()))	{
			total += GetNodeVal(from->GetNode())->GetLogGG();
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogGG(link->Out());
		}
		return total;
	}

	CalDateCopy* CreateNodeVal(const Link* link)	{
		return new CalDateCopy(chrono->GetCalibratedNodeDate(link->GetNode()),chi,chi2,alpha,beta,scale);
	}

	private:

	CalibratedChronogram* chrono;
	Var<PosReal>* chi;
	Var<PosReal>* chi2;
	double  alpha;
	double beta;
	Var<PosReal>* scale;

};



class BDCalibratedChronogram : public CalibratedChronogram, public Rnode {

	public:

	static int NumErrorCount;

	BDCalibratedChronogram(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* inchi, Var<PosReal>* inchi2, double inalpha, double inbeta, CalibrationSet* incalibset, int inpriortype, double inSofta = 0.025, double inLowerP =0.1, double inLowerC = 1)	{

		DAGnode::SetName("BD Chrono");

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		alpha = inalpha;
		beta = inbeta;
		chi = inchi;
		chi2 = inchi2;
		calibset = incalibset;
		priortype = inpriortype;
		if (priortype > 3)	{
			cerr << "in bd cal chrono: check soft bounds\n";
			exit(1);
		}
		ImproperLowerBound = false;
		if (priortype == 5)	{
			priortype = 4;
			ImproperLowerBound = true;
		}
		Softa = inSofta;
		LowerP = inLowerP;
		LowerC = inLowerC;
		// 1 : improper lower bounds
		// 2 : cauchy lower bounds
		if (priortype == 4)	{
			cerr << "soft bounds\n";
		}

		if ((priortype == 2) || (priortype == 4))	{
			p = 0.1;
			c = 1;
			A = 0.5 + atan(p/c) / Pi;
		}

		Ntaxa = GetTree()->GetSize(GetRoot());
		// cerr << "Ntaxa : " << Ntaxa << '\n';

		double rootmin = -1;
		double rootmax = -1;
		if (isCalibrated(GetRoot()->GetNode()))	{
			rootmin = GetLowerCalibration(GetRoot()->GetNode());
			rootmax = GetUpperCalibration(GetRoot()->GetNode());
		}

		if (isCalibrated(GetRoot()->GetNode()) && (priortype > 1))	{
			scale = new ChronoScale(rate,alpha,beta,rootmin,rootmax,false);
		}
		else	{
			scale = new ChronoScale(rate,alpha,beta,rootmin,rootmax,true);
		}
		// scale = new Gamma(alpha,beta);

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());
		// logggtree = new LogGGTree(this,chi,chi2,alpha,beta,scale);

		copytree = new CalCopyTree(this,chi,chi2,alpha,beta,scale);

		RecursiveRegisterProb(GetRoot());

		/*
		Register(chi);
		Register(chi2);
		Register(scale);
		Register(alpha);
		Register(beta);
		Register(rate);
		*/


		RecursiveSetCalibrations(GetRoot());

		RecursiveYoungerLimit(GetRoot());
		RecursiveOlderLimit(GetRoot(),-1);
		double d = RecursiveDrawAges(GetRoot());

		// cerr << "age  " << d << '\n';

		RecursiveNormalizeTree(GetRoot(),d,false);
		RecursiveUpdateBranches(GetRoot());

		scale->setval(d);
		map<DAGnode*,int> tmpmap;
		scale->FullCorrupt(tmpmap);
		scale->FullUpdate();

		double tmp = RecursiveGetLogProb(GetRoot());
		if (fabs(tmp) > 1e-6)	{
			cerr << "error : nodes out of calib\n";
			exit(1);
		}
	}

	~BDCalibratedChronogram()	{
        /*
		RecursiveDeleteBranch(GetRoot());
		RecursiveDeleteNode(GetRoot());
        */
        delete copytree;
		delete scale;
	}

	void Sample()	{
		CalibratedChronogram::Sample();
	}

	void drawSample()	{
		// CalibratedChronogram::Sample();
	}

	double Move(double tuning){
		return CalibratedChronogram::Move(tuning);
	}

	double GetLogProb()	{
		return Rnode::GetLogProb();
	}

	double ProposeMove(double)	{
		cerr << "error in BDCalibratedChronogram: in propose move\n";
		exit(1);
		return 0;
	}

	double GetTotalLogGG()	{
		return copytree->GetTotalLogGG();
	}

	double logProb()	{
		double total = 0;
		if (priortype == 2)	{
			total += logCauchyCalibPrior();
		}
		else if (priortype == 3)	{
			total += logRootBoundedCalibPrior();
		}
		else if (priortype == 4)	{
			total += logSoftBoundedCalibPrior();
		}

		total += logJacobian();
		// total += RecursiveGetLogProb(GetRoot());
		total += GetTotalLogGG();

		total += logAbsBDUncalibratedPrior();
		return total;
	}

	void RecursiveRegisterProb(const Link* from)	{
		// Register(GetNodeVal(from->GetNode()));
		// Register(logggtree->GetNodeVal(from->GetNode()));
		Register(copytree->GetNodeVal(from->GetNode()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegisterProb(link->Out());
		}
	}

	double logSoftBoundedCalibPrior()	{
		// only for lower bounds:  Cauchy
		return RecursiveLogSoftBoundedCalibPrior(GetRoot());
	}

	double RecursiveLogSoftBoundedCalibPrior(const Link* from)	{

		double total = 0;
		double a = Softa;
		double t = GetNodeVal(from->GetNode())->val() * scale->val();

		if ((isLowerCalibrated(from->GetNode())) && (isUpperCalibrated(from->GetNode())))	{
			double t_l = GetLowerCalibration(from->GetNode());
			double t_u = GetUpperCalibration(from->GetNode());

			double theta_l = (1 - 2*a) * t_l / a / (t_u - t_l);
			double theta_u = theta_l / t_l;
			if (t > t_u)	{
				total -= log(a) + log(theta_u) - theta_u * (t - t_u);
			}
			else if (t < t_l)	{
				total -= log(a) + log(theta_u) + (theta_l - 1) * log(t / t_l);
			}
			else	{
				total -= log(1-2*a) - log(t_u - t_l);
			}
			if (std::isnan(total))	{
				cerr << "error in log calib: log prob nan\n";
				cerr << "lower and upper bound\n";
				cerr << (t<t_l) << '\t' << (t>t_u) << '\n';
				cerr << theta_u << '\t' << theta_l << '\t' << t_u << '\t' << t_l << '\n';
				exit(1);
			}
		}
		if ((isLowerCalibrated(from->GetNode())) && (!isUpperCalibrated(from->GetNode())))	{

			double t_l = GetLowerCalibration(from->GetNode());

			if (ImproperLowerBound)	{ // improper
				double b = 0.1;
				double theta_l = log(b) / log(1-b);
				if (t < t_l)	{
					total -= log(a) + log(theta_l / t_l) + (theta_l - 1) *log(t / t_l);
				}
				else	{
					total -= log(1 - a) + log(theta_l / t_l);
				}
			}
			else	{
				double p = LowerP;
				double c = LowerC;
				double A = 0.5 + atan(p/c) / Pi;
				// double A = 0.5 + cos(p/c) / sin(p/c) / Pi;
				double theta = (1-a) / a / Pi / A / c / (1 + (p*p/c/c));
				double f = 0;
				double x = 0;
				if (t < t_l)	{
					total -= log(a) - theta * log(t_l) + log(theta) + (theta-1) * log(t);
				}
				else	{
					f = (t - t_l * (1 + p)) / c / t_l;
					x = log(1 - a) - log(A) - log(Pi) - log(c) - log(t_l) - log(1 + f*f);
					total -= x;
				}
				if (std::isnan(total))	{
					cerr << a << '\t' << p << '\t' << c << '\t' << A << '\t' << Pi << '\t' << theta << '\t' << t << '\t' << t_l  << '\t' << f << '\n';
					cerr << log(theta) << '\n';
					cerr << log(1 + f*f) << '\n';
					cerr << log(1 - a) << '\n';
					cerr << x << '\n';

					exit(1);
				}
			}
		}
		else if ((!isLowerCalibrated(from->GetNode())) && (isUpperCalibrated(from->GetNode())))	{

			double t_u = GetUpperCalibration(from->GetNode());

			double theta_u = (1-a)/a/t_u;
			if (t > t_u)	{
				total -= log(a) + log(theta_u) - theta_u * (t - t_u);
			}
			else	{
				total -= log(1 - a) + log(theta_u);
			}
		}
		if (std::isinf(total))	{
			cerr << "error in log calib: log prob inf\n";
			exit(1);
		}
		if (std::isnan(total))	{
			cerr << "error in log calib: log prob nan\n";
			exit(1);
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveLogCauchyCalibPrior(link->Out());
		}
		return total;
	}

	double logCauchyCalibPrior()	{
		// only for lower bounds:  Cauchy
		return RecursiveLogCauchyCalibPrior(GetRoot());
	}

	double RecursiveLogCauchyCalibPrior(const Link* from)	{
		double total = 0;
		if ((isLowerCalibrated(from->GetNode())) && (! isUpperCalibrated(from->GetNode())))	{
			double t = GetNodeVal(from->GetNode())->val() * scale->val();
			double t_l = GetLowerCalibration(from->GetNode());
			double f = (t - t_l * (1 + p)) / c / t_l;
			double x = - log(A) - log(Pi) - log(c) - log(t_l) - log(1 + f*f);
			total += x;
			if (std::isnan(total))	{
				cerr << "error in BDCalibChrono: logCalibPrior\n";
				exit(1);
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveLogCauchyCalibPrior(link->Out());
		}
		return total;
	}

	double logRootBoundedCalibPrior()	{
		// only for lower bounds:  Cauchy
		return RecursiveLogRootBoundedCalibPrior(GetRoot());
	}

	double RecursiveLogRootBoundedCalibPrior(const Link* from)	{
		double total = 0;
		if (isCalibrated(from->GetNode()))	{
			double upper = GetUpperCalibration(from->GetNode());
			if (upper == -1)	{
				upper = scale->val();
			}
			double lower = GetLowerCalibration(from->GetNode());
			if (lower == -1)	{
				lower = 0;
			}
			double t = GetNodeVal(from->GetNode())->val() * scale->val();
			if ((t > upper) || (t < lower))	{
				// cerr << "error in log Root bounded calib prior: out of bound\n";
				// cerr << t << '\t' << upper << '\t' << lower << '\t' << scale->val() << '\t' << GetUpperCalibration(from->GetNode()) << '\t' << GetLowerCalibration(from->GetNode()) << '\n';
				// exit(1);
				// total -= 1e6;
			}
			else	{
				total -= log(upper-lower);
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveLogRootBoundedCalibPrior(link->Out());
		}
		return total;
	}

	double logJacobian()	{
		return (Ntaxa-2)* log(scale->val());
	}

	typedef pair<double,bool>  doublet;
	bool compare_doublet(doublet d1, doublet d2) {return d1.first < d2.first;}

	void RecursiveGetAgeList(const Link* from, list<doublet>& agelist)	{
		if (! from->isLeaf())	{
			bool tmp = from->isRoot() || isCalibrated(from->GetNode());
			double temp = GetNodeVal(from->GetNode())->val();
			agelist.push_back(doublet(temp, tmp));
		}
		int degree = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetAgeList(link->Out(),agelist);
			degree++;
		}
		if (degree > 2)	{
			cerr << "error : BD accepts only bifurcating trees\n";
			exit(1);
		}
	}


	double logAbsBDUncalibratedPrior()	{


		// sort calibrated nodes

		double Scale0 = alpha / beta;

		list<doublet> agelist;
		RecursiveGetAgeList(GetRoot(),agelist);
		/*
		if (agelist.size() != Ntaxa-1)	{
			cerr << "error in BD log prob: agelist size : " << agelist.size() << '\n';
			exit(1);
		}
		*/

		agelist.sort();
		// agelist.sort(compare_doublet);
		// sort(agelist.begin(),agelist.end(),compare_doublet);

		double t0 = scale->val() / Scale0;
		double total = 0;

		// total += logggtree->GetTotalLogGG();
		// total += copytree->GetTotalLogGG();

		int N = calibset->GetNCalib();
		if (isCalibrated(GetRoot()->GetNode()))	{
			N--;
		}
		/*
		if (bk)	{
			N--;
		}
		*/

		double T[N+2];
		int I[N+2];
		I[0] = 0;
		T[0] = 0;
		int k = 1;
		int j = 1;
		for (list<doublet>::iterator i=agelist.begin(); i!= agelist.end(); i++)	{
			if (i->second) {
				if (k == N+2)	{
					cerr << "error in BD log prob\n";
					exit(1);
				}
				I[k] = j;
				T[k] = i->first * t0;
				k++;
			}
			j++;
		}

		/*
		if (k != N+2)	{
			cerr << "error : non matching number: " << k << '\t' << N+2 << '\n';
			exit(1);
		}
		if (I[N+1] != Ntaxa -1)	{
			cerr << "error : " << I[N+1]  << '\t' << Ntaxa -1 << '\n';
			exit(1);
		}
		if (T[N+1] != t0)	{
			cerr << "error for T : " << T[N+1] << '\t' << t0 << '\n';
			exit(1);
		}
		*/

		for (int i=0; i<N+1; i++)	{
			int d = I[i+1] - I[i] - 1;
			for (int k=2; k<=d; k++)	{
				total += log(k);
			}
			double tmp = BD_logDeltaG(chi->val(),chi2->val(),T[i],T[i+1],t0);
			if (std::isinf(tmp) || std::isnan(tmp))	{
				cerr << "chi2 : " << *chi2 << '\n';
				exit(1);
			}
			total -= d * tmp;
		}

		return total;
	}

	double BD_logDeltaG(double p1, double p2, double t1, double t2, double t0)	{

		double ret = 0;
		if (fabs(t2-t1) < 1e-10)	{
			return 0;
		}
		if (fabs(p1) > 1e-6)	{
			double e0 = exp(-p1*t0);
			double lognu = log(p2) + log(1 - e0) - log(p2*(1-e0) + p1*e0);

			ret = (log(p1) + log(p2)) - lognu;
			ret -= p1*t1;
			ret += log(1 - exp(-p1*(t2-t1)));
			ret -= log(p2*(1 - exp(-p1*t1)) + p1*exp(-p1*t1));
			ret -= log(p2*(1 - exp(-p1*t2)) + p1*exp(-p1*t2));
			if (std::isinf(ret) || std::isnan(ret))	{
				cerr << "numerical error in BN_logDeltaG\n";
				cerr << p1 << '\t' << p2 << '\n';
				cerr << log(p2) << '\n';
				cerr << t1 << '\t' << t2 << '\t' << t0 << '\n';
				cerr << p1*t1 << '\n';
				cerr << log(1 - exp(-p1*(t2-t1))) << '\n';
				cerr << log(p2*(1 - exp(-p1*t1)) + p1*exp(-p1*t1)) << '\n';
				cerr << log(p2 + (p1-p2)*exp(-p1*t2)) << '\n';
				exit(1);
				ret = 0;
				NumErrorCount++;
			}
		}
		else	{
			ret = (1 + p2*t0) * t2 / t0 / (1 + p2*t2) ;
			ret -= (1 + p2*t0) * t1 / t0 / (1 + p2*t1) ;
			if (std::isinf(ret) || std::isnan(ret))	{
				NumErrorCount++;
				ret = 1;
				cerr << "numerical error in BN_logDeltaG\n";
				exit(1);
			}
			ret = log(ret);
		}
		return ret;
	}

	private:

	unsigned int Ntaxa;
	Var<PosReal>* chi;
	Var<PosReal>* chi2;
	CalCopyTree* copytree;

	// Cauchy lower bound
	int priortype;
	double p,c,A;

	double Softa;
	double LowerP;
	double LowerC;
	bool ImproperLowerBound;

};

int BDCalibratedChronogram::NumErrorCount = 0;


#endif

