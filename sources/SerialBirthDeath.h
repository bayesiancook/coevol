
#ifndef SERIALBD_H
#define SERIALBD_H

#include "CalibratedChronogram.h"
#include <algorithm>

class SerialBDCalDateCopy : public Dvar<PosReal>	{

	public:

	SerialBDCalDateCopy(CalibratedNodeDate* innodedate, Var<PosReal>* indivrate, Var<PosReal>* inmu, Var<PosReal>* inpsi, Var<UnitReal>* inrho, Var<PosReal>* inScale)	{
		nodedate = innodedate;
		divrate = indivrate;
		mu = inmu;
		psi = inpsi;
		rho = inrho;
		Scale = inScale;

		Register(nodedate);
		Register(divrate);
		Register(mu);
		Register(psi);
		Register(rho);
		Register(Scale);
	}

	void specialUpdate()	{
		setval(nodedate->val());
	}

	private:

	CalibratedNodeDate* nodedate;
	Var<PosReal>* divrate;
	Var<PosReal>* mu;
	Var<PosReal>* psi;
	Var<UnitReal>* rho;
	Var<PosReal>* Scale;
};

class SerialBDCalCopyTree : public NodeValPtrTree<SerialBDCalDateCopy>	{

	public:

	SerialBDCalCopyTree(CalibratedChronogram* inchrono, Var<PosReal>* indivrate, Var<PosReal>* inmu, Var<PosReal>* inpsi, Var<UnitReal>* inrho, Var<PosReal>* inscale)	{
		chrono = inchrono;
		divrate = indivrate;
		mu = inmu;
		psi = inpsi;
		rho = inrho;
		scale = inscale;
		RecursiveCreate(GetRoot());
	}

	~SerialBDCalCopyTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return chrono->GetTree();}

	protected:

	SerialBDCalDateCopy* CreateNodeVal(const Link* link)	{
		return new SerialBDCalDateCopy(chrono->GetCalibratedNodeDate(link->GetNode()),divrate,mu,psi,rho,scale);
	}

	private:

	CalibratedChronogram* chrono;
	Var<PosReal>* divrate;
	Var<PosReal>* mu;
	Var<PosReal>* psi;
	Var<UnitReal>* rho;
	Var<PosReal>* scale;

};


class SerialBirthDeath : public CalibratedChronogram, public Rnode {

	public:

	SerialBirthDeath(Tree* intree, Var<PosReal>* inrate, double inalpha, double inbeta, Var<PosReal>* indivrate, Var<PosReal>* inmu, Var<PosReal>* inpsi, Var<UnitReal>* inrho, CalibrationSet* incalibset, double incutoff = 0, int inNextant = 0, bool origin = false)	{

		DAGnode::SetName("SerialBirthDeath");

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		divrate = indivrate;
		mu = inmu;
		psi = inpsi;
		rho = inrho;
		calibset = incalibset;
		cutoff = incutoff;
		Nextant = inNextant;

		Ntaxa = GetTree()->GetSize(GetRoot());
		// cerr << "Ntaxa : " << Ntaxa << '\n';

		scale = new ChronoScale(rate,inalpha,inbeta,-1,-1,true);

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		copytree = new SerialBDCalCopyTree(this,divrate,mu,psi,rho,scale);

		RecursiveRegisterProb(GetRoot());
		RecursiveSetCalibrations(GetRoot());

		double maxage = RecursiveSetNodesValues(GetRoot());
		// RecursiveEqualizeLeafNodes(GetRoot(),maxage);
		RecursiveNormalizeTree(GetRoot(),maxage,true);

		RecursiveUpdateBranches(GetRoot());

		scale->setval(maxage);
		map<DAGnode*,int> tmpmap;
		scale->FullCorrupt(tmpmap);
		scale->FullUpdate();

		double tmp = RecursiveGetLogProb(GetRoot());
		if (fabs(tmp) > 1e-6)	{
			cerr << "error : nodes out of calib\n";
			cerr << tmp << '\n';
			cerr << scale->GetLogProb() << '\n';
			exit(1);
		}
	}

	~SerialBirthDeath()	{
		RecursiveDeleteBranch(GetRoot());
		RecursiveDeleteNode(GetRoot());
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
		cerr << "error in SerialCoalescent: in propose move\n";
		exit(1);
		return 0;
	}

	void RecursiveRegisterProb(const Link* from)	{
		Register(copytree->GetNodeVal(from->GetNode()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegisterProb(link->Out());
		}
	}

	double pSurvival(double end)	{

		double t = end;
		// get the parameters
		double birth = divrate->val() + mu->val();
		double death = mu->val();
		double p     = psi->val();
		double r     = rho->val();

		double a = (birth-death-p);
		double c1 = sqrt( a*a + 4*birth*p );
		double c2 = - (a-2.0*birth*r)/c1;

		double oneMinusC2 = 1.0-c2;
		double onePlusC2  = 1.0+c2;

		double e = exp(-c1*t);

		double p0 = ( birth + death + p + c1 * ( e * oneMinusC2 - onePlusC2 ) / ( e * oneMinusC2 + onePlusC2 ) ) / (2.0 * birth);

		if (p0 < 0)	{
			cerr << "error: negative extinction prob\n";
			cerr << p0 << '\n';
			exit(1);
		}
		if (std::isinf(p0))	{
			cerr << "error in serial bd: p0 is inf\n";
			exit(1);
		}

		if (std::isnan(p0))	{
			cerr << "error in serial bd: p0 is nan\n";
			exit(1);
		}

		return 1.0 - p0;
	}

	double q( double t )	{

		// get the parameters
		double birth = divrate->val() + mu->val();
		double death = mu->val();
		double p     = psi->val();
		double r     = rho->val();

		double a = (birth-death-p);
		double c1 = sqrt( a*a + 4*birth*p );
		double c2 = - (a-2.0*birth*r)/c1;

		double oneMinusC2 = 1.0-c2;
		double onePlusC2  = 1.0+c2;

		double ct = -c1*t;
		double b = onePlusC2+oneMinusC2*exp(ct);
		double tmp = exp( ct ) / (b*b);

		if (std::isinf(tmp))	{
			cerr << "error in serial bd q function: inf\n";
			exit(1);
		}
		if (std::isnan(tmp))	{
			cerr << "error in serial bd q function: nan\n";
			exit(1);
		}

		return  tmp;    
	}

	// conditioning on mrca

	double RecursiveGetBDLogProb(const Link* from)	{

		double ret = 0;

		double age = GetNodeVal(from->GetNode())->val() * scale->val();

		if (from->isLeaf())	{

			double precision = 1.0;
			// fossil
			if (age > precision)	{
				if (age < cutoff)	{
					return log(0);
					/*
					cerr << "error: fossils sampled within cutoff time\n";
					cerr << age << '\t' << cutoff << '\n';
					cerr << from->GetNode()->GetName() << '\n';
					exit(1);
					*/
				}
				ret =  log(psi->val()) - log(q(age));
			}

			// extant taxon
			else	{
				if (cutoff && (! Nextant))	{
					ret = log(pSurvival(cutoff)) - log(q(cutoff));
				}
				else	{
					ret = log(4.0 * rho->val());
				}
			}
		}
		else	{
			if (age < cutoff)	{
				return log(0);
			}
			else	{
				ret = log(divrate->val() + mu->val()) + log(q(age));
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ret += RecursiveGetBDLogProb(link->Out());
		}
		return ret;
	}

	double logProb()	{

		double tot = RecursiveGetBDLogProb(GetRoot());

		double birth = divrate->val() + mu->val();
		double death = mu->val();
		double div = divrate->val();

		// conditioning on mrca
		tot -= log(birth);
		tot += log(q(scale->val()));
		tot -= 2 * log(pSurvival(scale->val()));

		// diversified sampling, lambert and stadler's style
		if (Nextant)	{
			double expo = exp(- div * cutoff);
			double tmp = log((birth * (1 - expo)) / (birth - death * expo));
			if (std::isnan(tmp))	{
				return log(0);
			}
			tot += (Nextant - Ntaxa) * tmp;
		}

		// Jacobian
		// tot += (Ntaxa-2 + GetNfossil(1.0)) * log(scale->val());

		/*
		if (isinf(tot))	{
			cerr << "error in SerialBDP: log prob is inf\n";
			exit(1);
		}
		*/
		if (std::isnan(tot))	{
			cerr << "error in SerialBDP: log prob is nan\n";
			cerr << birth << '\t' << death << '\t' << psi->val()  << '\n';
			exit(1);
		}

		return tot;
	}

	private:

	unsigned int Ntaxa;
	Var<PosReal>* divrate;
	Var<PosReal>* mu;
	Var<PosReal>* psi;
	Var<UnitReal>* rho;
	SerialBDCalCopyTree* copytree;
	double cutoff;
	int Nextant;
};

#endif

