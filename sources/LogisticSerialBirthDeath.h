
#ifndef LOGISTICSERIALBD_H
#define LOGISTICSERIALBD_H

#include "CalibratedChronogram.h"
#include <algorithm>

class LogBDCalDateCopy : public Dvar<PosReal>	{

	public:

	LogBDCalDateCopy(CalibratedNodeDate* innodedate, Var<PosReal>* indivrate, Var<PosReal>* inextrate, Var<PosReal>* inK0, Var<UnitReal>* inmassext, Var<PosReal>* inK1, Var<PosReal>* inT0, Var<PosReal>* inT1, Var<PosReal>* inpsi, Var<PosReal>* inScale)	{
		nodedate = innodedate;
		divrate = indivrate;
		extrate = inextrate;
		K0 = inK0;
		massext = inmassext;
		K1 = inK1;
		T0 = inT0;
		T1 = inT1;
		psi = inpsi;
		Scale = inScale;

		Register(nodedate);
		Register(divrate);
		Register(extrate);
		Register(K0);
		Register(massext);
		Register(K1);
		Register(T0);
		Register(T1);
		Register(psi);
		Register(Scale);
	}

	void specialUpdate()	{
		setval(nodedate->val());
	}

	private:

	CalibratedNodeDate* nodedate;
	Var<PosReal>* divrate;
	Var<PosReal>* extrate;
	Var<PosReal>* K0;
	Var<UnitReal>* massext;
	Var<PosReal>* K1;
	Var<PosReal>* T0;
	Var<PosReal>* T1;
	Var<PosReal>* psi;
	Var<PosReal>* Scale;
};

class LogBDCalCopyTree : public NodeValPtrTree<LogBDCalDateCopy>	{

	public:

	LogBDCalCopyTree(CalibratedChronogram* inchrono, Var<PosReal>* indivrate, Var<PosReal>* inextrate, Var<PosReal>* inK0, Var<UnitReal>* inmassext, Var<PosReal>* inK1, Var<PosReal>* inT0, Var<PosReal>* inT1, Var<PosReal>* inpsi, Var<PosReal>* inScale)	{
		chrono = inchrono;
		divrate = indivrate;
		extrate = inextrate;
		K0 = inK0;
		massext = inmassext;
		K1 = inK1;
		T0 = inT0;
		T1 = inT1;
		psi = inpsi;
		Scale = inScale;
		RecursiveCreate(GetRoot());
	}

	~LogBDCalCopyTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return chrono->GetTree();}

	protected:

	LogBDCalDateCopy* CreateNodeVal(const Link* link)	{
		return new LogBDCalDateCopy(chrono->GetCalibratedNodeDate(link->GetNode()),divrate,extrate,K0,massext,K1,T0,T1,psi,Scale);
	}

	private:

	CalibratedChronogram* chrono;
	Var<PosReal>* divrate;
	Var<PosReal>* extrate;
	Var<PosReal>* K0;
	Var<UnitReal>* massext;
	Var<PosReal>* K1;
	Var<PosReal>* T0;
	Var<PosReal>* T1;
	Var<PosReal>* psi;
	Var<PosReal>* Scale;

};


class LogisticSerialBirthDeath: public CalibratedChronogram, public Rnode {

	public:

	LogisticSerialBirthDeath(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* indivrate, Var<PosReal>* inextrate, Var<PosReal>* inK0, Var<UnitReal>* inmassext, Var<PosReal>* inK1, Var<PosReal>* inT0, Var<PosReal>* inT1, Var<PosReal>* inpsi, CalibrationSet* incalibset, double incutoff = 0, int inNextant = 0)	{

		DAGnode::SetName("SerialBirthDeath");

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		divrate = indivrate;
		extrate = inextrate;
		K0 = inK0;
		massext = inmassext;
		K1 = inK1;
		T0 = inT0;
		T1 = inT1;
		psi = inpsi;
		calibset = incalibset;
		if (! calibset)	{
			calibset = new CalibrationSet(tree);
		}
		cutoff = incutoff;
		Nextant = inNextant;

		Ntaxa = GetTree()->GetSize(GetRoot());
		// cerr << "Ntaxa : " << Ntaxa << '\n';

		scale = new ChronoScale(rate,1,1,-1,-1,false);

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		copytree = new LogBDCalCopyTree(this,divrate,extrate,K0,massext,K1,T0,T1,psi,scale);

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
			exit(1);
		}
	}

	~LogisticSerialBirthDeath()	{
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
		cerr << "error in SerialBirthDeath: in propose move\n";
		exit(1);
		return 0;
	}

	void RecursiveRegisterProb(const Link* from)	{
		// Register(GetNodeVal(from->GetNode()));
		// Register(logggtree->GetNodeVal(from->GetNode()));
		Register(copytree->GetNodeVal(from->GetNode()));
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegisterProb(link->Out());
		}
	}

	void PushNodeAges(vector<pair<double,const Link*> >& timenodes, const Link* from)	{

		double temp = GetNodeVal(from->GetNode())->val() * scale->val();
		timenodes.push_back(pair<double,const Link*>(temp,from));
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			PushNodeAges(timenodes,link->Out());
		}		
	}

	double GetT0()	{
		return T0->val() + scale->val();
	}

	double GetDiversityRightBeforeKT()	{
		return K0->val() / (1 + (K0->val() - 1) * exp(- divrate->val() * (GetT0() - T1->val())));
	}

	double GetDiversityRightAfterKT()	{
		return massext->val() * GetDiversityRightBeforeKT();
	}

	double GetMRightAfterKT()	{
		double N = GetDiversityRightAfterKT();
		return (K1->val() - N) / N;
	}

	double GetDiversity(double time)	{

		if (time < T1->val())	{
			double m1 = GetMRightAfterKT();
			return K1->val() / (1 + m1 * exp(-divrate->val() * (T1->val() - time)));
		}
		return K0->val() / (1 + (K0->val() - 1) * exp(-divrate->val() * (GetT0() - time)));
	}

	double GetCarryingCapacity(double time)	{

		if (time < T1->val())	{
			return K1->val();
		}
		return K0->val();
	}

	double GetRate(double time)	{
		return (divrate->val() + extrate->val()) / GetDiversity(time) - divrate->val() / GetCarryingCapacity(time);
	}

	double GetCumulRate(double t1, double t2)	{

		if ((t1 < T1->val()) && (t2 > T1->val()))	{
			cerr << "error: cumulative rate across KT boundary\n";
			exit(1);
		}

		if (t2 < t1)	{
			cerr << "error: cumul rate negative time interval\n";
			exit(1);
		}

		if (t2 > GetT0())	{
			cerr << "error: cumul rate t2 > t0\n";
			exit(1);
		}

		double tot = 0;
		if (t2 < T1->val())	{

			tot += extrate->val() * (t2 - t1);
			tot += (divrate->val() + extrate->val()) / divrate->val() * GetMRightAfterKT() * exp(-rate->val() * (T1->val() - t2)) * (1 - exp(-rate->val() * (t2 - t1)));

			tot /= K1->val();
		}

		else	{

			tot += extrate->val() * (t2 - t1);
			tot += (divrate->val() + extrate->val()) / divrate->val() * (K0->val() - 1) * exp(-rate->val() * (GetT0() - T1->val() - t2)) * (1 - exp(-rate->val() * (t2 - t1)));

			tot /= K0->val();
		}
			
		return tot;
	}

	virtual int CheckBounds()	{
		CalibratedChronogram::CheckBounds();
		
		vector<pair<double,const Link*> > timenodes;
		PushNodeAges(timenodes,GetRoot());
		timenodes.push_back(pair<double,const Link*>(T1->val(),0));
		// timenodes.push_back(pair<double,const Link*>(T0->val(),0));
		sort(timenodes.begin(),timenodes.end());
		size_t ntimes = timenodes.size();

		int ntaxa = 1;
		double currenttime = 0;
		size_t k = 1;
		double tot = 0;

		while ((k < ntimes) && (timenodes[k].first < cutoff))	{
			const Link* link = timenodes[k].second;
			if (! link)	{
				cerr << "error: cutoff is larger than KT\n";
				exit(1);
			}
			if (! link->isLeaf())	{
				cerr << "in serial coal: internal node within diversified sampling zone\n";
			}
			ntaxa++;
			k++;
		}

		return 0;
	}

	double logProb()	{

		vector<pair<double,const Link*> > timenodes;
		PushNodeAges(timenodes,GetRoot());
		timenodes.push_back(pair<double,const Link*>(T1->val(),0));
		// timenodes.push_back(pair<double,const Link*>(T0->val(),0));
		sort(timenodes.begin(),timenodes.end());
		size_t ntimes = timenodes.size();

		/*
		if (T0->val() < timenodes[ntimes-1].first)	{
			cerr << "T0 too small\n";
			return -100000;
		}
		*/

		int ntaxa = 1;
		double currenttime = 0;
		size_t k = 1;
		double tot = 0;

		while ((k < ntimes) && (timenodes[k].first < cutoff))	{
			const Link* link = timenodes[k].second;
			if (! link)	{
				cerr << "error: cutoff is larger than KT\n";
				exit(1);
			}
			if (! link->isLeaf())	{
				// cerr << "in serial coal: internal node within diversified sampling zone\n";
				return log(0);
				// -1000000;
			}
			ntaxa++;
			k++;
		}

		currenttime = cutoff;
		while (k < ntimes)	{
			const Link* link = timenodes[k].second;
			double nexttime = timenodes[k].first;
			if (! link)	{
				tot += - 0.5 * ntaxa * (ntaxa-1) * GetCumulRate(currenttime,nexttime);
			}
			else if (link->isLeaf())	{
				tot += - 0.5 * ntaxa * (ntaxa-1) * GetCumulRate(currenttime,nexttime);
				ntaxa++;
			}
			else	{

				if (ntaxa < 2)	{
					cerr << "error in SerialBDP::logProb: invalid number of lineages: " << ntaxa << "\n";
					exit(1);
				}
				tot += log(GetRate(nexttime)) - 0.5 * ntaxa * (ntaxa-1) * GetCumulRate(currenttime,nexttime);
				// tot += log(0.5 * ntaxa * (ntaxa-1) * rate) - 0.5 * ntaxa * (ntaxa-1) * rate * deltatime;
				// choose 2 among ntaxa
				// tot -= log(0.5 * ntaxa * (ntaxa-1));

				ntaxa--;
			}

			currenttime = nexttime;
			k++;
		}

		if (ntaxa != 1)	{
			cerr << "error in SerialBDP: backward processing does not end with 1 lineage : " << ntaxa << '\n';
			exit(1);
		}

		if (isinf(tot))	{
			cerr << "error in SerialBDP: log prob is inf\n";
			exit(1);
		}
		if (isnan(tot))	{
			cerr << "error in SerialBDP: log prob is nan\n";
			exit(1);
		}

		// to account for Hastings ratio
		tot += (Ntaxa-2) * log(scale->val());

		return tot;
	}

	private:

	unsigned int Ntaxa;
	Var<PosReal>* divrate;
	Var<PosReal>* extrate;
	Var<PosReal>* K0;
	Var<UnitReal>* massext;
	Var<PosReal>* K1;
	Var<PosReal>* T0;
	Var<PosReal>* T1;
	Var<PosReal>* psi;
	LogBDCalCopyTree* copytree;
	double cutoff;
	int Nextant;
};

#endif
