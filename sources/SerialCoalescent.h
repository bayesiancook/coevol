
#ifndef SERIALCOAL_H
#define SERIALCOAL_H

#include "CalibratedChronogram.h"
#include <algorithm>

class CoalCalDateCopy : public Dvar<PosReal>	{

	public:

	CoalCalDateCopy(CalibratedNodeDate* innodedate, Var<PosReal>* incoalrate, Var<PosReal>* inScale)	{
		nodedate = innodedate;
		coalrate = incoalrate;
		Scale = inScale;

		Register(nodedate);
		Register(coalrate);
		Register(Scale);
	}

	void specialUpdate()	{
		setval(nodedate->val());
	}

	private:

	CalibratedNodeDate* nodedate;
	Var<PosReal>* coalrate;
	Var<PosReal>* Scale;
};

class CoalCalCopyTree : public NodeValPtrTree<CoalCalDateCopy>	{

	public:

	CoalCalCopyTree(CalibratedChronogram* inchrono, Var<PosReal>* incoalrate, Var<PosReal>* inscale)	{
		chrono = inchrono;
		coalrate = incoalrate;
		scale = inscale;
		RecursiveCreate(GetRoot());
	}

	~CoalCalCopyTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return chrono->GetTree();}

	protected:

	CoalCalDateCopy* CreateNodeVal(const Link* link)	{
		return new CoalCalDateCopy(chrono->GetCalibratedNodeDate(link->GetNode()),coalrate,scale);
	}

	private:

	CalibratedChronogram* chrono;
	Var<PosReal>* coalrate;
	Var<PosReal>* scale;

};


class SerialCoalescent: public CalibratedChronogram, public Rnode {

	public:

	SerialCoalescent(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* incoalrate, CalibrationSet* incalibset, double incutoff = 0)	{

		DAGnode::SetName("SerialCoalescent");

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		coalrate = incoalrate;
		calibset = incalibset;
		if (! calibset)	{
			calibset = new CalibrationSet(tree);
		}
		cutoff = incutoff;

		Ntaxa = GetTree()->GetSize(GetRoot());
		// cerr << "Ntaxa : " << Ntaxa << '\n';

		scale = new ChronoScale(rate,1,1,-1,-1,false);

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		copytree = new CoalCalCopyTree(this,coalrate,scale);

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

	~SerialCoalescent()	{
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

	void PushNodeAges(vector<pair<double,const Link*> >& timenodes, const Link* from)	{

		double temp = GetNodeVal(from->GetNode())->val() * scale->val();
		timenodes.push_back(pair<double,const Link*>(temp,from));
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			PushNodeAges(timenodes,link->Out());
		}		
	}

	double logProb()	{

		vector<pair<double,const Link*> > timenodes;
		PushNodeAges(timenodes,GetRoot());
		sort(timenodes.begin(),timenodes.end());
		size_t ntimes = timenodes.size();

		int ntaxa = 1;
		double currenttime = 0;
		size_t k = 1;
		double tot = 0;

		double rate = coalrate->val();

		while ((k < ntimes) && (timenodes[k].first < cutoff))	{
			const Link* link = timenodes[k].second;
			if (! link->isLeaf())	{
				return log(0);
			}
			ntaxa++;
			k++;
		}

		currenttime = cutoff;
		int count = 0;
		while (k < ntimes)	{
			const Link* link = timenodes[k].second;
			double nexttime = timenodes[k].first;
			double deltatime = nexttime - currenttime;
			if (deltatime < 0)	{
				cerr << "in serial coal: negative time : " << deltatime << '\n';
				exit(1);
			}
			if (link->isLeaf())	{
				tot += - 0.5 * ntaxa * (ntaxa-1) * rate * deltatime;
				ntaxa++;
			}
			else	{

				if (ntaxa < 2)	{
					cerr << "error in SerialBDP::logProb: invalid number of lineages: " << ntaxa << "\n";
					exit(1);
				}
				tot += log(rate) - 0.5 * ntaxa * (ntaxa-1) * rate * deltatime;
				count++;
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

		if (std::isinf(tot))	{
			cerr << "error in SerialBDP: log prob is inf\n";
			exit(1);
		}
		if (std::isnan(tot))	{
			cerr << "error in SerialBDP: log prob is nan\n";
			exit(1);
		}

		// Jacobian
		// implicit here: a uniform prior probability of sampling the fossils within the allowed intervals
		// tot += (Ntaxa - 2 + GetNfossil(1.0)) * log(scale->val());

		return tot;
	}

	private:

	unsigned int Ntaxa;
	Var<PosReal>* coalrate;
	CoalCalCopyTree* copytree;
	double cutoff;
};

#endif

