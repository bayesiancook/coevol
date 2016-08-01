
#ifndef COALCHRONO
#define COALCHRONO

#include "CalibratedChronogram.h"
#include <algorithm>
#include <utility>
#include <list>

class CoalDateCopy : public Dvar<PosReal>	{

	public:

	CoalDateCopy(CalibratedNodeDate* innodedate, Var<PosReal>* ind, Var<PosReal>* inr, Var<PosReal>* inN0, Var<PosReal>* inN1, Var<PosReal>* inN2, double inT0)	{
		nodedate = innodedate;
		d = ind;
		r = inr;
		N0 = inN0;
		N1 = inN1;
		N2 = inN2;
		T0 = inT0;

		Register(nodedate);
		Register(r);
		Register(d);
		Register(N0);
		Register(N1);
		Register(N2);
	}

	void specialUpdate()	{
		setval(nodedate->val());
	}

	private:

	CalibratedNodeDate* nodedate;
	Var<PosReal>* d;
	Var<PosReal>* r;
	Var<PosReal>* N0;
	Var<PosReal>* N1;
	Var<PosReal>* N2;
	double T0;
};

class CoalCopyTree : public NodeValPtrTree<CoalDateCopy>	{

	public:

	CoalCopyTree(CalibratedChronogram* inchrono, Var<PosReal>* ind, Var<PosReal>* inr, Var<PosReal>* inN0, Var<PosReal>* inN1, Var<PosReal>* inN2, double inT0)	{
		chrono = inchrono;
		d = ind;
		r = inr;
		N0 = inN0;
		N1 = inN1;
		N2 = inN2;
		T0 = inT0;
		RecursiveCreate(GetRoot());
	}

	~CoalCopyTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return chrono->GetTree();}

	protected:

	CoalDateCopy* CreateNodeVal(const Link* link)	{
		return new CoalDateCopy(chrono->GetCalibratedNodeDate(link->GetNode()),d,r,N0,N1,N2,T0);
	}

	private:

	CalibratedChronogram* chrono;
	Var<PosReal>* d;
	Var<PosReal>* r;
	Var<PosReal>* N0;
	Var<PosReal>* N1;
	Var<PosReal>* N2;
	double T0;

};



class CoalCalibratedChronogram : public CalibratedChronogram, public Rnode {

	public:

	CoalCalibratedChronogram(Tree* intree, Var<PosReal>* inrate, Var<PosReal>* ind, Var<PosReal>* inr, Var<PosReal>* inN0, Var<PosReal>* inN1, Var<PosReal>* inN2, double inT0, CalibrationSet* incalibset)	{

		DAGnode::SetName("Coal Chrono");

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		calibset = incalibset;

		scale = new Gamma(rate,rate);

		d = ind;
		r = inr;
		N0 = inN0;
		N1 = inN1;
		N2 = inN2;
		T0 = inT0;
		Ntaxa = GetTree()->GetSize(GetRoot());
		cerr << "Ntaxa : " << Ntaxa << '\n';

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());
		copytree = new CoalCopyTree(this,d,r,N0,N1,N2,T0);

		RecursiveRegisterProb(GetRoot());
		RecursiveSetCalibrations(GetRoot());

		RecursiveYoungerLimit(GetRoot());
		RecursiveOlderLimit(GetRoot(),-1);
		double age = RecursiveDrawAges(GetRoot());
		cerr << "ancestor age : " << age << '\n';

		RecursiveNormalizeTree(GetRoot(),T0,false);
		RecursiveUpdateBranches(GetRoot());

		scale->setval(T0);
		scale->Clamp();
		map<DAGnode*,int> tmpmap;
		scale->FullCorrupt(tmpmap);
		scale->FullUpdate();

		double tmp = RecursiveGetLogProb(GetRoot());
		if (fabs(tmp) > 1e-6)	{
			cerr << "error : nodes out of calib\n";
			exit(1);
		}
	}

	~CoalCalibratedChronogram()	{
		RecursiveDeleteBranch(GetRoot());
		RecursiveDeleteNode(GetRoot());
	}

	void Sample()	{
		CalibratedChronogram::Sample();
	}

	void drawSample()	{
		double age = RecursiveDrawAges(GetRoot());
		cerr << "in sample : recursive draw age : " << age << '\n';
		RecursiveNormalizeTree(GetRoot(),T0,false);
		RecursiveUpdateBranches(GetRoot());
		scale->Corrupt(true);
		scale->setval(T0);
		scale->Clamp();
		scale->Update();
		// CalibratedChronogram::Sample();
	}

	virtual double GetRootAge()	{
		return GetNodeVal(GetRoot()->GetNode())->val() * T0;
	}

	/*
	double Move(double tuning){
		return CalibratedChronogram::Move(tuning);
	}
	*/

	double Move(double tuning){
		int n = 0;
		double tmp = CoalMoveTime(tuning, GetRoot(),n);
		return tmp / n;
	}


	double RootMoveTime(double tuning){
		RootSetMin();
		double ret = GetNodeVal(GetRoot()->GetNode())->Move(tuning);
		return ret;
	}


	void RootSetMin(){
		Link* from = GetRoot();
		double min = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next()){
			if(GetAbsoluteTime(link->Out()) > min){
				min = GetAbsoluteTime(link->Out());
			}
		}
		GetNodeDate(from->GetNode())->SetMinMax(min, 3);
	}

	double CoalMoveTime(double tuning, Link* from, int& n){
		if(from->isLeaf()){
			return 0;
		}
		else{
			double retour = 0;
			for(Link* link=from->Next(); link!=from; link=link->Next()){
				retour += CoalMoveTime(tuning,link->Out(),n);
			}
			if(!from->isRoot()){
				retour += MoveTime(tuning, from);
				n++;
			}
			else	{
				retour += RootMoveTime(tuning);
				n++;
			}
			return retour;
		}
	}

	double GetLogProb()	{
		return Rnode::GetLogProb();
	}

	double ProposeMove(double)	{
		cerr << "error in CoalCalibratedChronogram: in propose move\n";
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
			cerr << "error : CoalChrono accepts only bifurcating trees\n";
			exit(1);
		}
	}


	double logProb()	{

		double rate = r->val() * T0;
		double A = (d->val() + r->val()) * T0 / N2->val() * (N2->val() - N1->val()) / N1->val() * exp(-rate);
		double B = d->val() * T0 / N2->val();
		double C = d->val() * T0 / N0->val();

		list<doublet> agelist;
		RecursiveGetAgeList(GetRoot(),agelist);

		agelist.sort();
		// agelist.sort(compare_doublet);
		// sort(agelist.begin(),agelist.end(),compare_doublet);

		double total = 0;
		int k = Ntaxa;

		double t1 = 0;
		for (list<doublet>::iterator i=agelist.begin(); i!= agelist.end(); i++)	{
			if (! k)	{
				cerr << "error in coal log prob\n";
				exit(1);
			}
			double t2 = i->first;
			double tmp = 0;
			if (t1 > 1)	{
				tmp = log(k * (k-1) * C) - k * (k-1) * C * (t2 - t1);
			}
			else if (t2 > 1)	{
				double tmp1 = - k * (k-1) * (B*(1 - t1) + A/rate*(exp(rate*1) - exp(rate*t1)));
				double tmp2 = log(k * (k-1) * C) - k * (k-1) * C * (t2 - 1);
				tmp = tmp1 + tmp2;
			}
			else	{
				tmp = log( k * (k-1) * (A * exp(rate * t2) + B)) - k * (k-1) * (B*(t2 - t1) + A/rate*(exp(rate*t2) - exp(rate*t1)));
			}
			total += tmp;
			// cerr << k << '\t' << rate << '\t' << t2 << '\t' << exp(rate * t2) << '\t' << tmp << '\n';
			t1 = t2;
			k--;
		}
		// cerr << '\n' << "total " << total << '\n';
		// exit(1);

		return total;
	}

	private:

	unsigned int Ntaxa;
	CoalCopyTree* copytree;
	Var<PosReal>* d;
	Var<PosReal>* r;
	Var<PosReal>* N0;
	Var<PosReal>* N1;
	Var<PosReal>* N2;
	double T0;
};

#endif

