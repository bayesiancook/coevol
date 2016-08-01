
#ifndef TIMELINE_H
#define TIMELINE_H

#include "Chronogram.h"

class TimeIntervals : public ValPtrArray<Dvar<PosReal> >	{

	public:

	TimeIntervals(Chronogram* inchronogram) : ValPtrArray<Dvar<PosReal> > (inchronogram->GetTree()->GetSize()-1) {
		chronogram = inchronogram;
		datelist = chronogram->GetDateList();
		Create();
	}

	Var<PosReal>* GetValFromLink(const Link* link)	{
		if (link->isLeaf())	{
			return GetVal(GetSize() - 1);
		}
		if (link->isRoot())	{
			return 0;
		}
		Var<PosReal>* date = chronogram->GetNodeDate(link->GetNode());
		if (date2bl.find(date) == date2bl.end())	{
			cerr << "error in time intervals: did not find date in hash table\n";
			exit(1);
		}
		return date2bl[date];
	}

	double GetDate(int i)	{
		return datelist[i]->val();
	}

	protected:

	Dvar<PosReal>* CreateVal(int i)	{
		double tmp = ((double) datelist[i]->val()) - ((double) datelist[i+1]->val());
		if (tmp < 0)	{
			cerr << "error in time intervals: negative duration\n";
			exit(1);
		}
		Dvar<PosReal>* bl = new BranchLength(datelist[i], datelist[i+1], chronogram->GetRate());
		if (date2bl.find(datelist[i]) != date2bl.end())	{
			cerr << "error in time intervals: twice in hash table\n";
			exit(1);
		}
		date2bl[datelist[i]] = bl;
		return bl;
	}

	private:

	Chronogram* chronogram;
	vector<Var<PosReal>*> datelist;
	map<Var<PosReal>*, Var<PosReal>* > date2bl;
};

class TimeLine : public  IIDArray<RealVector>	{

	public:

	TimeLine(Chronogram* inchronogram, TimeIntervals* intimeintervals, Var<CovMatrix>* insigma) : IIDArray<RealVector>(inchronogram->GetTree()->GetSize()) {
		chronogram = inchronogram;
		timeintervals = intimeintervals;
		sigma = insigma;
		Create();
	}

	Rvar<RealVector>* GetValFromLink(const Link* link)	{
		Var<PosReal>* bl = timeintervals->GetValFromLink(link);
		if (! bl)	{
			return GetVal(0);
		}
		if (bl2val.find(bl) == bl2val.end())	{
			cerr << "error in timeline: did not find entry in hashtable\n";
			exit(1);
		}
		return bl2val[bl];
	}

	protected:

	Rvar<RealVector>* CreateVal(int i)	{
		MultiNormal* tmp = 0;
		if (! i)	{
			tmp = new MultiNormal(sigma);
			tmp->ClampAtZero();
			if (bl2val.find(0) != bl2val.end())	{
				cerr << "error in tineline create val: root\n";
				exit(1);
			}
			bl2val[0] = tmp;
		}
		else	{
			Var<PosReal>* bl = timeintervals->GetVal(i-1);
			tmp = new MultiNormal(sigma,GetVal(i-1),bl);
			if (bl2val.find(bl) != bl2val.end())	{
				cerr << "error in time line: twice in hash table\n";
				exit(1);
			}
			bl2val[bl] = tmp;
		}
		return tmp;
	}

	private:

	Chronogram* chronogram;
	TimeIntervals* timeintervals;
	Var<CovMatrix>* sigma;
	map<Var<PosReal>*, Rvar<RealVector>* > bl2val;
};

#endif
