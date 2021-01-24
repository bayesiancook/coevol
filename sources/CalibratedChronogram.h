
#ifndef CALIBCHRONO_H
#define CALIBCHRONO_H

#include "Chronogram.h"
#include "MultiVariateTreeProcess.h"


class Calibration	{

	public:

	Calibration() : taxon1(""), taxon2(""), older(-1), younger(-1) {}

	Calibration(string intax1, string intax2, double inolder, double inyounger) : taxon1(intax1), taxon2(intax2), older(inolder), younger(inyounger) {}

	~Calibration()	{}

	double GetOlderLimit() const {return older;}
	double GetYoungerLimit() const {return younger;}

	string GetTaxon1() const {return taxon1;}
	string GetTaxon2() const {return taxon2;}

	void ToStream(ostream& os) const {
		os << taxon1 << '\t' << taxon2 << '\t' << older << '\t' << younger << '\n';
	}

	void FromStream(istream& is)	{
		is >> taxon1 >> taxon2 >> older >> younger;
	}

	friend ostream& operator<<(ostream& os, const Calibration& cal)	{
		cal.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, Calibration& cal)	{
		cal.FromStream(is);
		return is;
	}

	private:
	// two names : taxon1 and taxon 2
	// a tree
	// two limits: older and younger
	// map : a pair of taxon -> two ages
	// map : a pair of taxon -> a node of the tree
	// function: a node of the tree -> two ages
	string taxon1;
	string taxon2;
	double older;
	double younger;

};

class ChronoScale : public Rvar<PosReal> {

	public:

	ChronoScale(Var<PosReal>* inroot, double inalpha, double inbeta, double inmin, double inmax, bool inactivelogprob = true)	{

		root = inroot;
		alpha = inalpha;
		beta = inbeta;
		min = inmin;
		max = inmax;
		if (min == -1)	{
			min = 0;
		}
		activelogprob = inactivelogprob;
		Register(root);
		Sample();
	}

	void drawSample()	{
		if (max == -1)	{
			double tmp = Random::Gamma(alpha,beta);
			setval(tmp);
		}
		else	{
			double tmp = Random::Uniform() * (max - min) + min;
			setval(tmp);
		}
	}

	double ProposeMove(double tuning)	{
		if (max == -1)	{
			/*
			cerr << "error in CalibratedNodeDate::ProposeMove: max is -1\n";
			exit(1);
			*/
			return PosReal::ProposeMove(tuning);
		}
		double tmp = val() + (max - min) * tuning * (Random::Uniform() - 0.5);
		int count = 0;
		while ((count < 1000) && ((max != -1) && (tmp > max)) || (tmp < min))	{
			if ((max != -1) && (tmp > max))	{
				tmp = 2*max - tmp;
			}
			if (tmp < min)	{
				tmp = 2*min - tmp;
			}
		}
		if (count == 1000)	{
			cerr << "error in chrono scale propose move\n";
			exit(1);
		}
		setval(tmp);
		return 0;
	}

	double logProb()	{
		if (! activelogprob)	{
			/*
			cerr << "log prob of root not active\n";
			cerr << "check soft bounds\n";
			exit(1);
			*/
			return 0;
		}
		if ((max != -1) && (val() > max))	{
			return log(0);
			// return -100000;
			/*
			cerr << "error in chrono scale: out of bound\n";
			cerr << min << '\t' << max << '\t' << val() << '\n';
			exit(1);
			*/
		}
		if (val() < min)	{
			return log(0);
			// return -100000;
			/*
			cerr << "error in chrono scale: out of bound\n";
			cerr << min << '\t' << max << '\t' << val() << '\n';
			exit(1);
			*/
		}

		if (alpha == -1)	{
			return  0;
		}
		return - beta * val()  + (alpha-1) * log(val());
	}

	protected:

	bool activelogprob;
	double alpha;
	double beta;
	Var<PosReal>* root;
	double min;
	double max;

};


class CalibrationSet	{

	public:

	CalibrationSet(Tree* intree) : tree(intree) {}

	bool IsCalibrated(const Node* node) const {
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		return (i != nodecalmap.end());
	}

	bool IsUpperCalibrated(const Node* node) const {
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		return ((i != nodecalmap.end())	&& (i->second.GetOlderLimit() != -1));
	}

	bool IsLowerCalibrated(const Node* node) const {
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		return ((i != nodecalmap.end())	&& (i->second.GetYoungerLimit() != -1));
	}

	double GetLowerCalibration(const Node* node) const {
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		if (i != nodecalmap.end())	{
			return i->second.GetYoungerLimit();
		}
		cerr << "error : get lower calibration called on node without lower calibration\n";
		exit(1);
		return -1;
	}

	double GetUpperCalibration(const Node* node) const {
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		if (i != nodecalmap.end())	{
			return i->second.GetOlderLimit();
		}
		cerr << "error : get upper calibration called on node without upper calibration\n";
		exit(1);
		return -1;
	}

	const Calibration& GetCalibration(const Node* node)	{
		map<const Node*,Calibration>::const_iterator i = nodecalmap.find(node);
		if (i == nodecalmap.end())	{
			cerr << "error : node does not belong to set of calibrated nodes\n";
			exit(1);
		}
		return i->second;
	}

	int GetNCalib()	{
		return nodecalmap.size();
	}

	void ToStream(ostream& os) const {
		int Ncalib = 0;
		for (map<const Node*,Calibration>::const_iterator i = nodecalmap.begin(); i!=nodecalmap.end(); i++)	{
			Ncalib++;
		}
		os << Ncalib << '\n';
		for (map<const Node*,Calibration>::const_iterator i = nodecalmap.begin(); i!=nodecalmap.end(); i++)	{
			os << i->second << '\n';
		}
	}

	void FromStream(istream& is)	{
		int Ncalib;
		is >> Ncalib;
		for (int i=0; i<Ncalib; i++)	{
			Calibration cal;
			is >> cal;
			const Link* link= tree->GetLCA(cal.GetTaxon1(),cal.GetTaxon2());
			if (! link)	{
				cerr << "error in calibration set: did not find common ancestor of " << cal << '\n';
				exit(1);
			}
			nodecalmap[link->GetNode()] = cal;
		}
	}

	private:

	map<const Node*,Calibration> nodecalmap;
	Tree* tree;
};


class FileCalibrationSet : public CalibrationSet	{

	public:

	FileCalibrationSet(string filename, Tree* intree) : CalibrationSet(intree) {
		ifstream is(filename.c_str());
		if (! is)	{
			cerr << "error : cannot open calibration file " << filename << "\n";
			exit(1);
		}
		FromStream(is);
	}
};

class CalibratedNodeDate : public NodeDate	{

	protected:

	double olderlimit;
	double youngerlimit;

	Var<PosReal>* scale;

	// CalibratedChronogram* chrono;

	public:

	CalibratedNodeDate(Var<PosReal>* inRate, Var<PosReal>* inScale) : NodeDate(inRate), olderlimit(-1), youngerlimit(-1), scale(inScale) {
		Register(scale);
		SetName("nodedate");
		// chrono = 0;
	}

	/*
	void SetChrono(CalibratedChronogram* inchrono)	{
		chrono = inchrono;
	}
	*/

	bool isCalibrated()	{
		return ((olderlimit != -1) || (youngerlimit != -1));
	}

	bool isUpperCalibrated()	{
		return (olderlimit != -1);
	}

	bool isLowerCalibrated()	{
		return (youngerlimit != -1);
	}

	void SetCalibrations(double inolderlimit, double inyoungerlimit)	{
		olderlimit = inolderlimit;
		youngerlimit = inyoungerlimit;
	}

	int CheckCalibration()	{
		if (((olderlimit != -1) && (value > olderlimit/scale->val())) || ((youngerlimit != -1) && (value < youngerlimit/scale->val())))	{
			cerr << "out of bounds : " << olderlimit << '\t' << youngerlimit << '\t' << value << '\t' << value * scale->val() << '\n';
			return 1;
		}
		return 0;
	}

	void SetMinMax(double inMin, double inMax){
		min = inMin;
		max = inMax;
		/*
		if (min < 0)	{
			cerr << "error in set min max: negative min\n";
			exit(1);
		}
		if (max < 0)	{
			cerr << "error in set min max: negative max\n";
			exit(1);
		}
		*/
		if ((olderlimit != -1) && (olderlimit/scale->val() < max))	{
			max = olderlimit/scale->val();
		}
		if ((youngerlimit != -1) && (youngerlimit/scale->val() > min))	{
			min = youngerlimit/scale->val();
		}
		if (max < min)	{
			cerr << "error : max < min : " << min << '\t' << max << '\n';
			cerr << inMin << '\t' << inMax << '\n';
			cerr << val() << '\n';
			cerr << youngerlimit << '\t' << olderlimit << '\t' << scale->val() << '\n';
			cerr << youngerlimit / scale->val() << '\t' << olderlimit / scale->val() << '\n';
			exit(1);
		}
	}

	double logProb()	{
		if ((olderlimit != -1) && (value > olderlimit/scale->val()))	{
			/*
			cerr << "date out of range\n";
			cerr << olderlimit << '\t' << youngerlimit << '\n';
			cerr << value * scale->val() << '\n';
			exit(1);
			*/
			return log(0);
			// return -100000;
		}
		if ((youngerlimit != -1) && (value < youngerlimit/scale->val()))	{
			/*
			cerr << "date out of range\n";
			cerr << olderlimit << '\t' << youngerlimit << '\n';
			cerr << value * scale->val() << '\n';
			exit(1);
			*/
			return log(0);
			// return -100000;
		}
		return 0;
	}

	// double Move(double tuning);

};

class CalibratedChronogram : public Chronogram	{

	public:

	CalibratedNodeDate* GetCalibratedNodeDate(const Node* node)	{
		CalibratedNodeDate* tmp = dynamic_cast<CalibratedNodeDate*>(GetNodeVal(node));
		if (! tmp)	{
			cerr << "error in Chronogram::GetNodeDate : dynamic cast\n";
			cerr << GetNodeVal(node) << '\n';
			cerr << GetNodeVal(node)->val() << '\n';
			exit(1);
		}
		return tmp;
	}

	bool isCalibrated(const Node* node) {
		return calibset->IsCalibrated(node);
		// return GetCalibratedNodeDate(node)->isCalibrated();
	}

	bool isUpperCalibrated(const Node* node) {
		return calibset->IsUpperCalibrated(node);
	}

	bool isLowerCalibrated(const Node* node) {
		return calibset->IsLowerCalibrated(node);
	}

	double GetLowerCalibration(const Node* node)	{
		return calibset->GetLowerCalibration(node);
	}

	double GetUpperCalibration(const Node* node)	{
		return calibset->GetUpperCalibration(node);
	}

	CalibratedChronogram()	{}

	CalibratedChronogram(Tree* intree, Var<PosReal>* inrate, double inalpha, double inbeta, CalibrationSet* incalibset, bool sample = true)	{

		SetWithRoot(false);
		tree = intree;
		rate = inrate;
		alpha = inalpha;
		beta = inbeta;
		calibset = incalibset;
		if (! calibset)	{
			calibset = new CalibrationSet(tree);
		}

		double rootmin = -1;
		double rootmax = -1;
		if (isCalibrated(GetRoot()->GetNode()))	{
			rootmin = GetLowerCalibration(GetRoot()->GetNode());
			if (rootmin == -1)	{
				rootmin = 0;
			}
			rootmax = GetUpperCalibration(GetRoot()->GetNode());
		}

		scale = new ChronoScale(rate,alpha,beta,rootmin,rootmax);
		// scale = new Gamma(alpha,beta);

		RecursiveCreateNode(GetRoot());
		RecursiveCreateBranch(GetRoot());

		// RecursiveSetChrono(GetRoot());

		if (!sample)	{

			RecursiveSetCalibrations(GetRoot());
			double maxage = RecursiveSetNodesValues(GetRoot());
			// RecursiveEqualizeLeafNodes(GetRoot(),maxage);
			RecursiveNormalizeTree(GetRoot(),maxage,true);
			RecursiveUpdateBranches(GetRoot());

			scale->setval(maxage);
			map<DAGnode*,int> tmpmap;
			scale->FullCorrupt(tmpmap);
			scale->FullUpdate();

		}
		else	{
			RecursiveSetCalibrations(GetRoot());

			RecursiveYoungerLimit(GetRoot());
			RecursiveOlderLimit(GetRoot(),-1);
			double d = RecursiveDrawAges(GetRoot());

			RecursiveNormalizeTree(GetRoot(),d,false);
			RecursiveUpdateBranches(GetRoot());

			scale->setval(d);
			map<DAGnode*,int> tmpmap;
			scale->FullCorrupt(tmpmap);
			scale->FullUpdate();

			double max = 0;
			for (int i=0; i<10000; i++)	{
				EmptyMove();
				double tmp = GetMinDeltaTime(GetRoot());
				if (max < tmp)	{
					max = tmp;
				}
			}
			max /= 2;
			while (GetMinDeltaTime(GetRoot()) < max)	{
				EmptyMove();
			}

			tmpmap.clear();
			scale->FullCorrupt(tmpmap);
			scale->FullUpdate();

			if (GetLogProb() != 0)	{
				cerr << "error in CalibratedChronogram: constructor: log prob : " << GetLogProb() << '\n';
				exit(1);
			}
		}
	}

    /*
	~CalibratedChronogram()	{
		RecursiveDeleteBranch(GetRoot());
		RecursiveDeleteNode(GetRoot());
	}
    */

	void drawSample()	{
		double d = RecursiveDrawAges(GetRoot());
		RecursiveNormalizeTree(GetRoot(),d,false);
		RecursiveUpdateBranches(GetRoot());
		scale->Corrupt(true);
		scale->setval(d);
		scale->Update();
	}

	/*
	void RecursiveSetChrono(const Link* from)	{

		GetCalibratedNodeDate(from->GetNode())->SetChrono(this);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetChrono(link->Out());
		}
	}
	*/

	virtual double GeneralMoveTime(double tuning, Link* from, int& n){
		if(from->isLeaf() && (! isCalibrated(from->GetNode()))){
			return 0;
		}
		else{
			double retour = 0;
			for(Link* link=from->Next(); link!=from; link=link->Next()){
				retour += GeneralMoveTime(tuning,link->Out(),n);
			}
			if(!from->isRoot()){
				retour += MoveTime(tuning, from);
				n++;
			}
			return retour;
		}
	}

	double ScaleEmptyMove()	{
		double total = 0;
		total += scale->Move(10);
		total += scale->Move(1);
		return total / 2;
	}

	ChronoScale* GetScale() {
		return scale;
	}

	virtual double GetRootAge()	{
		return scale->val();
	}

	double GetLogProb()	{
		return RecursiveGetLogProb(GetRoot());
	}

	virtual int CheckBounds()	{
		return RecursiveCheckBounds(GetRoot());
	}

	void Newick(ostream& os)	{

		NewickTreeToStream(os,GetRoot());
		os << ";\n";
	}

	int CheckLeafTimeConsistency(double factor)	{
		return RecursiveCheckLeafTimeConsistency(GetRoot(),factor,1.0);
	}

	int GetNfossil(double cutoff)	{
		return RecursiveGetNfossil(GetRoot(), cutoff);
	}

	int GetHowManyCross(double time)	{
		return RecursiveGetHowManyCross(GetRoot(),time);
	}

	protected:

	Rvar<PosReal>* CreateNodeVal (const Link* link){
		return new CalibratedNodeDate(rate,scale);
	}

	double RecursiveGetLogProb(const Link* from)	{
		double total = GetNodeVal(from->GetNode())->GetLogProb();
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogProb(link->Out());
		}
		return total;
	}

	void RecursiveSetCalibrations(const Link* from)	{
		if (calibset->IsCalibrated(from->GetNode()))	{
			double older = calibset->GetCalibration(from->GetNode()).GetOlderLimit();
			double younger = calibset->GetCalibration(from->GetNode()).GetYoungerLimit();
			GetCalibratedNodeDate(from->GetNode())->SetCalibrations(older,younger);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetCalibrations(link->Out());
		}
	}

	double RecursiveYoungerLimit(const Link* from)	{
		if (from->isLeaf())	{
			youngerlimit[from] = 0;
		}
		else	{
			double max = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				double tmp = RecursiveYoungerLimit(link->Out());
				if (max < tmp)	{
					max = tmp;
				}
			}
			if (calibset->IsCalibrated(from->GetNode()))	{
				double tmp = calibset->GetCalibration(from->GetNode()).GetYoungerLimit();
				if ((tmp != -1) && (max < tmp))	{
					max = tmp;
				}
			}
			youngerlimit[from] = max;
		}
		return youngerlimit[from];
	}

	void RecursiveOlderLimit(const Link* from, double min)	{
		if (from->isLeaf())	{
			olderlimit[from] = 0;
		}
		else	{
			if (from->isRoot())	{
				if (calibset->IsCalibrated(from->GetNode()))	{
					min = calibset->GetCalibration(from->GetNode()).GetOlderLimit();
				}
				if (min == -1)	{
					if (youngerlimit[from] == -1)	{
						cerr << "error in set calibrations : root is not constrained\n";
						exit(1);
					}
					min = 2 * youngerlimit[from];
				}
			}
			else	{
				if (min == -1)	{
					cerr << "error in set calib: " << min << '\n';
					exit(1);
				}
				if (calibset->IsCalibrated(from->GetNode()))	{
					double tmp = calibset->GetCalibration(from->GetNode()).GetOlderLimit();
					if (tmp != -1)	{
						if (min > tmp)	{
							min = tmp;
						}
					}
				}
			}
			olderlimit[from] = min;
			if (olderlimit[from] < youngerlimit[from])	{
				cerr << "error in set calibrations : " << olderlimit[from] << '\t' << youngerlimit[from] << '\n';
				exit(1);
			}
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				RecursiveOlderLimit(link->Out(),min);
			}
		}
	}

	double RecursiveDrawAges(const Link* from)	{
		double age = 0;
		if (!from->isLeaf())	{
			double min = 0;
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				double tmp = RecursiveDrawAges(link->Out());
				if (min < tmp)	{
					min = tmp;
				}
			}
			if (min < youngerlimit[from])	{
				min = youngerlimit[from];
			}
			if (min < 0)	{
				cerr << "error in recursive draw ages : min : " << min << '\n';
				exit(1);
			}
			if (olderlimit[from] < min)	{
				cerr << "error in recursive draw ages : min : " << min << " and max : " << olderlimit[from] << '\n';
				exit(1);
			}
			age = (olderlimit[from] - min) * Random::Uniform() + min;
		}

		GetNodeVal(from->GetNode())->setval(age);
		return age;
	}

	int RecursiveCheckBounds(const Link* from)	{
		int tot = 0;
		if (calibset->IsCalibrated(from->GetNode()))	{
			int tmp = GetCalibratedNodeDate(from->GetNode())->CheckCalibration();
			if (tmp)	{
				if (from->isRoot())	{
					cerr << "root out of bounds\n";
					exit(1);
				}
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveCheckBounds(link->Out());
		}
		return tot;
	}

	int RecursiveGetNfossil(const Link* from, double cutoff)	{

		int tot = 0;
		if (from->isLeaf() && (GetNodeVal(from->GetNode())->val() * scale->val() > cutoff))	{
			tot++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			tot += RecursiveGetNfossil(link->Out(),cutoff);
		}
		return tot;
	}

	int RecursiveCheckLeafTimeConsistency(const Link* from, double factor, double previoustime)	{

		int ret = 1;
		if (from->isLeaf())	{
			ret *= (GetNodeDate(from->GetNode())->val() * factor < previoustime);
		}
		else	{
			double time = GetNodeDate(from->GetNode())->val();
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				ret *= RecursiveCheckLeafTimeConsistency(link->Out(),factor,time);
			}
		}
		return ret;
	}

	int RecursiveGetHowManyCross(const Link* from, double time)	{

		int n = 0;
		if (! from->isRoot())	{
			double mintime = GetNodeDate(from->GetNode())->val() * scale->val();
			double maxtime = GetNodeDate(from->Out()->GetNode())->val() * scale->val();
			if ((mintime < time) && (maxtime > time))	{
				n++;
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			n += RecursiveGetHowManyCross(link->Out(),time);
		}
		return n;
	}

	string GetNodeName(const Link* link)	{

		ostringstream s;
		s << link->GetNode()->GetName();
		s << "_";
		// s << GetRootAge() * ;
		s << GetNodeDate(link->GetNode())->val();
		return s.str();
	}

	string GetBranchName(const Link* link)	{

		ostringstream s;
		// s << GetRootAge() * ;
		s << GetBranchLength(link->GetBranch())->val();
		return s.str();
	}

	void NewickTreeToStream(ostream& os, const Link* from) {

		if (!from->isLeaf())	{
			os << '(';
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				NewickTreeToStream(os,link->Out());
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
		// if (from->isLeaf())	{
			os << GetNodeName(from);
		// }
		if (!from->isRoot())	{
			string brval = GetBranchName(from);
			if (brval != "")	{
				os << ':' << brval;
			}
		}
	}

	protected:

	map<const Link*,double> olderlimit;
	map<const Link*,double> youngerlimit;

	CalibrationSet* calibset;
	ChronoScale* scale;

	double alpha;
	double beta;
};


class AllInternalNodesMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;

	public:

	AllInternalNodesMove(CalibratedChronogram* inchrono, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
	}

	double Move(double tuning_modulator){
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = 1.0 / e;

		int consistent = chrono->CheckLeafTimeConsistency(f);
		if (!consistent)	{
			return 0;
		}

		/*
		double logprobbefore = chrono->GetLogProb();
		double scalelogprobbefore = chrono->GetScale()->GetLogProb();
		*/

		Corrupt(true);
		chrono->GetScale()->ScalarMultiplication(e);
		chrono->MultiplyLeafTimes(chrono->GetRoot(),f);
		double logratio = Update();

		/*
		double logprobafter = chrono->GetLogProb();
		double scalelogprobafter = chrono->GetScale()->GetLogProb();

		double delta = logprobafter - logprobbefore;
		double scaledelta = scalelogprobafter - scalelogprobbefore;
		cerr << logratio << '\t' << delta << '\t' << scaledelta << '\t' << delta + scaledelta << '\n';
		*/

		// log hastings
		logratio += (chrono->GetTree()->GetSize() - 1) * m;
		// logratio += (1 - chrono->GetNfossil(1.0)) * m;

		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
			/*
			double logprobrestore = chrono->GetLogProb();
			double scalelogprobrestore = chrono->GetScale()->GetLogProb();
			if (fabs(logprobbefore - logprobrestore) > 1e-8)	{
				cerr << "restore failure\n";
				exit(1);
			}
			if (fabs(scalelogprobbefore - scalelogprobrestore) > 1e-8)	{
				cerr << "restore failure\n";
				exit(1);
			}
			*/
		}
		return ((double) accepted);
	}
};

class ChronoRootOnlyMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;

	public:

	ChronoRootOnlyMove(CalibratedChronogram* inchrono, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
	}

	double Move(double tuning_modulator){

		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = 1.0 / e;
		double age = chrono->GetScale()->val();
		double left = chrono->GetNodeVal(chrono->GetRoot()->Next()->Out()->GetNode())->val() * age;
		double right = chrono->GetNodeVal(chrono->GetRoot()->Next()->Next()->Out()->GetNode())->val() * age;

		double min = (left < right) ? right : left;

		if (age * e < min)	{
			return 0;
		}

		Corrupt(true);

		age *= e;
		chrono->MultiplyInternalTimes(chrono->GetRoot(),f);
		chrono->GetScale()->setval(age);

		double logratio = Update();

		logratio += m;
		// logratio += (1 - (chrono->GetTree()->GetSize() - 2 + chrono->GetNfossil(1.0))) * m;

		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return ((double) accepted);
	}
};

class RootMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;

	public:

	RootMove(CalibratedChronogram* inchrono, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
	}

	double Move(double tuning_modulator){
		cerr << "in root move: incorrect log prob (should call chrono log prob)\n";
		exit(1);
		Corrupt(true);
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = exp(-m);
		// double checkupdate = - chrono->GetScale()->logProb() - chrono->CheckLogProb();
		chrono->GetScale()->ScalarMultiplication(e);
		chrono->MultiplyTimes(chrono->GetRoot(),f);
		chrono->GetNodeDate(chrono->GetRoot()->GetNode())->setval(1.0);
		double logratio = Update();
		// checkupdate += chrono->GetScale()->logProb() + chrono->CheckLogProb();
		cerr << "delta log prob : " << logratio << '\n';
		// cerr << "chrono log prob ; " << chrono->CheckLogProb() << '\n';
		logratio -= (chrono->GetNinternalNode() - 2) * m;
		cerr << "hastings : " << (chrono->GetNinternalNode() - 2) * m << '\n';
		cerr << "final log ratio : " << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		cerr << "accepted : " << accepted << '\n';
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return ((double) accepted);
	}
};

class RootProcessMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;
	MultiVariateTreeProcess* process;
	int index;

	public:

	RootProcessMove(CalibratedChronogram* inchrono, MultiVariateTreeProcess* inprocess, int inindex, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		process = inprocess;
		index = inindex;
		process->RecursiveRegister(this,process->GetRoot());
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = exp(-m);
		chrono->GetScale()->ScalarMultiplication(e);
		chrono->MultiplyTimes(chrono->GetRoot(),f);
		chrono->GetNodeDate(chrono->GetRoot()->GetNode())->setval(1.0);
		process->RootNeighborPiecewiseTranslation(m,index,1);
		double logratio = Update();
		logratio -= (chrono->GetNinternalNode() - 2) * m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return ((double) accepted);
	}
};

	/*
	double CalibratedNodeDate::Move(double tuning)	{

		cerr << "cal node date move\n";
		cerr << chrono << '\n';
		if (! isClamped())	{
			cerr << "before value  " << val() << '\n';
			cerr << "before log prob \n";
			double before = chrono->GetLogProb();
			Corrupt(true);
			double logHastings = ProposeMove(tuning);
			cerr << "after value  " << val() << '\n';
			cerr << "update\n";
			double deltaLogProb = Update();
			double logRatio = deltaLogProb + logHastings;
			cerr << "after log prob\n";
			double after = chrono->GetLogProb();

			cerr << before << '\t' << after << '\t' << after - before << '\t' << deltaLogProb << '\t' << logHastings << '\n';

			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
			}
			return (double) accepted;
		}
		return 1;
	}
	*/

#endif
