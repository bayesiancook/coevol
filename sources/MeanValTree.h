
#ifndef MEANVALTREE
#define MEANVALTREE

#include "ValTree.h"
#include "Var.h"
#include "MultiVariateTreeProcess.h"
#include <cmath>
#include <list>

class MeanBranchTree : public NewickTree {

	public:

	MeanBranchTree(Tree* intree, bool inwithroot = false) : tree(intree), withRoot(inwithroot), size(0)	{
        Reset();
	}

	const Tree* GetTree() const {return tree;}

	const Link* GetRoot() const {return GetTree()->GetRoot();}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Add(BranchValPtrTree<Rvar<PosReal> >* sample)	{
		RecursiveAdd(sample,GetRoot());
		size++;
	}

	void Add(BranchValPtrTree<Dvar<PosReal> >* sample)	{
		RecursiveAdd(sample,GetRoot());
		size++;
	}

	void Add(BranchValPtrTree<Rvar<UnitReal> >* sample)	{
		RecursiveAdd(sample,GetRoot());
		size++;
	}

	void Add(BranchValPtrTree<Dvar<UnitReal> >* sample)	{
		RecursiveAdd(sample,GetRoot());
		size++;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	double GetMean(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = mean.find(branch);
		return i->second;
	}

	double GetVar(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = var.find(branch);
		return i->second;
	}

	string GetNodeName(const Link* link) const {
		return link->GetNode()->GetName();
	}

	string GetBranchName(const Link* link) const {
		ostringstream s;
		s << GetMean(link->GetBranch());
		return s.str();
	}

	private:

	void RecursiveAdd(BranchValPtrTree<Rvar<PosReal> >* sample, const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, link->Out());
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
		}
		if (withRoot)	{
			if (! sample->WithRoot())	{
				cerr << "error in mean branch tree: sample has no root value\n";
				throw;
			}
			double tmp = sample->GetBranchVal(0)->val();
			mean[0] += tmp;
			var[0] += tmp * tmp;
		}
	}

	void RecursiveAdd(BranchValPtrTree<Dvar<PosReal> >* sample, const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, link->Out());
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
		}
		if (withRoot)	{
			if (! sample->WithRoot())	{
				cerr << "error in mean branch tree: sample has no root value\n";
				throw;
			}
			double tmp = sample->GetBranchVal(0)->val();
			mean[0] += tmp;
			var[0] += tmp * tmp;
		}
	}

	void RecursiveAdd(BranchValPtrTree<Rvar<UnitReal> >* sample, const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, link->Out());
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
		}
		if (withRoot)	{
			if (! sample->WithRoot())	{
				cerr << "error in mean branch tree: sample has no root value\n";
				throw;
			}
			double tmp = sample->GetBranchVal(0)->val();
			mean[0] += tmp;
			var[0] += tmp * tmp;
		}
	}

	void RecursiveAdd(BranchValPtrTree<Dvar<UnitReal> >* sample, const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, link->Out());
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
		}
		if (withRoot)	{
			if (! sample->WithRoot())	{
				cerr << "error in mean branch tree: sample has no root value\n";
				throw;
			}
			double tmp = sample->GetBranchVal(0)->val();
			mean[0] += tmp;
			var[0] += tmp * tmp;
		}
	}

	void RecursiveReset(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			mean[link->GetBranch()] = 0;
			var[link->GetBranch()] = 0;
		}
		if (withRoot)	{
			mean[0] = 0;
			var[0] = 0;
		}
	}

	void RecursiveNormalise(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			mean[link->GetBranch()] /= size;
			var[link->GetBranch()] /= size;
			var[link->GetBranch()] -= mean[link->GetBranch()] * mean[link->GetBranch()];
		}
	}

	map<const Branch*,double> mean;
	map<const Branch*,double> var;
	Tree* tree;
	bool withRoot;
	int size;

};

class MeanExpNormTree : public NewickTree {

	public:

	MeanExpNormTree(Tree* intree, bool inlogit, bool inprintlog, bool inprintmean, bool inprintmed, bool inprintci, bool inprintstdev, bool inwithleaf, bool inwithinternal, double inmeanreg = 0, double instdevreg = 0) : tree(intree), logit(inlogit), printlog(inprintlog), printmean(inprintmean), printmed(inprintmed), printci(inprintci), printstdev(inprintstdev), withleaf(inwithleaf), withinternal(inwithinternal), meanreg(inmeanreg), stdevreg(instdevreg) {
		ppleafroot = 0;
		threshold = 0;
		withpp = false;
		withdepth = false;

        logscale = 1.0;

        if (meanreg)    {
            cerr << "in MeanExpNormTree: check interplay meanreg, _GetMeanLog, _GetVarLog and log scale\n";
            exit(1);
        }

		Reset();
	}

	void ActivatePP(double inthreshold)	{
		threshold = inthreshold;
		withpp = true;
		Reset();
	}

	void SetLog() {printlog = true;}

    void SetLogScale(double inscale)    {
        logscale = log(inscale);
    }

	void SetPrintMean(bool inprintmean = true) {printmean = inprintmean;}

    void SetPrintMedian(bool inprintmed = true) {printmed = inprintmed;}

	void SetPrintCI(bool inprintci = true) {
		printci = inprintci;
		if (printci && printstdev)	{
			printstdev = false;
		}
	}

	void SetPrintStdev(bool inprintstdev = true) {
		printstdev = inprintstdev;
		if (printci && printstdev)	{
			printci = false;
		}
	}

	void SetWithLeaf(bool inwithleaf)	{
		withleaf = inwithleaf;
	}

	void SetWithInternal(bool inwithinternal)	{
		withinternal = inwithinternal;
	}

	bool isFixed(const Node* node)	const {
		if (printlog)	{
			return (fabs(GetVar(node)) < 1e-5);
		}
		return (sqrt(fabs(GetVar(node))) / GetMean(node) < 1e-5);
	}

	const Tree* GetTree() const {return tree;}

	const Link* GetRoot() const {return GetTree()->GetRoot();}

	void CutoffFromBelow(double cutoff)	{
		RecursiveCutoffFromBelow(GetRoot(),cutoff);
	}

	void SetWithDepth(LengthTree* inlengthtree)	{
		withdepth = true;
		lengthtree = inlengthtree;
	}

	double GetDepth(const Link* from)	{
		if (from->isLeaf())	{
			return 0;
		}
		double mean = 0;
		int n = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			mean += GetDepth(link->Out());
			mean += lengthtree->GetBranchVal(link->GetBranch())->val();
			n++;
		}
		return mean / n;
	}

	double GetMin95(const Node* node) const {
		if (meanreg)	{
			double tmp = _GetMeanLog(node) - 1.959964 * sqrt(_GetVarLog(node));
			if (printlog)	{
				return tmp / logscale;
			}
			else if (logit)	{
				return exp(tmp) / (1 + exp(tmp));
			}
			return exp(tmp);
		}
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 2.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return printlog ? *i / logscale : (logit? exp(*i) / (1 + exp(*i)) : exp(*i));
	}

	double GetMax95(const Node* node) const {
		if (meanreg)	{
			double tmp = _GetMeanLog(node) + 1.959964 * sqrt(_GetVarLog(node));
			if (printlog)	{
				return tmp;
			}
			else if (logit)	{
				return exp(tmp) / (1 + exp(tmp));
			}
			return exp(tmp);
		}
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 97.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return printlog ? *i / logscale : (logit ? exp(*i) / (1 + exp(*i)) : exp(*i));
	}

	double _GetMean(const Node* node) const	{
		map<const Node*, double>::const_iterator i = mean.find(node);
		return i->second;
	}

	double _GetMedian(const Node* node) const	{
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 50));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return exp(*i);
	}
	
	double _GetVar(const Node* node) const	{
		map<const Node*, double>::const_iterator i = var.find(node);
		return i->second;
	}

	double _GetMeanLog(const Node* node) const	{
		map<const Node*, double>::const_iterator i = meanlog.find(node);
		if (meanreg)	{
			return meanreg * i->second + stdevreg;
		}
		return i->second / logscale;
	}

	double _GetVarLog(const Node* node) const	{
		map<const Node*, double>::const_iterator i = varlog.find(node);
		map<const Node*, double>::const_iterator j = meanlog.find(node);
		if (meanreg)	{
			return meanreg * meanreg * i->second;
		}
		return i->second / logscale / logscale;
	}

	double GetMeanTime(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = meantime.find(branch);
		return i->second;
	}

	double GetVarTime(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = vartime.find(branch);
		return i->second;
	}

	double GetMean(const Node* node) const {
		return printlog ? _GetMeanLog(node) : _GetMean(node);
	}

	double GetMedian(const Node* node) const {
		return printlog ? log(_GetMedian(node)) / logscale : _GetMedian(node);
	}
	
	double GetVar(const Node* node) const {
		return printlog ? _GetVarLog(node) : _GetVar(node);
	}

	string GetNodeName(const Link* link) const {
		ostringstream s;
		bool empty = true;
		if (link->isLeaf())	{
			s << link->GetNode()->GetName();
			empty = false;
		}
		if (printmean)	{
			if (! empty)	{
				s << '_';
			}
			s << GetMean(link->GetNode());
			empty = false;
		}
		if (printmed)	{
			if (! empty)	{
				s << '_';
			}
			s << GetMedian(link->GetNode());
			empty = false;
		}
		if (((! printmean) && (! printmed)) || (! isFixed(link->GetNode())))	{
			if (printstdev)	{
				if (! empty)	{
					s << '_';
				}
				s << sqrt(GetVar(link->GetNode()));
				empty = false;
			}
			if (printci)	{
				if (! empty)	{
					s << '_';
				}
				// s.precision(2); s.setf ( ios::floatfield|ios::showpoint|ios::fixed );
				// s << ((double) ((int) (100 * GetMin95(link->GetNode())))) / 100;
				s << GetMin95(link->GetNode());
				s << '_';
				s << GetMax95(link->GetNode());
				// s << ((double) ((int) (100 * GetMax95(link->GetNode())))) / 100;
				empty = false;
			}
		}
		return s.str();
	}

	string GetBranchName(const Link* link) const {
		ostringstream s;
		s << GetMeanTime(link->GetBranch());
		return s.str();
	}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		ppleafroot /= size;
		RecursiveNormalise(GetTree()->GetRoot());
	}

	double GetPPLeafRoot()	{
		return ppleafroot;
	}

	void AddLogLinearCombination(NodeVarTree<RealVector>* sample, LengthTree* chronogram, const vector<double>& alpha, double beta)   {
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		RecursiveAddLogLinearCombination(sample, chronogram, GetTree()->GetRoot(), alpha, beta);
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}

    /*
	void Add(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index, Var<Real>* Offset = 0)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		double offset = 0;
		if (Offset)	{
			offset = Offset->val();
		}
		RecursiveAdd(sample, chronogram, index, offset, GetTree()->GetRoot());
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}
    */

	void Add(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index, double offset = 0)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		RecursiveAdd(sample, chronogram, index, offset, GetTree()->GetRoot());
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}

	void Add(NodeVarTree<RealVector>* sample, LengthTree* chronogram, double* slopes, double offset)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		RecursiveAdd(sample, chronogram, slopes, offset, GetTree()->GetRoot());
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}

	void AddGC(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index)	{
		RecursiveAddGC(sample, chronogram, index, GetTree()->GetRoot());
		size++;
	}

	void Add(NodeVarTree<Real>* sample, LengthTree* chronogram, double offset = 0)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot(), offset);
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}

	void Add(BranchVarTree<PosReal>* sample, LengthTree* chronogram, bool mean = false)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot(), mean);
		size++;
	}

	void Add(BranchVarTree<UnitReal>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void LooTabulate(ostream& os, ContinuousData* contdata)	{
		RecursiveLooTabulate(os,contdata,GetTree()->GetRoot());
	}

    void TabulateHeader(ostream& os)    {
        os << "#taxon1\ttaxon2";
        if (withdepth)	{
            os << "\tdepth";
        }
        else	{
            os << "\tbranch_delta_t";
        }
        if (printmean)	{
            os << "\tmean";
        }
        if (printmed)   {
            os << "\tmedian";
        }
        if (printstdev)	{
            os << "\tstdev";
        }
        if (printci)	{
            os <<"\tmin95\tmax95";
        }
        os << '\n';
    }

	void Tabulate(ostream& os)	{
        TabulateHeader(os);
		RecursiveTabulate(os,GetTree()->GetRoot());
	}

	void TabulatePP(ostream& os)	{
		RecursiveTabulatePP(os,GetTree()->GetRoot());
	}

	void TabulateDistribution(ostream& os)	{
		RecursiveTabulateDistribution(os,GetTree()->GetRoot());
	}

	void CheckRootToTip()	{
		cerr << "check root to tip\n";
		CheckRootToTip(GetRoot());
		cerr << '\n';
		cerr << "ok\n";
	}

	void CheckRootToTip(const Link* from, double p = 0)	{
		if (from->isLeaf())	{
			cerr << p  << '\t';
		}
		else	{
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				CheckRootToTip(link->Out(), p + GetMeanTime(link->GetBranch()));
			}
		}
	}

	private:

	void RecursiveCutoffFromBelow(const Link* from, double cutoff)	{

		if (mean[from->GetNode()] < cutoff)	{
			mean[from->GetNode()] = 0;
		}
		
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCutoffFromBelow(link->Out(),cutoff);
		}
	}

	void RecursiveTabulatePP(ostream& os, Link* from)	{
		if ((from->isLeaf() && withleaf) || ((! from->isLeaf()) && (withinternal)))	{
			if (pp[from->GetBranch()] > 0.95)	{
				os << '*' << '\t';
			}
			else	{
				os << '.' << '\t';
			}
			os << ((double) ((int) (1000 * pp[from->GetBranch()]))) / 10.0 << '\t';
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulatePP(os,link->Out());
		}
	}

	void RecursiveTabulate(ostream& os, Link* from)	{
		if ((from->isLeaf() && withleaf) || ((! from->isLeaf()) && (withinternal)))	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from);
			if (withdepth)	{
				os << '\t' << GetDepth(from);
			}
			else	{
				if (from->GetBranch())	{
					os << '\t' << GetMeanTime(from->GetBranch());
				}
				else	{
					os << '\t' << 0;
				}
			}
			if (printmean)	{
				os << '\t' << GetMean(from->GetNode());
			}
            if (printmed)   {
                os << '\t' << GetMedian(from->GetNode());
            }
			if (printstdev)	{
				if (! isFixed(from->GetNode()))	{
					os << '\t' << sqrt(GetVar(from->GetNode()));
				}
				else	{
					os << '\t' << 0;
				}
			}
			if (printci)	{
				os << '\t' << GetMin95(from->GetNode()) << '\t' << GetMax95(from->GetNode());
			}
			os << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(os,link->Out());
		}
	}

	void RecursiveTabulateDistribution(ostream& os, Link* from)	{
		/*
		if (fabs(var[from->GetNode()]) > 1e-6)	{
			os << meanlog[from->GetNode()] / log(10.0)  << '\t' << sqrt(varlog[from->GetNode()]) / log(10.0)  << '\t';
		}
		else	{
			os << meanlog[from->GetNode()] / log(10.0) << '\t' << 0 << '\t';
		}
		*/

		if ((from->isLeaf() && withleaf) || ((! from->isLeaf()) && (withinternal)))	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			// os << meanlog[from->GetNode()] << '\t';
			map<const Node*, list<double> >::const_iterator f = dist.find(from->GetNode());
			list<double> l = f->second;
			for (list<double>::const_iterator i = l.begin(); i != l.end(); i++)	{
				if (printlog)	{
					os << *i / log(10.0) << '\t';
				}
				else	{
					os << exp(*i) << '\t';
				}
			}
			os << '\n';
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulateDistribution(os,link->Out());
		}
	}

	void RecursiveLooTabulate(ostream& os, ContinuousData* contdata, Link* from)	{
		if (from->isLeaf())	{
			if (contdata->isMissing(from->GetNode()->GetName()))	{
			// if (varlog[from->GetNode()] > 1e-3)	{
				os << meanlog[from->GetNode()]  << '\t' << sqrt(varlog[from->GetNode()])  << '\n';
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveLooTabulate(os,contdata,link->Out());
		}
	}

	void RecursiveAddGC(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index, Link* from)	{
		double tmp1 = exp((* sample->GetNodeVal(from->GetNode()))[index+1]) + exp((* sample->GetNodeVal(from->GetNode()))[index+2]);
		double tmp2 = exp((* sample->GetNodeVal(from->GetNode()))[index+0]) + exp((* sample->GetNodeVal(from->GetNode()))[index+3]);
		double temp = tmp1 / (tmp2 + tmp1);
		double tmp = log(temp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddGC(sample, chronogram, index, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(NodeVarTree<RealVector>* sample, LengthTree* chronogram, double* slopes, double offset, Link* from)	{
		double* nodeval = sample->GetNodeVal(from->GetNode())->GetArray();
		int dim = sample->GetNodeVal(from->GetNode())->GetDim();
		// compute linear combination, based on nodeval, slopes and offset
		// store it into tmp
		double tmp = offset;
		for (int i=0; i<dim; i++)	{
			tmp += slopes[i] * nodeval[i];
		}

		if (from->isRoot())	{
			meanroot = tmp;
		}
		else if (from->isLeaf())	{
			meanleaf += tmp;
			leafsize++;
		}
		double temp = logit ? exp(tmp) / (1 + exp(tmp)) : exp(tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
		dist[from->GetNode()].push_front(tmp);

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, slopes, offset, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAddLogLinearCombination(NodeVarTree<RealVector>* sample, LengthTree* chronogram, Link* from, const vector<double>& alpha, double beta)  {

        double tmp = beta;
        for (unsigned int i=0; i<alpha.size(); i++) {
            tmp += alpha[i] * (* sample->GetNodeVal(from->GetNode()))[i];
        }

		if (from->isRoot())	{
			meanroot = tmp;
		}
		else if (from->isLeaf())	{
			meanleaf += tmp;
			leafsize++;
		}

		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		dist[from->GetNode()].push_front(tmp);

        double exptmp = exp(tmp);
		double temp = logit ? exptmp / (1.0 + exptmp) : exptmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddLogLinearCombination(sample, chronogram, link->Out(), alpha, beta);
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}
	
	void RecursiveAdd(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index, double offset, Link* from)	{
		double tmp = (* sample->GetNodeVal(from->GetNode()))[index];
		tmp += offset;
		if (from->isRoot())	{
			meanroot = tmp;
		}
		else if (from->isLeaf())	{
			meanleaf += tmp;
			leafsize++;
		}
		double temp = logit ? exp(tmp) / (1 + exp(tmp)) : exp(tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
		dist[from->GetNode()].push_front(tmp);

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, index, offset, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(NodeVarTree<Real>* sample, LengthTree* chronogram, Link* from, double offset)	{
		double tmp = sample->GetNodeVal(from->GetNode())->val();
		tmp += offset;
		if (from->isRoot())	{
			meanroot = tmp;
		}
		else if (from->isLeaf())	{
			meanleaf += tmp;
			leafsize++;
		}
		double temp = logit ? exp(tmp) / (1 + exp(tmp)) : exp(tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;

		dist[from->GetNode()].push_front(tmp);

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out(), offset);
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(BranchVarTree<PosReal>* sample, LengthTree* chronogram, Link* from, bool branchmean)	{
		// double temp = sample->GetBranchVal(from->GetBranch())->val();
		double temp = (from->isRoot()) ? 1.0 : ((double) sample->GetBranchVal(from->GetBranch())->val());
		double length = 1.0;
		if (! from->isRoot())	{
			length = chronogram->GetBranchVal(from->GetBranch())->val();
		}
		if (branchmean && (!length))	{
			// cerr << "null branch length\n";
			length = 1e-6;
		}
		if (branchmean)	{
			if (!from->isRoot())	{
				temp /= length;
			}
			else	{
				temp = 1;
			}
		}
		if (std::isnan(temp))	{
			cerr << "error in recursive add: nan\n";
			if (from->isRoot())	{
				cerr << "root\n";
				cerr << "integral: " << temp << '\n';
			}
			else	{
				cerr << "branch length : " << chronogram->GetBranchVal(from->GetBranch())->val() << '\n';
				cerr << "integral: " << temp << '\n';
				cerr << chronogram->GetTree()->GetLeftMost(from) << '\t' << chronogram->GetTree()->GetRightMost(from) << '\n';
			}
			exit(1);
		}
		double tmp = logit ? log(temp / (1.0 - temp)) : log(temp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;

		if (withpp)	{
			if (temp > threshold)	{
				pp[from->GetBranch()]++;
			}
		}

		dist[from->GetNode()].push_front(tmp);

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out(), branchmean);
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(BranchVarTree<UnitReal>* sample, LengthTree* chronogram, Link* from)	{
		double temp = (from->isRoot()) ? 1.0 : ((double) sample->GetBranchVal(from->GetBranch())->val());
		double tmp = logit ? log(temp / (1.0 - temp)) : log(temp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;

		dist[from->GetNode()].push_front(tmp);

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveReset(Link* from)	{
		mean[from->GetNode()] = 0;
		var[from->GetNode()] = 0;
		meanlog[from->GetNode()] = 0;
		varlog[from->GetNode()] = 0;
		dist[from->GetNode()].clear();
		if (withpp)	{
			pp[from->GetBranch()] = 0;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			meantime[link->GetBranch()] = 0;
			vartime[link->GetBranch()] = 0;
		}
	}

	void RecursiveNormalise(Link* from)	{
		mean[from->GetNode()] /= size;
		var[from->GetNode()] /= size;
		var[from->GetNode()] -= mean[from->GetNode()] * mean[from->GetNode()];
		meanlog[from->GetNode()] /= size;
		varlog[from->GetNode()] /= size;
		varlog[from->GetNode()] -= meanlog[from->GetNode()] * meanlog[from->GetNode()];

		dist[from->GetNode()].sort();

		if (withpp)	{
			pp[from->GetBranch()] /= size;
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] -= meantime[link->GetBranch()] * meantime[link->GetBranch()];
		}
	}

	private:

	const Tree* tree;
	int size;
	bool logit;

    double logscale;

	int leafsize;
	double meanleaf;
	double meanroot;
	double ppleafroot;

	bool printlog;
	bool printmean;
    bool printmed;
	bool printci;
	bool printstdev;

	bool withleaf;
	bool withinternal;

	bool withdepth;
	LengthTree* lengthtree;

	double meanreg;
	double stdevreg;

	map<const Branch*,double> pp;
	bool withpp;
	double threshold;

	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;
	map<const Node*,double> mean;
	map<const Node*,double> var;
	map<const Node*,double> meanlog;
	map<const Node*,double> varlog;

	map<const Node*, list<double> > dist;

};

#endif

