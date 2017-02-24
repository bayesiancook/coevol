
#ifndef MEANVALTREE
#define MEANVALTREE

#include "ValTree.h"
#include "Var.h"
#include "MultiVariateTreeProcess.h"
#include <cmath>
#include <list>

class MeanBranchTree : public NewickTree {

	public:

	MeanBranchTree(Tree* intree, bool inwithroot = false) : tree(intree), withRoot(inwithroot)	{

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

	/*
	MeanExpNormTree(Tree* intree, bool inlogit) : tree(intree), logit(inlogit), printlog(false), printmean(false), printci(true), printstdev(false) {
		Reset();
	}
	*/

	MeanExpNormTree(Tree* intree, bool inlogit, bool inprintlog, bool inprintmean, bool inprintci, bool inprintstdev, bool inwithleaf, bool inwithinternal, double inmeanreg = 0, double instdevreg = 0) : tree(intree), logit(inlogit), printlog(inprintlog), printmean(inprintmean), printci(inprintci), printstdev(inprintstdev), withleaf(inwithleaf), withinternal(inwithinternal), meanreg(inmeanreg), stdevreg(instdevreg) {
		Reset();
		ppleafroot = 0;
		threshold = 0;
		withpp = false;
	}

	void ActivatePP(double inthreshold)	{
		threshold = inthreshold;
		withpp = true;
		Reset();
	}

	void SetLog() {printlog = true;}

	void SetPrintMean(bool inprintmean = true) {printmean = inprintmean;}

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
				return tmp;
			}
			else if (logit)	{
				return pow(10,tmp) / (1 + pow(10,tmp));
			}
			return pow(10,tmp);
		}
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 2.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return printlog ? *i : (logit? pow(10,*i) / (1 + pow(10,*i)) : pow(10,*i));
	}

	double GetMax95(const Node* node) const {
		if (meanreg)	{
			double tmp = _GetMeanLog(node) + 1.959964 * sqrt(_GetVarLog(node));
			if (printlog)	{
				return tmp;
			}
			else if (logit)	{
				return pow(10,tmp) / (1 + pow(10,tmp));
			}
			return pow(10,tmp);
		}
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 97.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return printlog ? *i : (logit ? pow(10,*i) / (1 + pow(10,*i)) : pow(10,*i));
	}

	double _GetMean(const Node* node) const	{
		map<const Node*, double>::const_iterator i = mean.find(node);
		return i->second;
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
		return i->second;
	}

	double _GetVarLog(const Node* node) const	{
		map<const Node*, double>::const_iterator i = varlog.find(node);
		map<const Node*, double>::const_iterator j = meanlog.find(node);
		if (meanreg)	{
			return meanreg * meanreg * i->second;
			// return meanreg * meanreg * i->second + j->second * j->second * stdevreg * stdevreg;
		}
		return i->second;
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
		if ((! printmean) || (! isFixed(link->GetNode())))	{
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
	
	
	void AddNe(NodeVarTree<RealVector>* sample, LengthTree* chronogram, double alpha[], double beta, int dim, int indice1, int indice2, Var<Real>* Offset = 0)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		double offset = 0;
		if (Offset)	{
			offset = Offset->val();
		}
		RecursiveAddNe(sample, chronogram, GetTree()->GetRoot(), alpha, beta, dim, indice1, indice2, offset);
		meanleaf /= leafsize;
		if (meanleaf > meanroot)	{
			ppleafroot++;
		}
		size++;
	}

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

	void Add(NodeVarTree<RealVector>* sample, LengthTree* chronogram, int index, double offset)	{
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

	void Add(NodeVarTree<Real>* sample, LengthTree* chronogram)	{
		meanleaf = 0;
		meanroot = 0;
		leafsize = 0;
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
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

	void Tabulate(ostream& os)	{
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
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			if (withdepth)	{
				os << GetDepth(from) << '\t';
			}
			else	{
				os << GetMeanTime(from->GetBranch()) << '\t';
			}
			if (printmean)	{
				os << GetMean(from->GetNode()) << '\t';
			}
			if (printstdev)	{
				if (! isFixed(from->GetNode()))	{
					os << sqrt(GetVar(from->GetNode())) << '\t';
				}
				else	{
					os << 0 << '\t';
				}
			}
			if (printci)	{
				os << GetMin95(from->GetNode()) << '\t' << GetMax95(from->GetNode()) << '\t';
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
	
	
	void RecursiveAddNe(NodeVarTree<RealVector>* sample, LengthTree* chronogram, Link* from, double alpha[], double beta, int dim, int indice1, int indice2, double offset)	{
		double tmp(0);
		tmp += alpha[0] * ((* sample->GetNodeVal(from->GetNode()))[0])/log(10);
		tmp += alpha[indice1] * ((* sample->GetNodeVal(from->GetNode()))[indice1])/log(10);
		tmp += alpha[indice2] * ((* sample->GetNodeVal(from->GetNode()))[indice2])/log(10);
		tmp += beta;
		tmp += offset;
		if (from->isRoot())	{
			meanroot = tmp;
		}
		else if (from->isLeaf())	{
			meanleaf += tmp;
			leafsize++;
		}
		double temp = logit ? pow(10,tmp) / (1 + pow(10,tmp)) : pow(10,tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
		dist[from->GetNode()].push_front(tmp);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddNe(sample, chronogram, link->Out(), alpha, beta, dim, indice1, indice2, offset);
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

	void RecursiveAdd(NodeVarTree<Real>* sample, LengthTree* chronogram, Link* from)	{
		double tmp = sample->GetNodeVal(from->GetNode())->val();
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
			RecursiveAdd(sample, chronogram, link->Out());
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
		if (isnan(temp))	{
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

	int leafsize;
	double meanleaf;
	double meanroot;
	double ppleafroot;

	bool printlog;
	bool printmean;
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




// a version with PopSize functions

/*
class MeanExpNormTree : public NewickTree {

	public:

	MeanExpNormTree(Tree* intree, bool inlogit = false) : tree(intree), logit(inlogit)  {
		Reset();
	}

	Tree* GetTree() const {return tree;}

	const Link* GetRoot() const {return GetTree()->GetRoot();}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(MultiVariateTreeProcess* sample, LengthTree* chronogram, int index)	{
		RecursiveAdd(sample, chronogram, index, GetTree()->GetRoot());
		size++;
	}

	void AddGC(MultiVariateTreeProcess* sample, LengthTree* chronogram, int index)	{
		RecursiveAddGC(sample, chronogram, index, GetTree()->GetRoot());
		size++;
	}

	void Add(NodeVarTree<Real>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->ToStream(os);
	}

	void PrintPopSizeTree2(ostream& os, double refpopsize1, double refpopsize2, string reftaxon1, string reftaxon2)	{
		double logomega1 = 0;
		double logomega2 = 0;
		GetPopSizeOffset2(GetTree()->GetRoot(),reftaxon1, reftaxon2, logomega1, logomega2);
		double alpha = (logomega2 - logomega1) / (log(refpopsize1) - log(refpopsize2));
		double beta = alpha * log(refpopsize1) + logomega1;
		RecursivePopSizeTree(os,GetTree()->GetRoot(), alpha, beta);
		tree->ToStream(os);
	}

	void PrintPopSizeTree(ostream& os, double alpha, double refpopsize, string reftaxon)	{
		double beta = GetPopSizeOffset(GetTree()->GetRoot(),alpha,refpopsize,reftaxon);
		RecursivePopSizeTree(os,GetTree()->GetRoot(), alpha, beta);
		tree->ToStream(os);
	}

	void Tabulate(ostream& os, bool leafonly)	{
		RecursiveTabulate(os,GetTree()->GetRoot(),leafonly);
	}

	private:

	void RecursiveTabulate(ostream& os, Link* from, bool leafonly)	{
		if ((! leafonly) || (from->isLeaf()))	{
			if (fabs(var[from->GetNode()]) > 1e-6)	{
				os << meanlog[from->GetNode()] / log(10.0)  << '\t' << sqrt(varlog[from->GetNode()]) / log(10.0)  << '\n';
			}
			else	{
				os << meanlog[from->GetNode()] / log(10.0) << '\t' << 0 << '\n';
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(os,link->Out(),leafonly);
		}
	}

	void RecursiveAddGC(MultiVariateTreeProcess* sample, LengthTree* chronogram, int index, Link* from)	{
		double tmp1 = exp((* sample->GetMultiNormal(from))[index+1]) + exp((* sample->GetMultiNormal(from))[index+2]);
		double tmp2 = exp((* sample->GetMultiNormal(from))[index+0]) + exp((* sample->GetMultiNormal(from))[index+3]);
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

	void RecursiveAdd(MultiVariateTreeProcess* sample, LengthTree* chronogram, int index, Link* from)	{
		double tmp = (* sample->GetMultiNormal(from))[index];
		double temp = logit ? exp(tmp) / (1 + exp(tmp)) : exp(tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, index, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(NodeVarTree<Real>* sample, LengthTree* chronogram, Link* from)	{
		double tmp = sample->GetNodeVal(from->GetNode())->val();
		double temp = logit ? exp(tmp) / (1 + exp(tmp)) : exp(tmp);
		meanlog[from->GetNode()] += tmp;
		varlog[from->GetNode()] += tmp * tmp;
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
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
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] -= meantime[link->GetBranch()] * meantime[link->GetBranch()];
		}
	}

	double GetPopSizeOffset(Link* from, double alpha, double refpopsize, string reftaxon)	{
		double beta = -1;
		if (from->isLeaf())	{
			string t = from->GetNode()->GetName();
			if (reftaxon == t)	{
				beta = alpha * log(refpopsize) + mean[from->GetNode()];
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = GetPopSizeOffset(link->Out(),alpha,refpopsize,reftaxon);
			if (tmp != -1)	{
				beta = tmp;
			}
		}
		return beta;
	}

	void GetPopSizeOffset2(Link* from, string reftaxon1, string reftaxon2, double& logomega1, double& logomega2)	{
		if (from->isLeaf())	{
			string t = from->GetNode()->GetName();
			ostringstream s;
			unsigned int k = 0;
			while ((k<t.length()) && (t[k] != '_'))	{
				s << t[k];
				k++;
			}

			if (reftaxon1 == s.str())	{
				logomega1 = mean[from->GetNode()];
			}
			if (reftaxon2 == s.str())	{
				logomega2 = mean[from->GetNode()];
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetPopSizeOffset2(link->Out(),reftaxon1, reftaxon2, logomega1, logomega2);
		}
	}

	void RecursivePopSizeTree(ostream& os, Link* from, double alpha, double beta)	{
		ostringstream s;
		if (from->isLeaf())	{
			string t = from->GetNode()->GetName();
			unsigned int i = 0;
			char c = ' ';
			while ((i<t.size()) && (c != '_'))	{
				c = t[i];
				if (c != '_')	{
					s << c;
				}
				i++;
			}
			s << '_';
		}
		s << exp((beta - mean[from->GetNode()]) / alpha);
		from->GetNode()->SetName(s.str());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivePopSizeTree(os, link->Out(),alpha,beta);
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		ostringstream s;
		if (from->isLeaf())	{
			string t = from->GetNode()->GetName();
			unsigned int i = 0;
			char c = ' ';
			while ((i<t.size()) && (c != '_'))	{
				c = t[i];
				if (c != '_')	{
					s << c;
				}
				i++;
			}
			s << '_';
		}
		s << mean[from->GetNode()];
		from->GetNode()->SetName(s.str());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetName(os, link->Out());
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	Tree* tree;
	double dalpha;
	int size;
	bool logit;

	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;
	map<const Node*,double> mean;
	map<const Node*,double> var;
	map<const Node*,double> meanlog;
	map<const Node*,double> varlog;

};



class MeanPopSizeTree	{

	public:

	MeanPopSizeTree(Tree* intree, double indalpha) : tree(intree) , dalpha(indalpha) {
		Reset();
	}

	Tree* GetTree() {return tree;}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(NodeVarTree<Real>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->ToStream(os);
	}

	private:

	void RecursiveAdd(NodeVarTree<Real>* sample, LengthTree* chronogram, Link* from)	{
		double tmp = sample->GetNodeVal(from->GetNode())->val();
		double temp = - tmp / dalpha / log(10.0);
		mean[from->GetNode()] += temp;
		var[from->GetNode()] += temp * temp;
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
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] -= meantime[link->GetBranch()] * meantime[link->GetBranch()];
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		ostringstream s;
		if (from->isLeaf())	{
			string t = from->GetNode()->GetName();
			unsigned int i = 0;
			char c = ' ';
			while ((i<t.size()) && (c != '_'))	{
				c = t[i];
				if (c != '_')	{
					s << c;
				}
				i++;
			}
			s << '_';
		}
		s << mean[from->GetNode()] / log(10.0);
		from->GetNode()->SetName(s.str());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetName(os, link->Out());
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	Tree* tree;
	double dalpha;
	int size;

	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;
	map<const Node*,double> mean;
	map<const Node*,double> var;

};


class MeanBranchTree	{

	public:

	MeanBranchTree(Tree* intree, bool inwithroot = false) : tree(intree), withRoot(inwithroot)	{

	}

	Tree* GetTree() {return tree;}

	Link* GetRoot() {return tree->GetRoot();}

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

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->ToStream(os);
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	private:

	void RecursiveAdd(BranchValPtrTree<Rvar<PosReal> >* sample, Link* from)	{
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

	void RecursiveAdd(BranchValPtrTree<Dvar<PosReal> >* sample, Link* from)	{
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

	void RecursiveReset(Link* from)	{
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

	void RecursiveNormalise(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			mean[link->GetBranch()] /= size;
			var[link->GetBranch()] /= size;
			var[link->GetBranch()] -= mean[link->GetBranch()] * mean[link->GetBranch()];
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetName(os, link->Out());
			ostringstream s;
			s <<  ((int) (100 * mean[link->GetBranch()])) + sqrt(var[link->GetBranch()]);
			link->GetBranch()->SetName(s.str());
		}
	}


	map<const Branch*,double> mean;
	map<const Branch*,double> var;
	Tree* tree;
	bool withRoot;
	int size;

};

class MeanRateBranchTree	{

	public:

	MeanRateBranchTree(Tree* intree, BranchValType inbval = INTEGRAL) : tree(intree), bval(inbval) {
		Reset();
	}

	Tree* GetTree() {return tree;}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(BranchVarTree<PosReal>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void Add(BranchVarTree<Real>* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->ToStream(os);
	}

	void Tabulate(ostream& os, bool leafonly, bool withtimes)	{
		RecursiveTabulate(os,GetTree()->GetRoot(),leafonly,withtimes);
	}

	private:

	void RecursiveTabulate(ostream& os, Link* from, bool leafonly, bool withtimes)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			if ((! leafonly) || (link->Out()->isLeaf()))	{
				if (withtimes)	{
					if (fabs(vartime[link->GetBranch()]) > 1e-6)	{
						os << meantime[link->GetBranch()] << '\t' << sqrt(vartime[link->GetBranch()])  << '\n';
					}
					else	{
						os << meantime[link->GetBranch()] << '\t' << 0 << '\n';
					}
				}
				if (fabs(var[link->GetBranch()]) > 1e-6)	{
					os << meanlog[link->GetBranch()] / log(10.0)  << '\t' << sqrt(varlog[link->GetBranch()]) / log(10.0)  << '\n';
				}
				else	{
					os << meanlog[link->GetBranch()] / log(10.0) << '\t' << 0 << '\n';
				}
			}
			RecursiveTabulate(os,link->Out(),leafonly,withtimes);
		}
	}


	void RecursiveAdd(BranchVarTree<PosReal>* sample, LengthTree* chronogram, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			if (bval == MEAN)	{
				tmp /= time;
			}
			double temp = log(tmp);
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
			meanlog[link->GetBranch()] += temp;
			varlog[link->GetBranch()] += temp * temp;
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveAdd(BranchVarTree<Real>* sample, LengthTree* chronogram, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out());
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			double tmp = sample->GetBranchVal(link->GetBranch())->val();
			if (bval == MEAN)	{
				tmp /= time;
			}
			double temp = log(tmp);
			mean[link->GetBranch()] += tmp;
			var[link->GetBranch()] += tmp * tmp;
			meanlog[link->GetBranch()] += temp;
			varlog[link->GetBranch()] += temp * temp;
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveReset(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			mean[link->GetBranch()] = 0;
			var[link->GetBranch()] = 0;
			meanlog[link->GetBranch()] = 0;
			varlog[link->GetBranch()] = 0;
			meantime[link->GetBranch()] = 0;
			vartime[link->GetBranch()] = 0;
		}
	}

	void RecursiveNormalise(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			mean[link->GetBranch()] /= size;
			var[link->GetBranch()] /= size;
			var[link->GetBranch()] -= mean[link->GetBranch()] * mean[link->GetBranch()];
			meanlog[link->GetBranch()] /= size;
			varlog[link->GetBranch()] /= size;
			varlog[link->GetBranch()] -= meanlog[link->GetBranch()] * meanlog[link->GetBranch()];
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] -= meantime[link->GetBranch()] * meantime[link->GetBranch()];
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetName(os, link->Out());
			ostringstream s;
			if (link->Out()->isLeaf())	{
				string t = link->Out()->GetNode()->GetName();
				unsigned int i = 0;
				char c = ' ';
				while ((i<t.size()) && (c != '_'))	{
					c = t[i];
					if (c != '_')	{
						s << c;
					}
					i++;
				}
				s << '_';
			}
			s << mean[link->GetBranch()];
			// s << mean[link->GetBranch()] << "_" << var[link->GetBranch()] ;
			link->Out()->GetNode()->SetName(s.str());
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	Tree* tree;
	BranchValType bval;
	int size;

	map<const Branch*,double> mean;
	map<const Branch*,double> var;
	map<const Branch*,double> meanlog;
	map<const Branch*,double> varlog;
	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;

};

*/

#endif

