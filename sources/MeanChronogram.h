
#ifndef MEANCHRONO_H
#define MEANCHRONO_H

#include "CalibratedChronogram.h"
#include <map>
#include <list>

class MeanChronogram : public NewickTree {

	public:

	MeanChronogram(Tree* intree, bool inprintci = true, bool inprintmean=false, bool inprintstdev = false, bool inprintmedian = false) : tree(intree), printci(inprintci), printmean(inprintmean), printstdev(inprintstdev), printmedian(inprintmedian), withleafdates(false)  {
		Reset();
	}

	void SetWithLeafDates(bool in)	{
		withleafdates = in;
	}

	Tree* GetTree() const {return tree;}

	const Link* GetRoot() const {return GetTree()->GetRoot();}

	double GetMeanDate(const Node* node) const	{
		map<const Node*, double>::const_iterator i = meandate.find(node);
		return i->second;
	}

	double GetVarDate(const Node* node) const	{
		map<const Node*, double>::const_iterator i = vardate.find(node);
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

	double GetMedian(const Node* node) const {
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 2));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return *i;
	}

	double GetMin95(const Node* node) const {
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 2.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return *i;
	}

	double GetMax95(const Node* node) const {
		map<const Node*, list<double> >::const_iterator f = dist.find(node);
		list<double> l = f->second;
		list<double>::const_iterator i = l.begin();
		int n = ((int) (((double) l.size()) / 100 * 97.5));
		for (int j=0; j<n; j++)	{
			i++;
		}
		return *i;
	}

	string GetNodeName(const Link* link) const {
		ostringstream s;
		if (link->isLeaf())	{
			s << link->GetNode()->GetName();
		}
		if ((! link->isLeaf()) || withleafdates)	{
			if (printci)	{
				s << GetMin95(link->GetNode()) << "_" << GetMax95(link->GetNode());
			}
			else	{
				if (printmean)	{
					s << GetMeanDate(link->GetNode());
				}
				else	{
					s << GetMedian(link->GetNode());
				}
				if (printstdev)	{
					s << "_";
					if (fabs(GetVarDate(link->GetNode())) > 1e-6)	{
						s << sqrt(GetVarDate(link->GetNode())) ;
					}
					else	{
						s << 0;
					}
				}
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
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(Chronogram* chronogram)	{
		double scale = 1;
		CalibratedChronogram* tmp = dynamic_cast<CalibratedChronogram*>(chronogram);
		if (tmp)	{
			scale = tmp->GetScale()->val();
		}
		RecursiveAdd(chronogram, scale, GetTree()->GetRoot());
		size++;
	}

	/*
	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->ToStream(os);
	}
	*/

    void TabulateHeader(ostream& os)    {
        os << "#taxon1\ttaxon2";
        if (printmean)	{
            os << "\tmean";
        }
        if (printmedian)	{
            os << "\tmedian";
        }
        if (printci)	{
            os << "\tmin95\tmax95";
        }
        if (printstdev)	{
            os << "\tstdev";
        }
        os << '\n';
    }

	void Tabulate(ostream& os)	{
        TabulateHeader(os);
		RecursiveTabulate(os,GetTree()->GetRoot());
	}

	void TabulateTimes(ostream& os)	{
		RecursiveTabulateTimes(os,GetTree()->GetRoot());
	}

	private:

	void RecursiveTabulate(ostream& os, const Link* from)	{
		if (! from->isLeaf())	{
		// if (fabs(vardate[from->GetNode()]) > 1e-6)	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from);
			if (printmean)	{
				os << '\t' << meandate[from->GetNode()];
			}
			if (printmedian)	{
				os << '\t' << GetMedian(from->GetNode());
			}
			if (printci)	{
				os << '\t' << GetMin95(from->GetNode()) << '\t' << GetMax95(from->GetNode());
			}
			if (printstdev)	{
				os << '\t' << sqrt(vardate[from->GetNode()]);
			}
			os << '\n';
		}
		/*
		else	{
			os << meandate[from->GetNode()] << '\t' << 0 << '\n';
		}
		*/
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(os,link->Out());
		}
	}

	void RecursiveTabulateTimes(ostream& os, const Link* from)	{
		os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
		os << meantime[from->GetBranch()] << '\t' << sqrt(vartime[from->GetBranch()]) << '\n';
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulateTimes(os,link->Out());
		}
	}


	void RecursiveAdd(Chronogram* chronogram, double scale, const Link* from)	{
		double date = chronogram->GetNodeVal(from->GetNode())->val() * scale;
		meandate[from->GetNode()] += date;
		vardate[from->GetNode()] += date*date;
		dist[from->GetNode()].push_front(date);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(chronogram, scale, link->Out());
			double time = chronogram->GetBranchTimeLength(link->Out()) * scale;
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveReset(const Link* from)	{
		meandate[from->GetNode()] = 0;
		vardate[from->GetNode()] = 0;
		dist[from->GetNode()].clear();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			meantime[link->GetBranch()] = 0;
			vartime[link->GetBranch()] = 0;
		}
	}

	void RecursiveNormalise(const Link* from)	{
		meandate[from->GetNode()]/=size;
		vardate[from->GetNode()]/=size;
		vardate[from->GetNode()] -= meandate[from->GetNode()] * meandate[from->GetNode()];
		dist[from->GetNode()].sort();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] -= meantime[link->GetBranch()] * meantime[link->GetBranch()];
		}
	}

	/*
	void RecursiveSetName(ostream& os, const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetName(os, link->Out());
			ostringstream s;
			s << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s.str());
		}
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
		}
		else	{
			if (fabs(vardate[from->GetNode()]) > 1e-6)	{
				s << meandate[from->GetNode()] << "_" << sqrt(vardate[from->GetNode()]) ;
			}
			else	{
				s << meandate[from->GetNode()] << "_" << 0;
			}
		}
		from->GetNode()->SetName(s.str());
	}
	*/

	Tree* tree;
	int size;

	bool printci;
	bool printmean;
	bool printmedian;
	bool printstdev;

	bool withleafdates;

	map<const Node*,double> meandate;
	map<const Node*,double> vardate;
	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;

	map<const Node*, list<double> > dist;
};

#endif

