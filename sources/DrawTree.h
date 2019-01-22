
#ifndef DRAWTREE_H
#define DRAWTREE_H

#include "Tree.h"
#include <cmath>

class FileTree : public Tree	{

	public:

	FileTree(string filename) : Tree(filename) , leafval(-1), internalval(-1), rescale(0) {
		Check(GetRoot());
	}

	double GetNodeMin(const Link* from);
	double GetNodeMax(const Link* from);
	double GetNodeVal(const Link* from);

	string GetLeafNodeName(const Link* from);

	bool HasNodeMinMax() const {
		return ((internalval == 2) && (leafval == 2));
	}

	bool HasNodeVal() const {
		return ((internalval == 1) && (leafval == 1));
	}

	bool HasLeafMinMax() const {
		return (leafval == 2);
	}

	bool HasLeafVal() const {
		return (leafval == 1);
	}

	bool HasInternalMinMax() const {
		return (internalval == 2);
	}

	bool HasInternalVal() const {
		return (internalval == 1);
	}

	void Tabulate(ostream& os)	{
		RecursiveTabulate(GetRoot(),os);
	}

	void SetRescalingFactor(double in)	{
		rescale = in;
	}

	protected:

	void RecursiveTabulate(const Link* from, ostream& os)	{
		os << GetLeafNodeName(GetLeftMostLink(from)) << '\t' << GetLeafNodeName(GetRightMostLink(from)) << '\t' << GetNodeVal(from) << '\n';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(link->Out(),os);
		}
	}

	void Check(const Link* from);

	int leafval;
	int internalval;
	double rescale;
};


class DrawTree	{

	// abstract class
	public:

	DrawTree() : prec(2), beamer(true), sizeX(6), sizeY(10), Z(0), fontsize(4), groupfontsize(6), header(true), shiftname(-0.02), texscale(1.0), withleafnames(true),  withnodetext(false) {}
	virtual ~DrawTree() {}

	void Draw(string targetfile);
	void Draw(string targetfile, double insizeX, double insizeY)	{
		SetScale(insizeX, insizeY);
		Draw(targetfile);
	}

	void SetDocumentStyle(string style)	{
		if (style == "beamer")	{
			beamer = true;
		}
		else if (style == "article")	{
			beamer = false;
		}
		else	{
			cerr << "in drawtree: does not recognize style\n";
			exit(1);
		}
	}

	void SetNodeText(bool innt)	{
		withnodetext = innt;
	}

	// virtual double GetNodeVal(const Link* from) = 0;

	void SetScale(double insizeX, double insizeY)	{
		sizeX = insizeX;
		sizeY = insizeY;
	}

	void SetGroupOffset(double inZ)	{
		Z = inZ;
	}

	double GetGroupOffset()	{
		return Z;
	}

	void SetTexScale(double inscale)	{
		texscale = inscale;
	}

	void SetWithLeafNames(bool in)	{
		withleafnames = in;
	}

	void SetFontSize(double infontsize)	{
		fontsize = infontsize;
	}

	void SetGroupFontSize(double infontsize)	{
		groupfontsize = infontsize;
	}

	void SetShiftName(double s)	{
		shiftname = s;
	}

	void SetWithHeader(bool h)	{
		header = h;
	}

	virtual double GetLength(const Link* link) const = 0;
	virtual string GetNodeName(const Link* link) const = 0;

	virtual string GetLeafNodeName(const Link* link) {
		return GetNodeName(link);
	}

	void WriteNodeText(const Link* from, ostream& os, double min, double max);

	virtual string GetPreLeafNodeName(const Link* link) {return "";}

	double GetDepth(const Link* link) {return depth[link];}
	virtual double GetDepth() {return GetDepth(GetRoot());}

	int GetSize(const Link* link) {return size[link];}

	virtual const Link* GetRoot() const = 0;

	virtual const Link* GetLCA(string tax1, string tax2) = 0;

	void SetGroups(string filename)	{

		ifstream is(filename.c_str());
		int N;
		is >> N;
		for (int i=0; i<N; i++)	{
			string tax1, tax2, name;
			int text,color;
			is >> name >> tax1 >> tax2 >> text >> color;
            cerr << name << '\t' << tax1 << '\t' << tax2 << '\t' << text << '\t' << color << '\n';
			const Link* link = GetLCA(tax1,tax2);
			if (link)	{
				groupname[link->GetNode()] = name;
				grouptext[link->GetNode()] = text;
				groupcolor[link->GetNode()] = color;
			}
			else	{
				cerr << "did not find common ancestor of: " << tax1 << '\t' << tax2 << '\n';
				exit(1);
			}
		}
	}	

	protected:

	// node size: number of spanning leaves
	void ComputeSize()	{
		size.clear();
		RecursiveComputeSize(GetRoot());
	}
	int RecursiveComputeSize(const Link* from);

	// node depth: max length from this node down to its leaves
	void ComputeDepth()	{
		depth.clear();
		RecursiveComputeDepth(GetRoot());
	}
	double RecursiveComputeDepth(const Link* from);

	double	RecursiveDraw(const Link* from, ostream& os, double X, double Y, double xscale, double yscale);
	double	RecursiveAfterDraw(const Link* from, ostream& os, double X, double Y, double xscale, double yscale);

	virtual void LocalAfterDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{}

	virtual void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawNodeText(from,os,x,y);
	}

	virtual void LocalDrawVerticalTrait(const Link* from, ostream& os, double x, double y1, double y2);
	virtual void LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx);

	void DrawNodeText(const Link* from, ostream& os, double x, double y);

	void WriteGroupNames(ostream& os);

	virtual void PrepareDrawing();
	virtual void FinishDrawing(ostream& os, double xscale, double yscale) {}

	inline double texapprox(double f)	{
		double tmp = ((double) ((int) (f * 10000))) / 10000;
		if (fabs(tmp) < 1e-6)	{
			tmp = 0;
		}
		return tmp;
	}

	inline double approx(double f)	{
		for (int i=0; i<prec; i++)	{
			f *= 10;
		}
		f = ((double) ((int) f));
		for (int i=0; i<prec; i++)	{
			f /= 10;
		}
		return f;
	}

	int prec;

	bool beamer;
	double sizeX;
	double sizeY;
	double Z;
	double  fontsize;
	double  groupfontsize;

	map<const Link*, int> size;
	map<const Link*, double> depth;

	bool header;

	double shiftname;

	double texscale;

	bool withleafnames;
	bool withnodetext;

	map<const Node*, string> groupname;
	map<const Node*, double> groupy;
	map<const Node*, int> grouptext;
	map<const Node*, int> groupcolor;


};

class FileDrawTree : public FileTree, public DrawTree	{

	public:

	FileDrawTree(string filename) : FileTree(filename) {}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	Link* GetRoot() const {return Tree::GetRoot();}
	const Link* GetLCA(string tax1, string tax2) {return Tree::GetLCA(tax1,tax2);}

	virtual string GetLeafNodeName(const Link* link)  {
		return FileTree::GetLeafNodeName(link);
	}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		return atof(from->GetBranch()->GetName().c_str());
	}

};


class BubbleTree : public virtual DrawTree	{

	public:

	BubbleTree() : maxnodeval(10.0), nodepower(1.0) , withbubbletext(false) {}

	virtual double GetNodeMin(const Link* from) = 0;
	virtual double GetNodeMax(const Link* from) = 0;
	virtual double GetNodeVal(const Link* from) = 0;

	virtual double GetCircleNodeVal(const Link* link);
	virtual double GetCircleNodeMin(const Link* link);
	virtual double GetCircleNodeMax(const Link* link);

	void SetNodeScale(double scale)	{
		maxnodeval = scale;
	}

	void SetNodePower(double p)	{
		nodepower = p;
	}

	void SetBubbleText(bool innt)	{
		withbubbletext = innt;
	}

	protected:

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
		DrawNodeBubble(from,os,x,y);
		if ((!from->Out()->isLeaf()) && withbubbletext)	{
			WriteBubbleText(from,os,x,y);
		}
	}

	void DrawNodeBubble(const Link* from, ostream& os, double x, double y);
	void WriteBubbleText(const Link* from, ostream& os, double min, double max);

	void PrepareDrawing();
	void FinishDrawing(ostream& os, double xscale, double yscale);

	void DrawBubbleScale(ostream& os, double x, double y, double min, double max, double step);

	double ComputeMaxNodeVal(const Link* from);
	double ComputeMinNodeVal(const Link* from);

	double  maxnodeval;
	double nodescale;
	double nodepower;

	bool withbubbletext;
};

class FileBubbleTree : public FileTree, public BubbleTree	{

	public:

	FileBubbleTree(string filename) : FileTree(filename) {
		if (! HasNodeMinMax())	{
			cerr << "error : filebubbletree should have min and max values associated to each node\n";
		}
	}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	virtual string GetLeafNodeName(const Link* link)  {
		return FileTree::GetLeafNodeName(link);
	}

	Link* GetRoot() const {return Tree::GetRoot();}
	const Link* GetLCA(string tax1, string tax2) {return Tree::GetLCA(tax1,tax2);}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		return atof(from->GetBranch()->GetName().c_str());
	}

	/*
	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
		BubbleTree::LocalDraw(from,os,x,y,xscale,yscale);
	}

	void PrepareDrawing()	{
		DrawTree::PrepareDrawing();
		BubbleTree::PrepareDrawing();
	}
	*/

	double GetNodeMin(const Link* from)	{
		return FileTree::GetNodeMin(from);
	}
	double GetNodeMax(const Link* from)	{
		return FileTree::GetNodeMax(from);
	}
	double GetNodeVal(const Link* from)	{
		return FileTree::GetNodeVal(from);
	}
};

class ChronoDrawTree : public virtual DrawTree	{

	public:

	ChronoDrawTree() : barwidth(0.04), withtimescale(true), onlygroups(false), maxtime(0) {}

	virtual double GetNodeMin(const Link* from) = 0;
	virtual double GetNodeMax(const Link* from) = 0;
	virtual double GetMinTime(const Link* from) {return GetNodeMin(from);}
	virtual double GetMaxTime(const Link* from) {return GetNodeMax(from);}

	double GetMaxTime()	{
		if (maxtime)	{
			return maxtime;
		}
		return GetMaxTime(GetRoot());
	}

	void SetMaxTime(double intime)	{
		maxtime = intime;
	}

	void SetBarWidth(double inw)	{
		barwidth = inw;
	}

	void SetTimeScale(bool in)	{
		withtimescale = in;
	}

	void SetOnlyGroups(bool in)	{
		onlygroups = in;
	}

	double GetDepth()	{
		if (maxtime)	{
			return maxtime;
		}
		return DrawTree::GetDepth(GetRoot());
	}

	protected:

	void FinishDrawing(ostream& os, double xscale, double yscale) {
		DrawTimeScale(os,xscale,yscale);
	}

	void DrawTimeScale(ostream& os, double xscale, double yscale);

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTimeInterval(from,os,x,y,xscale,yscale);
	}

	void DrawTimeInterval(const Link* from, ostream& os, double x, double y, double xscale, double yscale);

	// double GetScale();
	double barwidth;
	bool withtimescale;
	bool onlygroups;
	double maxtime;

};

class FileChronoDrawTree : public FileTree, public ChronoDrawTree	{

	public:

	FileChronoDrawTree(string filename) : FileTree(filename) {
		if (! HasInternalMinMax())	{
			cerr << "error : filechronodrawtree should have min and max values associated to each internal node\n";
		}
	}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	virtual string GetLeafNodeName(const Link* link)  {
		return FileTree::GetLeafNodeName(link);
	}

	Link* GetRoot() const {return Tree::GetRoot();}
	const Link* GetLCA(string tax1, string tax2) {return Tree::GetLCA(tax1,tax2);}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		double x = atof(from->GetBranch()->GetName().c_str());
		return x;
	}

	double GetNodeMin(const Link* from)	{
		return FileTree::GetNodeMin(from);
	}
	double GetNodeMax(const Link* from)	{
		return FileTree::GetNodeMax(from);
	}
	/*
	double GetMinTime(const Link* from);
	double GetMaxTime(const Link* from);
	*/

	protected:

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
		ChronoDrawTree::LocalDraw(from,os,x,y,xscale,yscale);
	}

	void PrepareDrawing()	{
		DrawTree::PrepareDrawing();
		ChronoDrawTree::PrepareDrawing();
	}

};

#endif


class HeatTree : public virtual DrawTree	{

	public:

	HeatTree() : maxnodeval(-1), minnodeval(-1), nodepower(1.0) , thickness(0.06), unibranch(false), withbranchval(false), withexternalnodeval(false), withexternalnodeci(false) {}

	void SetMinMax(double min, double max)	{
		minnodeval = min;
		maxnodeval = max;
	}

	void SetBranchVal(string infile);
	void SetExternalNodeVal(string infile);
	void SetExternalNodeCI(string infile);

	virtual const Link* GetLCA(string tax1, string tax2) = 0;

	virtual double GetNodeVal(const Link* from) = 0;

	void SetNodePower(double p)	{
		nodepower = p;
	}

	void SetThickness(double inth)	{
		thickness = inth;
	}

	void SetUniBranch(bool inb)	{
		unibranch = inb;
	}

	/*
	void SetNodeText(bool innt)	{
		withnodetext = innt;
	}
	*/

	protected:

	void LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx);
	// void WriteNodeText(const Link* from, ostream& os, double min, double max);

	void PrepareDrawing();
	void FinishDrawing(ostream& os, double xscale, double yscale);
	void MakeScale(ostream& os);

	void DrawCI(const Link* from, ostream& os, double x, double y);

	void LocalAfterDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		if (withexternalnodeci)	{
			if (! from->isLeaf())	{
				DrawCI(from,os,x,y);
			}
		}
	}

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
	}

	string GetColorCode(double z);

	// void DrawHeatScale(ostream& os, double x, double y, double min, double max, double length);

	double ComputeMaxNodeVal(const Link* from);
	double ComputeMinNodeVal(const Link* from);

	double  maxnodeval;
	double  minnodeval;
	double nodepower;
	double thickness;
	// bool withnodetext;
	bool unibranch;
	bool withbranchval;
	bool withexternalnodeval;
	bool withexternalnodeci;
	map<const Branch*, string> branchval;
	// map<const Branch*, double> branchval;
	map<const Node*, double> extnodeval;
	map<const Node*, double> extnodemin;
	map<const Node*, double> extnodemax;

};

class FileHeatTree : public FileTree, public HeatTree	{

	public:

	FileHeatTree(string filename) : FileTree(filename) {
		if (! HasNodeVal())	{
			cerr << "error in fileheattree: only one value should be associated with each node\n";
			cerr << "leaf: " << HasLeafVal() << '\n';
			cerr << "internal : " << HasInternalVal() << '\n';
			exit(1);
		}
	}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}

	virtual string GetPreLeafNodeName(const Link* link) {

		ostringstream s;

		if (withexternalnodeci)	{
			map<const Node*, double>::const_iterator imin = extnodemin.find(link->GetNode());
			map<const Node*, double>::const_iterator imax = extnodemax.find(link->GetNode());

			s.precision(2);

			Link* link2 = (Link*) link;
			double z = GetNodeVal(link2);
			s << fixed << z << " ";

			s << "(";
			if (imin == extnodemin.end())	{
				s << "?";
			}
			else	{
				double tmp = imin->second;
				s << fixed << tmp;
			}
			s << " , ";
			if (imin == extnodemax.end())	{
				s << "?";
			}
			else	{
				double tmp = imax->second;
				s << fixed << tmp;
			}
			s << ")";
		}

		// s << FileTree::GetLeafNodeName(link);
		return s.str();
	}

	virtual string GetLeafNodeName(const Link* link) {

		ostringstream s;

		/*
		if (withexternalnodeci)	{
			map<const Node*, double>::const_iterator imin = extnodemin.find(link->GetNode());
			map<const Node*, double>::const_iterator imax = extnodemax.find(link->GetNode());

			s.precision(2);
			// s.setw(2);
			// s.setf ( ios::floatfield|ios::showpoint|ios::fixed );
			// os << setw(10) << GetTotalLength();

			s << "(";
			if (imin == extnodemin.end())	{
				s << "?";
			}
			else	{
				double tmp = imin->second;
				s << fixed << tmp;
				// s << ((double) ((int) ((100) * tmp))) / 100;
			}
			s << " , ";
			if (imin == extnodemax.end())	{
				s << "?";
			}
			else	{
				double tmp = imax->second;
				s << fixed << tmp;
				// s << ((double) ((int) ((100) * tmp))) / 100;
			}
			s << ") \\,  ";
		}
		*/
		s << FileTree::GetLeafNodeName(link);
		if (withnodetext)	{
			s << ' ' << ((double) ((int) (10 * GetNodeVal(link))) ) / 10;
		}

		return s.str();
	}

	Link* GetRoot() const {return Tree::GetRoot();}
	const Link* GetLCA(string tax1, string tax2) {return Tree::GetLCA(tax1,tax2);}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		double x = atof(from->GetBranch()->GetName().c_str());
		return x;
	}

	double GetExternalNodeMin(const Link* from)	{
		return extnodemin[from->GetNode()];
	}

	double GetExternalNodeMax(const Link* from)	{
		return extnodemax[from->GetNode()];
	}

	double GetNodeVal(const Link* from) {
		if (withexternalnodeval)	{
			return extnodeval[from->GetNode()];
		}
		return FileTree::GetNodeVal(from);
	}
};
