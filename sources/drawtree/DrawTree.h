
#ifndef DRAWTREE_H
#define DRAWTREE_H

#include "Tree.h"
#include <cmath>

class FileTree : public Tree	{

	public:

	FileTree(string filename) : Tree(filename) , leafval(-1), internalval(-1) {
		Check(GetRoot());
	}

	double GetNodeMin(const Link* from);
	double GetNodeMax(const Link* from);
	double GetNodeVal(const Link* from);

	string GetLeafNodeName(const Link* from) const;

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

	protected:

	void Check(const Link* from);

	int leafval;
	int internalval;

	
};


class DrawTree	{

	// abstract class
	public:

	DrawTree() : prec(2), beamer(true), sizeX(6), sizeY(10), fontsize(4), header(true), shiftname(0.04) {}
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

	void SetScale(double insizeX, double insizeY)	{
		sizeX = insizeX;
		sizeY = insizeY;
	}

	void SetFontSize(double infontsize)	{
		fontsize = infontsize;
	}

	void SetShiftName(double s)	{
		shiftname = s;
	}

	void SetWithHeader(bool h)	{
		header = h;
	}

	virtual double GetLength(const Link* link) const = 0;
	virtual string GetNodeName(const Link* link) const = 0;
	virtual string GetLeafNodeName(const Link* link)  const {
		return GetNodeName(link);
	}

	double GetDepth(const Link* link) {return depth[link];}
	int GetSize(const Link* link) {return size[link];}

	virtual const Link* GetRoot() const = 0;

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

	virtual void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawNodeText(from,os,x,y);
	}

	virtual void LocalDrawVerticalTrait(const Link* from, ostream& os, double x, double y1, double y2);
	virtual void LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx);

	void DrawNodeText(const Link* from, ostream& os, double x, double y);

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
	double  fontsize;

	map<const Link*, int> size;
	map<const Link*, double> depth;

	bool header;

	double shiftname;


};

class FileDrawTree : public FileTree, public DrawTree	{

	public:

	FileDrawTree(string filename) : FileTree(filename) {}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	Link* GetRoot() const {return Tree::GetRoot();}

	virtual string GetLeafNodeName(const Link* link)  const {
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
	virtual string GetLeafNodeName(const Link* link)  const {
		return FileTree::GetLeafNodeName(link);
	}

	Link* GetRoot() const {return Tree::GetRoot();}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		return atof(from->GetBranch()->GetName().c_str());
	}

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
		BubbleTree::LocalDraw(from,os,x,y,xscale,yscale);
	}

	void PrepareDrawing()	{
		DrawTree::PrepareDrawing();
		BubbleTree::PrepareDrawing();
	}

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

	ChronoDrawTree() : barwidth(0.04), withtimescale(true) {}

	virtual double GetNodeMin(const Link* from) = 0;
	virtual double GetNodeMax(const Link* from) = 0;
	virtual double GetMinTime(const Link* from) {return GetNodeMin(from);}
	virtual double GetMaxTime(const Link* from) {return GetNodeMax(from);}

	void SetBarWidth(double inw)	{
		barwidth = inw;
	}

	void SetTimeScale(bool in)	{
		withtimescale = in;
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

};

class FileChronoDrawTree : public FileTree, public ChronoDrawTree	{

	public:

	FileChronoDrawTree(string filename) : FileTree(filename) {
		if (! HasInternalMinMax())	{
			cerr << "error : filechronodrawtree should have min and max values associated to each internal node\n";
		}
	}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	virtual string GetLeafNodeName(const Link* link)  const {
		return FileTree::GetLeafNodeName(link);
	}

	Link* GetRoot() const {return Tree::GetRoot();}

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

	HeatTree() : maxnodeval(10.0), minnodeval(0), nodepower(1.0) , thickness(0.06), withnodetext(false), unibranch(false) {}

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

	void SetNodeText(bool innt)	{
		withnodetext = innt;
	}

	protected:

	void LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx);
	void WriteNodeText(const Link* from, ostream& os, double min, double max);

	void PrepareDrawing();
	// void FinishDrawing(ostream& os, double xscale, double yscale);

	string GetColorCode(double z);

	// void DrawHeatScale(ostream& os, double x, double y, double min, double max, double length);

	double ComputeMaxNodeVal(const Link* from);
	double ComputeMinNodeVal(const Link* from);

	double  maxnodeval;
	double  minnodeval;
	double nodepower;
	double thickness;
	bool withnodetext;
	bool unibranch;
};

class FileHeatTree : public FileTree, public HeatTree	{

	public:

	FileHeatTree(string filename) : FileTree(filename) {
		if (! HasNodeVal())	{
			cerr << "error in fileheattree: only one value should be associated with each node\n";
			exit(1);
		}
	}

	string GetNodeName(const Link* link) const {return Tree::GetNodeName(link);}
	virtual string GetLeafNodeName(const Link* link)  const {
		return FileTree::GetLeafNodeName(link);
	}

	Link* GetRoot() const {return Tree::GetRoot();}

	double GetLength(const Link* from)	const {
		if (! from->GetBranch())	{
			return 0;
		}
		double x = atof(from->GetBranch()->GetName().c_str());
		return x;
	}

	double GetNodeVal(const Link* from)	{
		return FileTree::GetNodeVal(from);
	}
};
