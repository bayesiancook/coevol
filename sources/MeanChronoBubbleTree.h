
#ifndef MEANCHRONOBUBBLE_H
#define MEANCHRONOBUBBLE_H

#include "DrawTree.h"
#include "MeanChronogram.h"
#include "MeanValTree.h"

class MeanChronoBubbleTree : public virtual BubbleTree, public virtual ChronoDrawTree	{

	public:

	MeanChronoBubbleTree(MeanChronogram* inchrono, MeanExpNormTree* inbubble, double inxscale, double inyscale, double innodescale, double innodepower, double inbarwidth, int infontsize, bool inwithbubbletext, bool inwithheader=true, double inshiftname = 0.04)	{
		chrono = inchrono;
		bubble = inbubble;
		SetWithHeader(inwithheader);
		SetNodeScale(innodescale);
		SetNodePower(innodepower);
		SetBarWidth(inbarwidth);
		SetScale(inxscale,inyscale);
		SetFontSize(infontsize);
		SetBubbleText(inwithbubbletext);
		SetShiftName(inshiftname);
	}

	string GetNodeName(const Link* link) const {return chrono->GetNodeName(link);}
	const Link* GetRoot() const {return chrono->GetRoot();}
	const Link* GetLCA(string tax1, string tax2) {
		cerr << "in MeanChronoBubbleTree::GetLCA\n";
		exit(1);
		return 0;
		// return BubbleTree::GetLCA(tax1,tax2);}
	}

	double GetLength(const Link* from)	const {
		return chrono->GetMeanTime(from->GetBranch());
	}

	double GetNodeMin(const Link* from)	{
		return bubble->GetMin95(from->GetNode());
	}

	double GetNodeMax(const Link* from)	{
		return bubble->GetMax95(from->GetNode());
	}

	double GetNodeVal(const Link* from)	{
		return bubble->GetMean(from->GetNode());
	}

	double GetMinTime(const Link* from)	{
		return chrono->GetMin95(from->GetNode());
	}

	double GetMaxTime(const Link* from)	{
		return chrono->GetMax95(from->GetNode());
	}

	protected:

	void LocalDraw(const Link* from, ostream& os, double x, double y, double xscale, double yscale)	{
		DrawTree::LocalDraw(from,os,x,y,xscale,yscale);
		BubbleTree::LocalDraw(from,os,x,y,xscale,yscale);
		ChronoDrawTree::LocalDraw(from,os,x,y,xscale,yscale);
	}

	void PrepareDrawing()	{
		DrawTree::PrepareDrawing();
		BubbleTree::PrepareDrawing();
		ChronoDrawTree::PrepareDrawing();
	}

	void FinishDrawing(ostream& os, double xscale, double yscale)	{
		DrawTree::FinishDrawing(os,xscale,yscale);
		BubbleTree::FinishDrawing(os,xscale,yscale);
		ChronoDrawTree::FinishDrawing(os,xscale,yscale);
	}

	private:
	MeanChronogram* chrono;
	MeanExpNormTree* bubble;


};

#endif


