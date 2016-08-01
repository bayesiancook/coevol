
#ifndef MEANDISCRETIZATION_H
#define	MEANDISCRETIZATION_H


#include "PureBrownianProcess.h"


class MeanDiscretization : public NewickTree {

private:
	Tree *tree;
	map<const Branch*,double> meandiscret;
	int size;



public:
	MeanDiscretization(Tree *intree) {
		tree = intree;
		Reset();
	}

	void Add(PureBrownianProcess *process) {
		RecursiveAdd(process, GetTree()->GetRoot());
		size++;
	}

	void RecursiveAdd(PureBrownianProcess *process, const Link* from) {
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveAdd(process, link->Out());
		int nSegments = process->GetBranchVal(link->GetBranch())->getNSegments();
		meandiscret[link->GetBranch()] += nSegments;
	}
	}

	Tree* GetTree() const {return tree;}
	const Link* GetRoot() const {return GetTree()->GetRoot();}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}
	void RecursiveReset(const Link* from)	{
		if(!from->isRoot())
			meandiscret[from->GetBranch()] = 0;

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}
	void RecursiveNormalise(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveNormalise(link->Out());
			meandiscret[link->GetBranch()] /= size;
		}
	}

	string GetNodeName(const Link* link) const {
		ostringstream s;
		if (link->isLeaf())	{
			s << link->GetNode()->GetName();
		}
		else	{
			s << "";
		}
		return s.str();
	}

	string GetBranchName(const Link* link) const {
		ostringstream s;
		s << GetMeanDiscret(link->GetBranch());
		return s.str();
	}

	double GetMeanDiscret(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = meandiscret.find(branch);
		return i->second;
	}
};




#endif	/* MEANDISCRETIZATION_H */

