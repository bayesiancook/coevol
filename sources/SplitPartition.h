
#ifndef SPLITPART_H
#define SPLITPART_H

#include "Partition.h"
#include "SplitTree.h"

class SplitBranchPartition : public BranchPartition	{

	public:

	SplitBranchPartition(BranchPartition* inmother, SplitTree* insplittree)	 : BranchPartition(insplittree) {
		mother = inmother;
		NComponent = mother->GetNComponent();
		RecursiveCreate(GetTree()->GetRoot());
	}

	SplitTree* GetSplitTree()	{
		SplitTree* tmp = dynamic_cast<SplitTree*>(GetTree());
		if (!tmp)	{
			cerr << "error in splitbranchpartition: null split tree pointer\n";
			exit(1);
		}
		return tmp;
	}

	protected:

	void RecursiveCreate(const Link* from)	{
		if (! from->isRoot())	{
			const Branch* branch = GetSplitTree()->GetSplitBranch(from->GetBranch())->GetMother();
			SetAlloc(from->GetBranch(),mother->GetAlloc(branch));
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out());
		}
	}

	BranchPartition* mother;

};

#endif

