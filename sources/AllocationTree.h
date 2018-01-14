
#ifndef ALLOCTREE_H
#define ALLOCTREE_H

#include "ContinuousData.h"

class AllocationTree {

	public:

	AllocationTree(Tree* intree, ContinuousData* indata, int inK, int inoffset)	{
		tree = intree;
		K = inK;
        offset = inoffset;
        data = indata;
        MakeAllocMap(tree->GetRoot());
        cerr << '\n';
	}

    int MakeAllocMap(const Link* from) {

        int ret = -1;
        if (from->isLeaf()) {
			int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax == -1)	{
                cerr << "error: did not find taxon: " << from->GetNode()->GetName() << '\n';
                exit(1);
            }
            int tmp = int(data->GetState(tax, 0)) + offset;
            if (tmp >= K)   {
                tmp = K-1;
            }
            allocmap[from->GetBranch()] = tmp;
            ret = tmp;
            cerr << from->GetNode()->GetName() << ":" << tmp;
        }
        else    {

            cerr << "(";
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                int tmp = MakeAllocMap(link->Out());
                if (ret == -1)  {
                    ret = tmp;
                }
                else    {
                    if (tmp != ret) {
                        ret = 0;
                    }
                }
                if (link->Next() != from)   {
                    cerr << ",";
                }
            }
            allocmap[from->GetBranch()] = ret;
            cerr << ":" << ret;
            cerr << ")";
        }
        return ret;
    }

	int GetBranchAllocation(const Branch* branch)	{
        return allocmap[branch];
    }

	private:
	Tree* tree;
    ContinuousData* data;
	int K;	
    int offset;
    map<const Branch*, int> allocmap;

};

#endif

