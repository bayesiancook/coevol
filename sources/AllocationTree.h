
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
        nshift = 0;
        if (data)   {
            cerr << "branch allocations: based on parsimony reconstruction using character data at the tips\n";
            MakeAllocMap(tree->GetRoot());
            cerr << '\n';
            cerr << "total number of shift events: " << nshift << '\n';
        }
        else    {
            cerr << "branch allocations: such as specified in tree file\n";
            RecursiveMakeBranchAllocations(tree->GetRoot());
        }
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
            if ((tmp < 0) || (tmp >= K))   {
                cerr << "error in MakeAllocMap: out of bound for taxon : " << from->GetNode()->GetName() << " trait : " << tmp << '\n';
                exit(1);
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
                        nshift++;
                    }
                }
                if (link->Next() != from)   {
                    cerr << ",";
                }
            }
            allocmap[from->GetBranch()] = ret;
            cerr << ")";
            cerr << ":" << ret;
        }
        return ret;
    }

	void RecursiveMakeBranchAllocations(const Link* from)	{
		
		if (! from->isRoot())	{
			int k = atoi(from->GetBranch()->GetName().c_str());
			if ((k<0) || (k >= K)) {
			    std::cerr << "error : allocation out of bound\n";
			    std::cerr << "k" << '\t' << "Ncond" << '\n';
			    exit(1);
			}
			allocmap[from->GetBranch()] = k;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMakeBranchAllocations(link->Out());
		}
	}

	int GetBranchAllocation(const Branch* branch)	{
        return allocmap[branch];
    }

	private:
	Tree* tree;
    ContinuousData* data;
	int K;	
    int offset;
    int nshift;
    map<const Branch*, int> allocmap;

};

#endif

