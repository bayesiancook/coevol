#ifndef ALLOCTREE_H
#define ALLOCTREE_H

#include "Tree.hpp"

class AllocationTree {
  public:
    AllocationTree(Tree* intree, int inK) {
        tree = intree;
        K = inK;
    }

    int GetBranchAllocation(const Branch* branch) {
        // if root : return 0
        if (branch == nullptr) {
            return 0;
        }

        int k = atoi(branch->GetName().c_str());

        /*
          if (k != 0)	{
          k=0;
          }
        */

        if (k >= K) {
            k = K - 1;
        }
        if (k < 0 || k >= K) {
            std::cerr << "error : allocation out of bound\n";
            std::cerr << "k" << '\t' << "K" << '\n';
            // std::cerr << tree->GetBranchVal(branch)->val() << '\t' << k <<
            //'\n';
            exit(1);
        }

        return k;
    }

    int GetComponentNumber() { return K; }

  private:
    Tree* tree;
    int K;
};

#endif
