#ifndef ONEMATRIXPHYLOPROCESS_H
#define ONEMATRIXPHYLOPROCESS_H

#include "PhyloProcess.hpp"
// #include "RandomSubMatrix.hpp"
// #include "diffsel/GTRSubMatrix.hpp"
// #include "RandomTypes.hpp"
#include "ValArray.hpp"

// all branches and all sites are described by one single General Time Reversible (GTR) substitution
// process

class OneMatrixPhyloProcess : public PhyloProcess {
  protected:
  public:
    OneMatrixPhyloProcess(LengthTree* intree, RandomSubMatrix* randmatrix,
                          SequenceAlignment* indata)
        : PhyloProcess(intree, indata) {
        matrix = randmatrix;
    }

    // the CreateRandomBranchSitePath function tells the PhyloProcess class
    // how to create the substitution process (and the associated substitution path)
    // for a given branch (accessible through link), and a given site
    //
    RandomBranchSitePath* CreateRandomBranchSitePath(const Link* link, int /*site*/) override {
        return new RandomBranchSitePath(this, tree->GetBranchLength(link->GetBranch()), nullptr,
                                        GetMatrix(), nullptr);
    }

    RandomSubMatrix* GetMatrix() { return matrix; }

  protected:
    RandomSubMatrix* matrix;
};

class SiteMatrixPhyloProcess : public PhyloProcess {
  protected:
  public:
    SiteMatrixPhyloProcess(LengthTree* intree, RandomSubMatrix** randmatrix,
                           SequenceAlignment* indata)
        : PhyloProcess(intree, indata) {
        matrix = randmatrix;
    }

    // the CreateRandomBranchSitePath function tells the PhyloProcess class
    // how to create the substitution process (and the associated substitution path)
    // for a given branch (accessible through link), and a given site
    //
    RandomBranchSitePath* CreateRandomBranchSitePath(const Link* link, int site) override {
        return new RandomBranchSitePath(this, tree->GetBranchLength(link->GetBranch()), nullptr,
                                        GetMatrix(site), nullptr);
    }

    RandomSubMatrix* GetMatrix(int site) { return matrix[site]; }

  protected:
    RandomSubMatrix** matrix;
};

class OneMatrixRASPhyloProcess : public OneMatrixPhyloProcess {
  protected:
  public:
    OneMatrixRASPhyloProcess(LengthTree* intree, VarArray<PosReal>* inrate,
                             RandomSubMatrix* randmatrix, SequenceAlignment* indata)
        : OneMatrixPhyloProcess(intree, randmatrix, indata) {
        rate = inrate;
    }

    RandomBranchSitePath* CreateRandomBranchSitePath(const Link* link, int site) override {
        return new RandomBranchSitePath(this, tree->GetBranchLength(link->GetBranch()),
                                        GetRate(site), matrix, nullptr);
    }

    Var<PosReal>* GetRate(int site) {
        if (rate == nullptr) {
            return nullptr;
        }
        return rate->GetVal(site);
    }

  protected:
    VarArray<PosReal>* rate;
};


#endif
