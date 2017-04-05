#ifndef BRANCHSITESUBPROCESS_H
#define BRANCHSITESUBPROCESS_H

#include "BranchSitePath.hpp"
#include "BranchSiteSubstitutionProcess.hpp"
#include "RandomSubMatrix.hpp"
#include "core/BaseType.hpp"
#include "core/Var.hpp"

// RandomBranchSitePath
// - represents a branch- and site-specific substitution process (the
// BranchSiteSubstitutionProcess
// superclass)
// - maintains a realisation of the substitution history along a branch/site (a
// path) conditional on
// the data (the BranchSitePath super class)
//
// - it is an abstract class (BranchSiteSubstitutionProcess has many pure
// virtual methods that need
// to be implemented in subclasses)
// - it is a random node of the model's graph (the Rnode superclass)
// - subclasses will be responsible for making the required connections between
// this Rnode and the
// relevant parameters of the process
//   on a model-specific basis

class PhyloProcess;

class RandomBranchSitePath : public virtual Rnode,
                             public BranchSitePath,
                             public BranchSiteSubstitutionProcess {
  public:
    RandomBranchSitePath(PhyloProcess *inprocess)
        : rate(nullptr),
          length(nullptr),
          matrix(nullptr),
          stationary(nullptr),
          myprocess(inprocess),
          propmatrix(nullptr),
          matswap(false) {
        SetName("path");
    }

    RandomBranchSitePath(PhyloProcess *inprocess, Var<PosReal> *inlength, Var<PosReal> *inrate,
                         RandomSubMatrix *inmatrix, Var<Profile> *instationary)
        : myprocess(inprocess), propmatrix(nullptr), matswap(false) {
        SetName("path");
        length = inlength;
        rate = inrate;
        matrix = inmatrix;
        stationary = instationary;
        active_flag = true;
        if (stationary != nullptr) {
            std::cerr << "in random branch site path: stat is non null\n";
            exit(1);
        }

        if ((matrix == nullptr) && (stationary == nullptr)) {
            std::cerr << "error in RandomBranchSitePath: should specify a matrix or "
                         "a stationary\n";
            exit(1);
        }

        Register(length);
        Register(rate);
        Register(matrix);
        Register(stationary);
        // Sample();
    }

    RandomBranchSitePath(PhyloProcess *inprocess, RandomTransitionMatrix *inmatrix,
                         Var<Profile> *instationary, Var<PosReal> *inlength)
        : myprocess(inprocess) {
        SetName("path");
        propmatrix = nullptr;
        matrix = nullptr;
        length = inlength;
        rate = nullptr;
        matswap = false;
        transitionmatrix = inmatrix;
        stationary = instationary;
        Register(transitionmatrix);
        if (stationary != nullptr) {
            std::cerr << "non null stat in random branch site path\n";
            exit(1);
        }
        // Register(stationary);
        // Register(length);
    }

    PhyloProcess *GetPhyloProcess() { return myprocess; }

    bool SampleBranchMapping();

    virtual AbstractTransitionMatrix *GetTransitionMatrix() { return transitionmatrix; }

    SubMatrix *GetSubMatrix() override {
        if (matswap) {
            if (propmatrix == nullptr) {
                std::cerr << "error in random branch site path: null prop matrix\n";
                exit(1);
            }
            return propmatrix;
        }
        if (matrix == nullptr) {
            std::cerr << "error in random branch site path: getsubmatrix\n";
            exit(1);
        }
        return matrix;
    }

    virtual void SwapMatrices() {
        if ((!matswap) && (propmatrix == nullptr)) {
            std::cerr << "error in random branch site path swap matrix: null prop matrix\n";
            exit(1);
        }
        matswap = !matswap;
        /*
          RandomSubMatrix* tmp = matrix;
          matrix = propmatrix;
          propmatrix = tmp;
        */
    }

    // void SetProposalMatrix(EmpiricalSubMatrix* inmatrix) {
    void SetProposalMatrix(SubMatrix *inmatrix) { propmatrix = inmatrix; }

    void BackwardPropagate(const double *down, double *up) {
        if (SampleBranchMapping()) {
            BranchSiteSubstitutionProcess::BackwardPropagate(down, up);
        } else {
            transitionmatrix->BackwardPropagate(down, up, 0);
        }
    }

    void ForwardPropagate(const double *up, double *down) {
        if (SampleBranchMapping()) {
            BranchSiteSubstitutionProcess::ForwardPropagate(up, down);
        } else {
            transitionmatrix->ForwardPropagate(up, down, 0);
        }
    }

    void GetFiniteTimeTransitionProb(int state, double *aux) {
        if (SampleBranchMapping()) {
            BranchSiteSubstitutionProcess::GetFiniteTimeTransitionProb(state, aux);
        } else {
            TransitionGetFiniteTimeTransitionProb(state, aux);
        }
    }

    int DrawStationary() {
        if (SampleBranchMapping()) {
            return BranchSiteSubstitutionProcess::DrawStationary();
        }
        return TransitionDrawStationary();
    }

    int DrawFiniteTime(int state) {
        if (SampleBranchMapping()) {
            return BranchSiteSubstitutionProcess::DrawFiniteTime(state);
        }
        return TransitionDrawFiniteTime(state);
    }

    void TransitionGetFiniteTimeTransitionProb(int state, double *aux) {
        const double *p = transitionmatrix->GetRow(state);
        for (int i = 0; i < GetNstate(); i++) {
            aux[i] = p[i];
        }
    }

    int TransitionDrawStationary() {
        const double *p = transitionmatrix->GetStationary();
        int newstate = Random::DrawFromDiscreteDistribution(p, GetNstate());
        return newstate;
    }

    int TransitionDrawFiniteTime(int state) {
        const double *p = transitionmatrix->GetRow(state);
        int newstate = Random::DrawFromDiscreteDistribution(p, GetNstate());
        return newstate;
    }

    int GetNstate() override {
        if (matrix != nullptr) {
            return GetSubMatrix()->GetNstate();
        }
        if (transitionmatrix != nullptr) {
            return transitionmatrix->GetNstate();
        } else {
            std::cerr << "error in RandomBranchSitePath: GetNstate\n";
            exit(1);
        }
    }

    double GetRate() override { return rate != nullptr ? ((double)rate->val()) : 1; }

    const double *GetStationary() override {
        if (stationary != nullptr) {
            std::cerr << "error : non null stationary in log prob path\n";
            exit(1);
            return stationary->GetArray();
        }
        if (SampleBranchMapping()) {
            return GetSubMatrix()->GetStationary();
        }
        return transitionmatrix->GetStationary();
    }
    // return stationary ? stationary->GetArray() : matrix->GetStationary();}

    double GetTime() override {
        //    return return length ? ((double) length->val()) : 0;}
        if (length != nullptr) {
            if (std::isnan(((double)(length->val())))) {
                std::cerr << "length is nan\n";
            }
            return length->val();
        }
        return 0;
    }

    bool isRoot() { return (length == nullptr); }

    bool isActivated() { return active_flag; }

    void SetActivated(bool inflag) { active_flag = inflag; }

    StateSpace *GetStateSpace();

    void SetUp(RandomBranchSitePath *inup);

    double GetTotalTime() override { return GetTime(); }
    void SetTotalTime(double intime) override {
        // what about corruption ?
        if (length == nullptr) {
            std::cerr << "error in RandomBranchSitePath::SetTotalTime\n";
            exit(1);
        }
        length->setval(intime);
    }

    int GetMaxTrial();

    virtual void Resample() {
        if (!SampleBranchMapping()) {
            std::cerr << "error in random branch site path : resample map called\n";
            exit(1);
        }
        if (isRoot()) {
            std::cerr << "error in RandomBranchSitePath::Resample : called on root\n";
            exit(1);
        } else {
            if (GetTime() == 0) {
                std::cerr << "error in resample : null bl\n";
                exit(1);
            }
            bool ok = ResampleAcceptReject(GetMaxTrial());
            if (!ok) {
                ResampleUniformized();
            }
            /*
              if (propmatrix) {
              ResampleUniformized();
              }
              else {
              bool ok = ResampleAcceptReject(GetMaxTrial());
              if (! ok) {
              ResampleUniformized();
              }
              }
            */
        }
    }

    double Move(double tuning) override {
        if (propmatrix == nullptr) {
            Resample();
            return 1;
        }
        return Rnode::Move(tuning);
    }

    void localRestore() override {
        if (propmatrix != nullptr) {
            stateup = bkstateup;
            statedown = bkstatedown;
            RestorePath();
            logprob = bklogprob;
            /*
              if (logprob != logProb()) {
              std::cerr << "error in local restore\n";
              std::cerr << logprob - logProb() << '\n';
              exit(1);
              }
            */
            updateFlag = true;
        } else {
            Rnode::localRestore();
        }
    }

    void localCorrupt(bool bk) override {
        if (propmatrix != nullptr) {
            if (bk) {
                bkstateup = stateup;
                bkstatedown = statedown;
                BackupPath();
                bklogprob = logprob;
            }
            updateFlag = false;
        } else {
            Rnode::localCorrupt(bk);
        }
    }

    void drawSample() override {
        std::cerr << "in random branch site path drawsample\n";
        exit(1);
        Resample();
    }

    double logProb() override;
    double PathLogProb();
    double NoPathLogProb();

    void ResampleRoot(int state) {
        stateup = statedown = state;
        Reset(state);
        localUpdate();
    }

    void Resample(int instateup, int instatedown) {
        stateup = instateup;
        statedown = instatedown;
        if (SampleBranchMapping()) {
            if (propmatrix != nullptr) {
                SwapMatrices();
                Resample();
                SwapMatrices();
            } else {
                Resample();
            }
        } else {
            Reset(stateup);
        }
        localUpdate();
    }

    void SetMatrix(RandomSubMatrix *inmatrix) { matrix = inmatrix; }

    double ProposeMove(double /*tuning*/) override {
        std::cerr << "error : in random branch site path propose move(tuning)\n";
        exit(1);
    }

    double ProposeMove(int instateup, int instatedown) {
        double logbefore = logProb();
        stateup = instateup;
        statedown = instatedown;
        Resample();
        double logafter = logProb();
        return logbefore - logafter;
    }

    double ProposeMoveRoot(int instate) {
        double logbefore = logProb();
        stateup = statedown = instate;
        Reset(instate);
        double logafter = logProb();
        return logbefore - logafter;
    }

  protected:
    bool ResampleAcceptReject(int maxtrial);
    void ResampleUniformized();
    // double   RecordResampleUniformizedLogProb();

    // RandomBranchSitePath* pathup;
    Var<PosReal> *rate;
    Var<PosReal> *length;

    RandomSubMatrix *matrix;
    RandomTransitionMatrix *transitionmatrix;

  protected:
    Var<Profile> *stationary;

    int stateup;
    int statedown;
    int bkstateup;
    int bkstatedown;
    PhyloProcess *myprocess;
    SubMatrix *propmatrix;
    // EmpiricalSubMatrix* propmatrix;
    bool active_flag;

    bool matswap;
    /*
      int uninsub;
      vector<int> unistate;
      vector<double> unitime;
    */
};

#endif  // RANDOMSITEPATH_H