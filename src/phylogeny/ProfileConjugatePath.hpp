#ifndef PROFILECONJUGATEPATH_H
#define PROFILECONJUGATEPATH_H

#include <map>
#include <utility>
#include "core/Conjugate.hpp"
#include "core/Move.hpp"
#include "phylogeny/AllocationTree.hpp"
#include "phylogeny/PhyloProcess.hpp"


class ProfilePathConjugate : public DSemiConjugatePrior<void> {
  public:
    ProfilePathConjugate(RandomSubMatrix* inmatrix) {
        matrix = inmatrix;
        Register(matrix);
        CreateSuffStat();
    }

    ~ProfilePathConjugate() override { DeleteSuffStat(); }

    void specialUpdate() override {}

    RandomSubMatrix* GetRandomSubMatrix() { return matrix; }

    virtual SubMatrix* GetMatrix() { return matrix; }
    virtual const double* GetStationary() { return matrix->GetStationary(); }

    int GetNstate() { return matrix->GetNstate(); }
    // StateSpace* GetStateSpace() {return matrix->GetStateSpace();}

    // suff stats:
    // for each state: total number of root occurrences
    // for each state: total waiting time
    // for each pair of state : total number of transitions

    void CreateSuffStat() {}

    void DeleteSuffStat() {}

    void ResetSufficientStatistic() override {
        rootcount.clear();
        paircount.clear();
        waitingtime.clear();
    }

    void SaveSufficientStatistic() override {
        std::cerr << "save sufficient\n";
        exit(1);
    }

    void RestoreSufficientStatistic() override {
        std::cerr << "restore sufficient\n";
        exit(1);
    }

    void IncrementRootCount(int state) { rootcount[state]++; }

    void AddWaitingTime(int state, double time) { waitingtime[state] += time; }

    void IncrementPairCount(int state1, int state2) {
        paircount[std::pair<int, int>(state1, state2)]++;
    }

    /*
      double SuffStatLogProb()	{
      double exptot = 1;
      double total = 0;
      // int totnsub = 0;

      const double* rootstat = GetStationary();
      for (std::map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
      double tmp = rootstat[i->first];
      for (int k=0; k<i->second; k++)	{
      exptot *= tmp;
      }
      // total += i->second * log(rootstat[i->first]);
      }

      SubMatrix& mat = *GetMatrix();
      for (std::map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
      total += i->second * mat(i->first,i->first);
      }
      for (std::map<std::pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end();
      i++)	{
      // total += i->second * log(mat(i->first.first, i->first.second));
      double tmp = mat(i->first.first, i->first.second);
      for (int k=0; k<i->second; k++)	{
      exptot *= tmp;
      }
      // totnsub += i->second;
      }
      total += log(exptot);

      return total;
      }
    */

    double SuffStatLogProb() override {
        double total = 0;
        int totnsub = 0;

        const double* rootstat = GetStationary();
        for (auto& i : rootcount) {
            total += i.second * log(rootstat[i.first]);
        }

        double totscalestat = 0;
        SubMatrix& mat = *GetMatrix();
        for (auto& i : waitingtime) {
            totscalestat += i.second * mat(i.first, i.first);
        }
        for (auto& i : paircount) {
            total += i.second * log(mat(i.first.first, i.first.second));
            totnsub += i.second;
        }
        total += totscalestat;

        return total;
    }

    double logProb() override {
        if (isActive()) {
            return SuffStatLogProb();
        }
        return 0;
    }

  protected:
    std::map<int, int> rootcount;
    std::map<std::pair<int, int>, int> paircount;
    std::map<int, double> waitingtime;

    RandomSubMatrix* matrix;
};


class ProfileConjugateRandomBranchSitePath : public virtual ConjugateSampling<void>,
                                             public virtual RandomBranchSitePath {
  public:
    ProfileConjugateRandomBranchSitePath(PhyloProcess* inprocess, ProfilePathConjugate* inpathconj,
                                         Var<PosReal>* inrate, Var<PosReal>* inlength)
        : RandomBranchSitePath(inprocess) {
        pathconj = inpathconj;

        if (pathconj == nullptr) {
            std::cerr << "error in ProfileConjugateRandomBranchSitePath\n";
            exit(1);
        }

        length = inlength;
        rate = inrate;
        matrix = pathconj->GetRandomSubMatrix();
        stationary = nullptr;
        active_flag = true;

        Register(length);
        Register(rate);
        /*
          Register(matrix);
          Register(stationary);
        */
        Register(pathconj);
        conjugate_up.insert(pathconj);
    }

  protected:
    void AddSufficientStatistic(SemiConjPrior* parent) override {
        if (parent != pathconj) {
            std::cerr << "error in ProfileConjugateRandomBranchSitePath::AddSufficientStatistic\n";
            exit(1);
        }

        if (isRoot()) {
            pathconj->IncrementRootCount(Init()->GetState());
        } else {
            Plink* link = Init();
            while (link != nullptr) {
                int state = link->GetState();
                pathconj->AddWaitingTime(state, GetAbsoluteTime(link));
                if (link != last) {
                    int newstate = link->Next()->GetState();
                    pathconj->IncrementPairCount(state, newstate);
                }
                link = link->Next();
            }
        }
    }

  private:
    ProfilePathConjugate* pathconj;
};

class ProfilePathConjugateArray {
  public:
    ProfilePathConjugateArray(int inNsite, int inK, RandomSubMatrix*** inmatrix) {
        Nsite = inNsite;
        K = inK;
        matrix = inmatrix;
        CreateArray();
    }

    ~ProfilePathConjugateArray() { DeleteArray(); }

    void CreateArray() {
        pathconjarray = new ProfilePathConjugate**[K];
        for (int k = 0; k < K; k++) {
            pathconjarray[k] = new ProfilePathConjugate*[Nsite];
            for (int i = 0; i < Nsite; i++) {
                pathconjarray[k][i] = new ProfilePathConjugate(matrix[k][i]);
            }
        }
    }

    void DeleteArray() {
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < Nsite; i++) {
                delete pathconjarray[k][i];
            }
            delete[] pathconjarray[k];
        }
        delete[] pathconjarray;
    }

    void ActivateSufficientStatistic() {
#ifdef _OPENMP
#pragma omp parallel for
#endif

        for (int i = 0; i < Nsite; i++) {
            for (int k = 0; k < K; k++) {
                pathconjarray[k][i]->ActivateSufficientStatistic();
            }
        }
    }

    void InactivateSufficientStatistic() {
#ifdef _OPENMP
#pragma omp parallel for
#endif

        for (int i = 0; i < Nsite; i++) {
            for (int k = 0; k < K; k++) {
                pathconjarray[k][i]->InactivateSufficientStatistic();
            }
        }
    }

    ProfilePathConjugate* GetProfilePathConjugate(int k, int site) {
        return pathconjarray[k][site];
    }

  protected:
    int Nsite;
    int K;
    ProfilePathConjugate*** pathconjarray;
    RandomSubMatrix*** matrix;
};


class ProfileConjugateSelectionPhyloProcess : public PhyloProcess {
  protected:
  public:
    ProfileConjugateSelectionPhyloProcess(AllocationTree* inalloctree, LengthTree* intree,
                                          SequenceAlignment* indata,
                                          ProfilePathConjugateArray* inpatharray)
        : PhyloProcess(intree, indata) {
        tree = intree;
        alloctree = inalloctree;
        patharray = inpatharray;
    }

    RandomBranchSitePath* CreateRandomBranchSitePath(const Link* link, int site) override {
        return new ProfileConjugateRandomBranchSitePath(
            this, patharray->GetProfilePathConjugate(
                      alloctree->GetBranchAllocation(link->GetBranch()), site),
            nullptr, tree->GetBranchVal(link->GetBranch()));
    }

  protected:
    ProfilePathConjugateArray* patharray;
    LengthTree* tree;
    AllocationTree* alloctree;
};


class ProfileConjugateMappingMove : public MCUpdate {
  public:
    ProfileConjugateMappingMove(PhyloProcess* inprocess, ProfilePathConjugateArray* inpathconjarray)
        : process(inprocess), pathconjarray(inpathconjarray) {}

    double Move(double /*tuning_modulator*/) override {
        pathconjarray->InactivateSufficientStatistic();
        process->Move(1);
        pathconjarray->ActivateSufficientStatistic();
        return 1;
    }

  protected:
    PhyloProcess* process;
    ProfilePathConjugateArray* pathconjarray;
};

class ProfileConjugateMove : public MCUpdate {
  public:
    ProfileConjugateMove(ProfilePathConjugateArray* inpathconjarray, bool intoggle)
        : pathconjarray(inpathconjarray), toggle(intoggle) {}

    double Move(double /*tuning_modulator*/) override {
        if (toggle) {
            pathconjarray->ActivateSufficientStatistic();
        } else {
            pathconjarray->InactivateSufficientStatistic();
        }
        return 1;
    }

  protected:
    ProfilePathConjugateArray* pathconjarray;
    bool toggle;
};

#endif
