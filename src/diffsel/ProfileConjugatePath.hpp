#ifndef PROFILECONJUGATEPATH_H
#define PROFILECONJUGATEPATH_H

#include "RandomBranchSitePath.hpp"
#include "Conjugate.hpp"
#include "RandomSubMatrix.hpp"
#include "PhyloProcess.hpp"

#include <utility>
#include <map>


class ProfilePathConjugate : public DSemiConjugatePrior<void>	{

public:

  ProfilePathConjugate(RandomSubMatrix* inmatrix) {
    matrix = inmatrix;
    Register(matrix);
    CreateSuffStat();

  }

  ~ProfilePathConjugate() {
    DeleteSuffStat();
  }

  void specialUpdate()	{
  }

  RandomSubMatrix* GetRandomSubMatrix() {return matrix;}

  virtual SubMatrix*    GetMatrix()		{return matrix;}
  virtual const double*     GetStationary()		{return matrix->GetStationary();}

  int GetNstate() {return matrix->GetNstate();}
  // StateSpace* GetStateSpace() {return matrix->GetStateSpace();}

  // suff stats:
  // for each state: total number of root occurrences
  // for each state: total waiting time
  // for each pair of state : total number of transitions

  void CreateSuffStat()	{
  }

  void DeleteSuffStat()	{
  }

  void ResetSufficientStatistic()	{
    rootcount.clear();
    paircount.clear();
    waitingtime.clear();
  }

  void SaveSufficientStatistic()	{
    cerr << "save sufficient\n";
    exit(1);
  }

  void RestoreSufficientStatistic()	{
    cerr << "restore sufficient\n";
    exit(1);
  }

  void IncrementRootCount(int state)	{
    rootcount[state]++;
  }

  void AddWaitingTime(int state, double time)	{
    waitingtime[state] += time;
  }

  void IncrementPairCount(int state1, int state2)	{
    paircount[ pair<int,int>(state1,state2)]++;
  }

  /*
    double SuffStatLogProb()	{
    double exptot = 1;
    double total = 0;
    // int totnsub = 0;

    const double* rootstat = GetStationary();
    for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
    double tmp = rootstat[i->first];
    for (int k=0; k<i->second; k++)	{
    exptot *= tmp;
    }
    // total += i->second * log(rootstat[i->first]);
    }

    SubMatrix& mat = *GetMatrix();
    for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
    total += i->second * mat(i->first,i->first);
    }
    for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
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

  double SuffStatLogProb()	{
    double total = 0;
    int totnsub = 0;

    const double* rootstat = GetStationary();
    for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
      total += i->second * log(rootstat[i->first]);
    }

    double totscalestat = 0;
    SubMatrix& mat = *GetMatrix();
    for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
      totscalestat += i->second * mat(i->first,i->first);
    }
    for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
      total += i->second * log(mat(i->first.first, i->first.second));
      totnsub += i->second;
    }
    total += totscalestat;

    return total;
  }

  double logProb()	{
    if (isActive())	{
      return SuffStatLogProb();
    }
    return 0;
  }

protected:

  map<int,int> rootcount;
  map< pair<int,int>, int> paircount;
  map<int,double> waitingtime;

  RandomSubMatrix* matrix;
};


class ProfileConjugateRandomBranchSitePath : public virtual ConjugateSampling<void>, public virtual RandomBranchSitePath	{

public:

  ProfileConjugateRandomBranchSitePath(PhyloProcess* inprocess, ProfilePathConjugate* inpathconj, Var<PosReal>* inrate, Var<PosReal>* inlength) : RandomBranchSitePath(inprocess)	{

    pathconj = inpathconj;

    if (!pathconj)	{
      cerr << "error in ProfileConjugateRandomBranchSitePath\n";
      exit(1);
    }

    length = inlength;
    rate = inrate;
    matrix = pathconj->GetRandomSubMatrix();
    stationary = 0;
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

  void AddSufficientStatistic(SemiConjPrior* parent) {
    if (parent != pathconj)	{
      cerr << "error in ProfileConjugateRandomBranchSitePath::AddSufficientStatistic\n";
      exit(1);
    }

    if (isRoot())	{
      pathconj->IncrementRootCount(Init()->GetState());
    }
    else	{
      Plink* link = Init();
      while (link)	{
        int state = link->GetState();
        pathconj->AddWaitingTime(state,GetAbsoluteTime(link));
        if (link != last)	{
          int newstate = link->Next()->GetState();
          pathconj->IncrementPairCount(state,newstate);
        }
        link = link->Next();
      }
    }
  }

private:

  ProfilePathConjugate* pathconj;
};

class ProfilePathConjugateArray	{

public:

  ProfilePathConjugateArray(int inNsite, int inK, RandomSubMatrix*** inmatrix) {
    Nsite = inNsite;
    K = inK;
    matrix = inmatrix;
    CreateArray();
  }

  ~ProfilePathConjugateArray()	{
    DeleteArray();
  }

  void CreateArray()	{
    pathconjarray = new ProfilePathConjugate**[K];
    for (int k=0; k<K; k++)	{
      pathconjarray[k] = new ProfilePathConjugate*[Nsite];
      for (int i=0; i<Nsite; i++)	{
        pathconjarray[k][i] = new ProfilePathConjugate(matrix[k][i]);
      }
    }

  }

  void DeleteArray()	{
    for (int k=0; k<K; k++)	{
      for (int i=0; i<Nsite; i++)	{
        delete pathconjarray[k][i];
      }
      delete[] pathconjarray[k];
    }
    delete[] pathconjarray;
  }

  void ActivateSufficientStatistic()	{

#ifdef _OPENMP
#pragma omp parallel for
#endif

    for (int i=0; i<Nsite; i++)	{
      for (int k=0; k<K; k++)	{
        pathconjarray[k][i]->ActivateSufficientStatistic();
      }
    }
  }

  void InactivateSufficientStatistic()	{

#ifdef _OPENMP
#pragma omp parallel for
#endif

    for (int i=0; i<Nsite; i++)	{
      for (int k=0; k<K; k++)	{
        pathconjarray[k][i]->InactivateSufficientStatistic();
      }
    }
  }

  ProfilePathConjugate* GetProfilePathConjugate(int k, int site)	{
    return pathconjarray[k][site];
  }

protected:

  int Nsite;
  int K;
  ProfilePathConjugate*** pathconjarray;
  RandomSubMatrix*** matrix;
};


class ProfileConjugateSelectionPhyloProcess : public PhyloProcess	{


protected:

public:

  ProfileConjugateSelectionPhyloProcess(AllocationTree* inalloctree, LengthTree* intree, SequenceAlignment* indata, ProfilePathConjugateArray* inpatharray) : PhyloProcess(intree,indata)	{
    tree = intree;
    alloctree = inalloctree;
    patharray = inpatharray;
  }

  virtual RandomBranchSitePath*   CreateRandomBranchSitePath(const Link* link, int site)	{
    return  new ProfileConjugateRandomBranchSitePath(this, patharray->GetProfilePathConjugate(alloctree->GetBranchAllocation(link->GetBranch()),site), 0, tree->GetBranchVal(link->GetBranch()));
  }

protected:
  ProfilePathConjugateArray* patharray;
  LengthTree* tree;
  AllocationTree* alloctree;

};


class ProfileConjugateMappingMove : public MCUpdate	{

public:

  ProfileConjugateMappingMove(PhyloProcess* inprocess, ProfilePathConjugateArray* inpathconjarray) : process(inprocess),  pathconjarray(inpathconjarray) {}

  double Move(double tuning_modulator=1)	{
    pathconjarray->InactivateSufficientStatistic();
    process->Move(1);
    pathconjarray->ActivateSufficientStatistic();
    return 1;
  }

protected:

  PhyloProcess* process;
  ProfilePathConjugateArray* pathconjarray;
};

class ProfileConjugateMove : public MCUpdate	{

public:

  ProfileConjugateMove(ProfilePathConjugateArray* inpathconjarray, bool intoggle) : pathconjarray(inpathconjarray), toggle(intoggle) {}

  double Move(double tuning_modulator=1)	{
    if (toggle)	{
      pathconjarray->ActivateSufficientStatistic();
    }
    else	{
      pathconjarray->InactivateSufficientStatistic();
    }
    return 1;
  }

protected:

  ProfilePathConjugateArray* pathconjarray;
  bool toggle;
};

#endif
