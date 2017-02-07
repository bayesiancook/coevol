#ifndef RNODE_H
#define RNODE_H

#include <set>
#include <map>
#include <cmath>// FIXME possibly better to remove and move log to cpp

#include "Random.hpp"
#include "MCMC.hpp"

class DAGnode;
class ProbModel;


class DAGnode {
public:
  static bool initmode;

  friend class ProbModel;
  friend class Rnode;
  friend class Dnode;

  DAGnode() : flag(false), name("") {}
  virtual  ~DAGnode();

  // methods associated to the DAG aspect
  int GetChildNumber() {return down.size();}

  void SetName(std::string inname) { name = inname;}
  std::string GetName() {return name;}

  virtual void Register(DAGnode* in);
  // virtual void RegisterChild(DAGnode* in) {in->Register(this);}

  void  RecursiveRegister(ProbModel* model);

  bool CheckUpdateFlags();
  int GetChildrenNumber();

  // methods associated to the Graphical Model aspect
  virtual void Corrupt(bool bk) = 0;
  virtual double Update() = 0;
  virtual void Restore() = 0;

  bool  isUpdated() {return flag;}
  virtual bool  isValueUpdated() {return flag;}

  virtual void NotifyCorrupt(bool bk) = 0;
  virtual double NotifyUpdate() = 0;
  virtual void NotifyRestore() = 0;

  // virtual void Initialise() = 0;
  virtual double FullUpdate(bool check = false) = 0;
  virtual void FullCorrupt(std::map<DAGnode*,int>& m) = 0;

  // template <class V> friend class PointerKeeper;

protected:
// public: // (VL): there was both a protected and a public ; the public does not seem necessary
  std::set<DAGnode*>  up;
  std::set<DAGnode*>  down;

  void  Detach();
  void  DeregisterFrom(DAGnode* parent);

  bool flag;
  std::string name;
};


class Rnode : public virtual DAGnode, public MH {
public:
  Rnode() : logprob(0), bklogprob(0), value_updated (false) {}

  virtual double GetFastLogProb() {return (flag ? logprob : logProb());}

  virtual double GetLogProb() { return logProb();}
  /*
    logprob = logProb();
    return logprob;
    }
  */

  // Metropolis Hastings move
  virtual double Move(double tuning = 1) {
    if (! isClamped()) {
      Corrupt(true);
      double logHastings = ProposeMove(tuning);
      double deltaLogProb = Update();
      double logRatio = deltaLogProb + logHastings;
      bool accepted = (log(Random::Uniform()) < logRatio);
      if (! accepted) {
        Corrupt(false);
        Restore();
      }
      return (double) accepted;
    }
    return 1;
  }

  virtual void Corrupt(bool bk);
  virtual double Update();
  virtual void Restore();

  virtual void NotifyCorrupt(bool bk);
  virtual double NotifyUpdate();
  virtual void NotifyRestore();

  virtual double FullUpdate(bool check = false);
  virtual void FullCorrupt(std::map<DAGnode*,int>& m);

  // virtual void  Initialise();

  /*
    virtual void Sample() {
    // this version is a bit contrieved, and computationally far from optimal
    // but this is necessary as long as corruption is dealt with in a forward manner
    Corrupt(false);
    if (! isClamped()) {
    drawSample();
    }
    Update();
    }
  */

  bool  isValueUpdated() {return value_updated;}

  virtual double localUpdate();

protected:
  virtual double logProb() = 0;

  virtual void RestoreBackup() {}

  virtual void localRestore();
  virtual void localCorrupt(bool bk);

  double  logprob;
  double bklogprob;
  bool value_updated;
};

class Dnode : public virtual DAGnode {
public:
  virtual void Corrupt(bool bk);
  virtual double Update();
  virtual void Restore();

  virtual double localUpdate();

  virtual void specialUpdate() = 0;

  virtual void FullCorrupt(std::map<DAGnode*,int>& m);
  virtual double FullUpdate(bool check = false);

protected:
  virtual void NotifyCorrupt(bool bk);
  virtual double NotifyUpdate();
  virtual void NotifyRestore();

  // virtual void  Initialise();

  virtual void localRestore();
  virtual void localCorrupt(bool bk);

};


class Mnode : public virtual DAGnode {
public:
  Mnode(bool inflag = true) {
    flag = inflag;
  }

  virtual void Corrupt(bool bk);
  virtual double Update();
  virtual void Restore();

protected:
  virtual void NotifyCorrupt(bool) {}
  virtual double NotifyUpdate() {return 0;}
  virtual void NotifyRestore() {}

  virtual double FullUpdate(bool) {return 0;}
  virtual void FullCorrupt(std::map<DAGnode*,int>&) {}

  // virtual void  Initialise() {}
};

#endif
