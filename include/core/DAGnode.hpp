#ifndef RNODE_H
#define RNODE_H

#include <set>
#include <map>
#include <cmath>// FIXME possibly better to remove and move log to cpp

#include "Random.hpp"
#include "MCMC.hpp"

class ProbModel;


class DAGnode {
public:
  static bool initmode;

  DAGnode() : flag(false), name("") {}
  virtual ~DAGnode();

  // Getters
  inline int GetChildNumber() { return down.size(); }
  inline void SetName(std::string inname) { name = inname; }
  inline std::string GetName() { return name; }
  inline bool isUpdated() { return flag; }
  inline virtual bool isValueUpdated() { return flag; }

  virtual void Register(DAGnode* in);
  void RecursiveRegister(ProbModel* model);

  bool CheckUpdateFlags();
  int GetChildrenNumber();

  // methods associated to the Graphical Model aspect
  virtual void Corrupt(bool bk) = 0;
  virtual double Update() = 0;
  virtual void Restore() = 0;

  virtual void NotifyCorrupt(bool bk) = 0;
  virtual double NotifyUpdate() = 0;
  virtual void NotifyRestore() = 0;

  virtual double FullUpdate(bool check = false) = 0;
  virtual void FullCorrupt(std::map<DAGnode*,int>& m) = 0;

protected:
  std::set<DAGnode*> up;
  std::set<DAGnode*> down;

  void Detach();
  void DeregisterFrom(DAGnode* parent);

  bool flag;
  std::string name;

};


class Rnode : public virtual DAGnode, public MH {
public:
  Rnode() : logprob(0), bklogprob(0), value_updated (false) {}

  // Metropolis Hastings move
  virtual double Move(double tuning = 1) ;

  virtual void Corrupt(bool bk);
  virtual double Update();
  virtual void Restore();

  virtual void NotifyCorrupt(bool bk);
  virtual double NotifyUpdate();
  virtual void NotifyRestore();

  virtual double FullUpdate(bool check = false);
  virtual void FullCorrupt(std::map<DAGnode*,int>& m);
  virtual double localUpdate();

  // Getters
  inline virtual double GetFastLogProb() { return (flag ? logprob : logProb()); }
  inline virtual double GetLogProb() { return logProb(); }
  inline bool isValueUpdated() { return value_updated; }

protected:
  virtual double logProb() = 0;

  virtual void RestoreBackup() {}

  virtual void localRestore();
  virtual void localCorrupt(bool bk);

  double logprob;
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

  virtual void localRestore();
  virtual void localCorrupt(bool bk);

};


class Mnode : public virtual DAGnode {
public:
  Mnode(bool inflag = true) { flag = inflag; }

  virtual void Corrupt(bool bk);
  virtual double Update();
  virtual void Restore();

protected:
  inline virtual void NotifyCorrupt(bool) {}
  inline virtual double NotifyUpdate() { return 0; }
  inline virtual void NotifyRestore() {}

  inline virtual double FullUpdate(bool) { return 0; }
  inline virtual void FullCorrupt(std::map<DAGnode*,int>&) {}

};

#endif
