#ifndef RNODE_H
#define RNODE_H

#include <map>
#include <set>
#include <string>
#include "MCMC.hpp"

class ProbModel;

class DAGnode {
  public:
    static bool initmode;

    DAGnode() : updateFlag(false), dotNodeFlag(false), dotVertexFlag(false), name("") {}
    virtual ~DAGnode();

    // Getters
    inline int GetChildNumber() { return down.size(); }
    inline void SetName(std::string inname) { name = inname; }
    inline std::string GetName() { return name; }
    inline bool isUpdated() { return updateFlag; }
    inline virtual bool isValueUpdated() { return updateFlag; }

    virtual void Register(DAGnode *parent);
    void RecursiveRegister(ProbModel *model);
    virtual std::set<std::string> getDotNodes();
    std::set<std::string> getDotVertices();

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
    virtual void FullCorrupt(std::map<DAGnode *, int> &m) = 0;

  protected:
    std::set<DAGnode *> up;
    std::set<DAGnode *> down;

    void Detach();
    void DeregisterFrom(DAGnode *parent);

    // VL refactoring
    bool parentsUpdated();

    bool updateFlag;
    bool dotNodeFlag;
    bool dotVertexFlag;
    std::string name;
};

class Rnode : public virtual DAGnode, public MH {
  public:
    Rnode() : logprob(0), bklogprob(0), value_updated(false) {}

    // Metropolis Hastings move
    double Move(double tuning = 1) override;

    std::set<std::string> getDotNodes() override;

    void Corrupt(bool bk) override;
    double Update() override;
    void Restore() override;

    void NotifyCorrupt(bool bk) override;
    double NotifyUpdate() override;
    void NotifyRestore() override;

    double FullUpdate(bool check = false) override;
    void FullCorrupt(std::map<DAGnode *, int> &m) override;
    virtual double localUpdate();

    // Getters
    inline virtual double GetFastLogProb() { return (updateFlag ? logprob : logProb()); }
    inline double GetLogProb() override { return logProb(); }
    inline bool isValueUpdated() override { return value_updated; }

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
    void Corrupt(bool bk) override;
    double Update() override;
    void Restore() override;

    virtual double localUpdate();
    virtual void specialUpdate() = 0;

    void FullCorrupt(std::map<DAGnode *, int> &m) override;
    double FullUpdate(bool check = false) override;

  protected:
    void NotifyCorrupt(bool bk) override;
    double NotifyUpdate() override;
    void NotifyRestore() override;

    virtual void localRestore();
    virtual void localCorrupt(bool bk);
};

class Mnode : public virtual DAGnode {
  public:
    Mnode(bool inflag = true) { updateFlag = inflag; }

    void Corrupt(bool bk) override;
    double Update() override;
    void Restore() override;

  protected:
    inline void NotifyCorrupt(bool /*bk*/) override {}
    inline double NotifyUpdate() override { return 0; }
    inline void NotifyRestore() override {}

    inline double FullUpdate(bool /*check*/) override { return 0; }
    inline void FullCorrupt(std::map<DAGnode *, int> & /*m*/) override {}
};

#endif