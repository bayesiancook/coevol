#ifndef PROBMODEL_H
#define PROBMODEL_H

#include <set>
#include <iostream>

#include "Move.hpp"


/// Probabilistic model
/**
 * A model is defined as a Graphical Model:
 * - a directed acyclic graph (DAG) with N nodes
 * - the nodes of the DAG are random variables (X_n)
 * - the joint probability law is given by
 *   \prod_n P(X_n | Pa(X_n))
 *   where  Pa(X_n) represents the set of all parents of node X_n
 *
 * A model implements MCMC:
 * - GetLogProb returns the log of the probability (or probability density) mentioned above
 * - Sample draws a model configuration from this joint probability
 * - Move resample the model's current configuration conditional on the data */
class ProbModel : public MCMC {
public:
  ProbModel();
  virtual ~ProbModel();

  /// obtain the set ("state") of all the nodes of the DAG by a recursive traversal from the root nodes to the tips
  void Register();
  void Register(DAGnode* var);
  void getDot();

  /// registers "var" among the root nodes (i.e. into the set "root")
  void RootRegister(DAGnode* var);

  // returns the log of the probability (or probability density) mentioned above
  virtual double GetLogProb() = 0;

  virtual void MakeScheduler() = 0;

  // resamples the model's current configuration conditional on the data
  // returns success rate
  inline virtual double Move(double tuning_modulator = 1) { return Move(tuning_modulator, 1, false, false);}

  virtual double Move(double tuning_modulator, int ncycle, bool verbose, bool check);

  void 	Corrupt();
  virtual double Update(bool check = false);

  bool 	CheckUpdateFlags();

  // save model configuration to stream
  virtual void ToStream(std::ostream& os) = 0;
  // get model configuration from stream
  virtual void FromStream(std::istream& is) = 0;

  // monitoring the run
  virtual void Trace(std::ostream&) {}
  virtual void TraceHeader(std::ostream&) {}
  virtual void Monitor(std::ostream& os, std::ostream& osdetail);
  // virtual void Monitor(std::ostream& os);

  virtual void Details(std::ostream&) {}
  virtual void test() {}

  std::set<DAGnode*> state;
  std::set<DAGnode*> root;

  MCScheduler scheduler;
};


#endif // PROBMODEL_H
