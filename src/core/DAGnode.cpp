#include <algorithm>
#include <cstdio>
using namespace std;

#include "core/ProbModel.hpp"
#include "core/DAGnode.hpp"
#include "core/Exception.hpp"


bool DAGnode::initmode = true; // (VL) FIXME


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * DAGnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
const double MCMC::MAXDIFF = 1e-4;

DAGnode::~DAGnode() {
  Detach();
}

int DAGnode::GetChildrenNumber() {
  int tot = 0;
  for (auto i=down.begin(); i!=down.end(); i++) {
    tot++;
  }
  return tot;
}

void DAGnode::Detach() {
  while (up.size()) {
    DeregisterFrom(*up.begin());
  }
  while(down.size()) {
    (*down.begin())->DeregisterFrom(this);
  }
  /*
    for (auto i=up.begin(); i!=up.end(); i++) {
    DeregisterFrom(*i);
    }
  */
}

void DAGnode::DeregisterFrom(DAGnode* parent) {
  if (parent) {
    parent->down.erase(this);
    up.erase(parent);
  }
}

void DAGnode::Register(DAGnode* parent) {
  if (parent) {
    parent->down.insert(this);
    up.insert(parent);
  }
}

void DAGnode::getDot() {
  printf("\t%s%p [label=%s]\n", name.c_str(), (void*)this, name.c_str());
  for (auto i:down)
    printf("\t%s%p -> %s%p\n", name.c_str(), (void*)this, i->GetName().c_str(), (void*)&(*i));
  for (auto i:down)
    i->getDot();
}

void DAGnode::RecursiveRegister(ProbModel* model) {
  auto i=up.begin();
  while ((i!=up.end()) && (*i)->flag)  i++;
  bool up_ok = (i == up.end());

  for (auto i: up) {
    up_ok &= i->flag;
  }
  if (up_ok) {
    model->Register(this);
    flag = true;
    for (auto i : down) {
      i->RecursiveRegister(model);
    }
  }
}

bool DAGnode::CheckUpdateFlags() {
  bool ret = flag;
  if (! flag) {
    cerr << "flag error : " << GetName() << '\n';
  }
  for (auto i: down) {
    ret &= i->CheckUpdateFlags();
  }
  return ret;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * Rnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

double Rnode::Move(double tuning) {
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

//-------------------------------------------------------------------------
// * corrupt
//-------------------------------------------------------------------------
void Rnode::Corrupt(bool bk) {
  value_updated = false;
  localCorrupt(bk);
  for (auto i=down.begin(); i!=down.end(); i++) {
    (*i)->NotifyCorrupt(bk);
  }
}

void Rnode:: NotifyCorrupt(bool bk) {
  localCorrupt(bk);
}

void Rnode::localCorrupt(bool bk) {
  if (bk) {
    bklogprob = logprob;
  }
  flag = false;
}

void Rnode::FullCorrupt(map<DAGnode*,int>& m) {
  localCorrupt(true);
  if (m.find(this) == m.end()) {
    m[this] = 1;
    for (auto i=down.begin(); i!=down.end(); i++) {
      (*i)->FullCorrupt(m);
    }
  }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Rnode::Update() {
  double ret = 0;
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      ret = localUpdate();
      value_updated = true;
      for (auto i=down.begin(); i!=down.end(); i++) {
        ret += (*i)->NotifyUpdate();
      }
    }
  }
  return ret;
}

double Rnode::NotifyUpdate() {
  double ret = 0;
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
      if (GetName() == "BD Chrono") {
      cerr << "BD::NotifyUpdate : " << up_ok << '\n';
      }
    */
    if (up_ok) {
      ret = localUpdate();
      if (! value_updated) {
        value_updated = true;
        for (auto i=down.begin(); i!=down.end(); i++) {
          ret += (*i)->NotifyUpdate();
        }
      }
    }
  }
  return ret;
}

double Rnode::localUpdate() {
  logprob =  logProb();
  flag = true;
  return logprob-bklogprob;
}

double Rnode::FullUpdate(bool check) {
  double ret = 0;
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      ret = localUpdate();
      value_updated = true;
      if ((fabs(ret)>MCMC::MAXDIFF) && check) {
        cerr << "NON ZERO CHECKSUM : " << GetName() << '\n';
        cerr << "number of parents : " << up.size() << '\n';
        throw CheckSumException(ret);
      }
      for (auto i=down.begin(); i!=down.end(); i++) {
        ret += (*i)->FullUpdate(check);
      }
    }
  }
  return ret;

}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Rnode::Restore() {
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      localRestore();
      value_updated = true;
      for (auto i=down.begin(); i!=down.end(); i++) {
        (*i)->NotifyRestore();
      }
    }
  }
}

void Rnode::NotifyRestore() {
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      localRestore();
      if (! value_updated) {
        value_updated = true;
        for (auto i=down.begin(); i!=down.end(); i++) {
          (*i)->NotifyRestore();
        }
      }
    }
  }
}

void Rnode::localRestore() {
  logprob = bklogprob;
  /*
    if (logprob != logProb()) {
    cerr << "error in local restore\n";
    exit(1);
    }
  */
  flag = true;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * Dnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// * corrupt
//-------------------------------------------------------------------------
void Dnode::Corrupt(bool bk) {
  localCorrupt(bk);
  for (auto i=down.begin(); i!=down.end(); i++) {
    (*i)->NotifyCorrupt(bk);
  }
}

void Dnode:: NotifyCorrupt(bool bk) {
  Corrupt(bk);
}

void Dnode::localCorrupt(bool) {
  flag = false;
}

void Dnode::FullCorrupt(map<DAGnode*,int>& m) {
  localCorrupt(true);
  if (m.find(this) == m.end()) {
    m[this] = 1;
    for (auto i=down.begin(); i!=down.end(); i++) {
      (*i)->FullCorrupt(m);
    }
  }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Dnode::Update() {
  double ret = 0;
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      ret = localUpdate();
      for (auto i=down.begin(); i!=down.end(); i++) {
        ret += (*i)->NotifyUpdate();
      }
    }
  }
  return ret;
}

double Dnode::NotifyUpdate() {
  return Update();
}

double Dnode::localUpdate() {
  specialUpdate();
  flag = true;
  return 0;
}

double Dnode::FullUpdate(bool check) {
  double ret = 0;
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      ret = localUpdate();
      for (auto i=down.begin(); i!=down.end(); i++) {
        ret += (*i)->FullUpdate(check);
      }
    }
  }
  return ret;
}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Dnode::Restore() {
  if (! flag) {
    auto i=up.begin();
    while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
    bool up_ok = (i == up.end());
    /*
      bool up_ok = true;
      for (auto i=up.begin(); i!=up.end(); i++) {
      up_ok &= (*i)->isValueUpdated();
      // up_ok &= (*i)->flag;
      }
    */
    if (up_ok) {
      localRestore();
      for (auto i=down.begin(); i!=down.end(); i++) {
        (*i)->NotifyRestore();
      }
    }
  }
}

void Dnode::NotifyRestore() {
  Restore();
}

void Dnode::localRestore() {
  flag = true;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * Mnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// * corrupt
//-------------------------------------------------------------------------
void Mnode::Corrupt(bool bk) {
  flag = false;
  for (auto i=down.begin(); i!=down.end(); i++) {
    (*i)->Corrupt(bk);
  }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Mnode::Update() {
  flag = true;
  double ret = 0;
  for (auto i=down.begin(); i!=down.end(); i++) {
    ret += (*i)->Update();
  }
  return ret;
}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Mnode::Restore() {
  flag = true;
  for (auto i=down.begin(); i!=down.end(); i++) {
    (*i)->Restore();
  }
}
