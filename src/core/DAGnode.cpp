#include "core/DAGnode.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "core/ProbModel.hpp"
#include "utils/Exception.hpp"
#include "utils/Random.hpp"
using namespace std;


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * DAGnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
bool DAGnode::initmode = true;  // (VL) FIXME

const double MCMC::MAXDIFF = 1e-4;

DAGnode::~DAGnode() { Detach(); }

int DAGnode::GetChildrenNumber() {
    int tot = 0;
    for (auto i = down.begin(); i != down.end(); i++) {
        tot++;
    }
    return tot;
}

void DAGnode::Detach() {
    while (!up.empty() != 0u) {
        DeregisterFrom(*up.begin());
    }
    while (!down.empty() != 0u) {
        (*down.begin())->DeregisterFrom(this);
    }
    /*
      for (auto i=up.begin(); i!=up.end(); i++) {
      DeregisterFrom(*i);
      }
    */
}

void DAGnode::DeregisterFrom(DAGnode* parent) {
    if (parent != nullptr) {
        parent->down.erase(this);
        up.erase(parent);
    }
}

void DAGnode::Register(DAGnode* parent) {
    if (parent != nullptr) {
        parent->down.insert(this);
        up.insert(parent);
    }
}

set<string> DAGnode::getDotNodes() {
    ostringstream stringStream;
    set<string> result;
    stringStream << "\tNode" << this << " [label=\"" << name << "\"]" << endl;
    ;
    result.insert(stringStream.str());
    for (auto i : down) {
        set<string> tmp = i->getDotNodes();
        result.insert(tmp.begin(), tmp.end());
    }
    return result;
}

set<string> DAGnode::getDotVertices() {
    ostringstream stringStream;
    set<string> result;
    for (auto i : down) {
        stringStream << "\tNode" << this << " -> Node" << &(*i) << endl;
        ;
        result.insert(stringStream.str());
        stringStream.str("");
    }
    for (auto i : down) {
        set<string> tmp = i->getDotVertices();
        result.insert(tmp.begin(), tmp.end());
    }
    return result;
}

void DAGnode::RecursiveRegister(ProbModel* model) {
    auto i = up.begin();
    while ((i != up.end()) && (*i)->updateFlag) {
        i++;
    }
    bool up_ok = (i == up.end());

    for (auto i : up) {
        up_ok &= static_cast<int>(i->updateFlag);
    }
    if (up_ok) {
        model->Register(this);
        updateFlag = true;
        for (auto i : down) {
            i->RecursiveRegister(model);
        }
    }
}

bool DAGnode::CheckUpdateFlags() {
    bool ret = updateFlag;
    if (!updateFlag) {
        cerr << "updateFlag error : " << GetName() << '\n';
    }
    for (auto i : down) {
        ret &= static_cast<int>(i->CheckUpdateFlags());
    }
    return ret;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * Rnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
set<string> Rnode::getDotNodes() {  // FIXME refactor with DAGnode version
    ostringstream stringStream;
    set<string> result;
    stringStream << "\tNode" << dynamic_cast<DAGnode*>(this) << " [label=\"" << name
                 << "\", shape=rectangle]" << endl;
    ;
    result.insert(stringStream.str());
    for (auto i : down) {
        set<string> tmp = i->getDotNodes();
        result.insert(tmp.begin(), tmp.end());
    }
    return result;
}


double Rnode::Move(double tuning) {
    if (!isClamped()) {
        Corrupt(true);
        double logHastings = ProposeMove(tuning);
        double deltaLogProb = Update();
        double logRatio = deltaLogProb + logHastings;
        bool accepted = (log(Random::Uniform()) < logRatio);
        if (!accepted) {
            Corrupt(false);
            Restore();
        }
        return (double)accepted;
    }
    return 1;
}

//-------------------------------------------------------------------------
// * corrupt
//-------------------------------------------------------------------------
void Rnode::Corrupt(bool bk) {
    value_updated = false;
    localCorrupt(bk);
    for (auto i : down) {
        i->NotifyCorrupt(bk);
    }
}

void Rnode::NotifyCorrupt(bool bk) { localCorrupt(bk); }

void Rnode::localCorrupt(bool bk) {
    if (bk) {
        bklogprob = logprob;
    }
    updateFlag = false;
}

void Rnode::FullCorrupt(map<DAGnode*, int>& m) {
    localCorrupt(true);
    if (m.find(this) == m.end()) {
        m[this] = 1;
        for (auto i : down) {
            i->FullCorrupt(m);
        }
    }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Rnode::Update() {
    double ret = 0;
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        if (i == up.end()) {
            ret = localUpdate();
            value_updated = true;
            for (auto i : down) {
                ret += i->NotifyUpdate();
            }
        }
    }
    return ret;
}

double Rnode::NotifyUpdate() {
    double ret = 0;
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
          if (GetName() == "BD Chrono") {
          cerr << "BD::NotifyUpdate : " << up_ok << '\n';
          }
        */
        if (up_ok) {
            ret = localUpdate();
            if (!value_updated) {
                value_updated = true;
                for (auto i : down) {
                    ret += i->NotifyUpdate();
                }
            }
        }
    }
    return ret;
}

double Rnode::localUpdate() {
    logprob = logProb();
    updateFlag = true;
    return logprob - bklogprob;
}

double Rnode::FullUpdate(bool check) {
    double ret = 0;
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            ret = localUpdate();
            value_updated = true;
            if ((fabs(ret) > MCMC::MAXDIFF) && check) {
                cerr << "NON ZERO CHECKSUM : " << GetName() << '\n';
                cerr << "number of parents : " << up.size() << '\n';
                throw CheckSumException(ret);
            }
            for (auto i : down) {
                ret += i->FullUpdate(check);
            }
        }
    }
    return ret;
}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Rnode::Restore() {
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            localRestore();
            value_updated = true;
            for (auto i : down) {
                i->NotifyRestore();
            }
        }
    }
}

void Rnode::NotifyRestore() {
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            localRestore();
            if (!value_updated) {
                value_updated = true;
                for (auto i : down) {
                    i->NotifyRestore();
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
    updateFlag = true;
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
    for (auto i : down) {
        i->NotifyCorrupt(bk);
    }
}

void Dnode::NotifyCorrupt(bool bk) { Corrupt(bk); }

void Dnode::localCorrupt(bool /*unused*/) { updateFlag = false; }

void Dnode::FullCorrupt(map<DAGnode*, int>& m) {
    localCorrupt(true);
    if (m.find(this) == m.end()) {
        m[this] = 1;
        for (auto i : down) {
            i->FullCorrupt(m);
        }
    }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Dnode::Update() {
    double ret = 0;
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            ret = localUpdate();
            for (auto i : down) {
                ret += i->NotifyUpdate();
            }
        }
    }
    return ret;
}

double Dnode::NotifyUpdate() { return Update(); }

double Dnode::localUpdate() {
    specialUpdate();
    updateFlag = true;
    return 0;
}

double Dnode::FullUpdate(bool check) {
    double ret = 0;
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            ret = localUpdate();
            for (auto i : down) {
                ret += i->FullUpdate(check);
            }
        }
    }
    return ret;
}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Dnode::Restore() {
    if (!updateFlag) {
        auto i = up.begin();
        while ((i != up.end()) && ((*i)->isValueUpdated())) {
            i++;
        }
        bool up_ok = (i == up.end());
        /*
          bool up_ok = true;
          for (auto i=up.begin(); i!=up.end(); i++) {
          up_ok &= (*i)->isValueUpdated();
          // up_ok &= (*i)->updateFlag;
          }
        */
        if (up_ok) {
            localRestore();
            for (auto i : down) {
                i->NotifyRestore();
            }
        }
    }
}

void Dnode::NotifyRestore() { Restore(); }

void Dnode::localRestore() { updateFlag = true; }


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// * Mnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// * corrupt
//-------------------------------------------------------------------------
void Mnode::Corrupt(bool bk) {
    updateFlag = false;
    for (auto i : down) {
        i->Corrupt(bk);
    }
}

//-------------------------------------------------------------------------
// * update
//-------------------------------------------------------------------------
double Mnode::Update() {
    updateFlag = true;
    double ret = 0;
    for (auto i : down) {
        ret += i->Update();
    }
    return ret;
}

//-------------------------------------------------------------------------
// * restore
//-------------------------------------------------------------------------
void Mnode::Restore() {
    updateFlag = true;
    for (auto i : down) {
        i->Restore();
    }
}
