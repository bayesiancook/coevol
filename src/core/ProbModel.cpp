#include <fstream>
#include <iostream>
#include <map>
#include <set>
using namespace std;

#include "core/ProbModel.hpp"


ProbModel::ProbModel() : scheduler(this) {}

ProbModel::~ProbModel() = default;

bool ProbModel::CheckUpdateFlags() {
    bool ret = true;
    for (auto i : root) {
        ret &= static_cast<int>(i->CheckUpdateFlags());
    }
    return ret;
}

void ProbModel::Register(DAGnode* var) { state.insert(var); }

void ProbModel::RootRegister(DAGnode* var) { root.insert(var); }

void ProbModel::Corrupt() {
    map<DAGnode*, int> m;
    for (auto i : root) {
        i->FullCorrupt(m);
    }
}

double ProbModel::Update(bool check) {
    Corrupt();
    double total = 0;
    for (auto i : root) {
        total += i->FullUpdate(check);
    }
    DAGnode::initmode = false;
    return total;
}

void ProbModel::Register() {
    if (!state.empty()) {
        cerr << "error : state is not empty\n";
        exit(1);
    }
    Corrupt();
    for (auto i : root) {
        i->RecursiveRegister(this);
    }
    cerr << "model size : " << state.size() << '\n';
}

void ProbModel::getDot() {
    ofstream myfile("tmp.dot");
    myfile << "digraph {" << endl;
    for (auto i : root) {
        for (auto j : i->getDotNodes()) {
            myfile << j;
        }
    }
    for (auto i : root) {
        for (auto j : i->getDotVertices()) {
            myfile << j;
        }
    }
    myfile << "}" << endl;
    myfile.close();
}

double ProbModel::Move(double tuning_modulator, int ncycle, bool verbose, bool check) {
    scheduler.Cycle(tuning_modulator, ncycle, verbose, check);
    return 1;
}

void ProbModel::Monitor(ostream& os, ostream& osdetail) {
    scheduler.ToStream(os, osdetail);
    Details(osdetail);
}
