#include <cmath>
#include <cstdio>
#include <list>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
#include "utils/Random.hpp"
using namespace std;


// =======================
//        CONSTANTS
// =======================
// #define REFERENCE_TEST2
double lambda = 4;
bool adaptive = false;
// =======================



template <class T1, class T2>
class MyDoubleMove : public MCUpdate {
    Rvar<T1>& managedNode1;
    Rvar<T2>& managedNode2;

    // move memory
    vector<double> values;  // list of all values
    double mean;            // on-line mean (approximation)
    int nbVal;              // number of values computed so far (accepted and refused)
    double M2;              // variance * nbVals (approximation)
    int accept;             // number of accepted proposals

  public:
    MyDoubleMove(Rvar<T1>& managedNode1, Rvar<T2>& managedNode2)
        : managedNode1(managedNode1),
          managedNode2(managedNode2),
          mean(0),
          nbVal(0),
          M2(0),
          accept(0) {}

    double Move(double) override {  // decided to ignore tuning modulator (ie, assume = 1)
        // if node is clamped print a warning message
        if (managedNode1.isClamped()) {
            printf("WARNING: Trying to move a clamped node!\n");
            exit(1);
        } else {
            // It seems important to put this BEFORE proposemove. Corrupt sets value_updated to
            // false on the node and its direct children (in Rnode default implementation). The
            // parameter determines if a backup of logprob should be kept.
            managedNode1.Corrupt(true);

            // double logHastings = managedNode1.ProposeMove(1.0);  // ProposeMove modifies the
            // actual value of the node and returns the log of the Hastings ratio (proposal ratio)

            if (nbVal > 100 and adaptive) {
                (T1&)managedNode1 = mean + Random::sNormal() * lambda * sqrt(M2 / nbVal);
            } else {
                double m = (Random::Uniform() - 0.5);
                managedNode1 += m;
            }

            double logHastings = 0;

            double logMetropolis = managedNode1.Update();

            bool accepted = log(Random::Uniform()) < logMetropolis + logHastings;
            if (!accepted) {
                managedNode1.Corrupt(false);
                managedNode1.Restore();
            } else {
                accept += 1;
            }
            values.push_back(managedNode1);


            nbVal += 1;
            double delta = managedNode1 - mean;
            mean += delta / nbVal;
            double delta2 = managedNode1 - mean;
            M2 += delta * delta2;


            // return something
            return (double)accepted;  // for some reason Move seems to return (double)accepted where
                                      // accepted is a
                                      // bool that says if the move was accepted
        }
    }

    void debug() {
        printf("New value %f, mean=%f, variance=%f, acceptance=%f%%\n", double(managedNode1), mean,
               M2 / nbVal, (accept * 100.0) / nbVal);
    }
};

class MyModel : public ProbModel {
  public:
    // graphical model
    Const<PosReal>* posOne;
    Const<Real>* one;
    Normal* a;
    Exponential* b;
    list<Normal> leaves;

    // moves
    MyDoubleMove<PosReal, Real>* myMove;

    MyModel()
        : posOne(new Const<PosReal>(1)),
          one(new Const<Real>(1)),
          a(new Normal(one, posOne)),
          b(new Exponential(posOne, Exponential::MEAN)) {
        for (int i = 0; i < 5; i++) {
            leaves.emplace_back(a, b);
            leaves.back().ClampAt(i < 3 ? 4 : 1);
        }
        RootRegister(one);
        RootRegister(posOne);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override {
#ifndef REFERENCE_TEST2
        myMove = new MyDoubleMove<PosReal, Real>(*b, *a);
        scheduler.Register(myMove, 1, "ab");
#else
        // scheduler.Register(mymove, 1, "p");
        scheduler.Register(new SimpleMove(a, 1.0), 1, "a");
        scheduler.Register(new SimpleMove(b, 1.0), 1, "b");
#endif
    }

    double GetLogProb() override { return a->GetLogProb() + b->GetLogProb(); }

    void drawSample() override {
        a->drawSample();
        b->drawSample();
    }

    void ToStream(ostream& os) override { os << a->val() << endl; }

    void FromStream(istream&) override {}
};


void printCaracs(vector<double> data, string name) {
    double mean = 0.0;
    for (auto i : data) {
        mean += i;
    }
    double variance = 0.0;
    for (auto i : data) {
        variance += i * i;
    }
    cout << "<" << name << "> Mean: " << mean / data.size()
         << " ; variance: " << (variance - (mean * mean / data.size())) / data.size() << endl;
}


int main() {
    MyModel model;
    vector<double> resultsA, resultsB;
    for (int i = 0; i < 1000000; i++) {
        model.Move(1.0);
        resultsA.push_back(model.a->val());
        resultsB.push_back(model.b->val());
    }
    printCaracs(resultsA, "a");  // expected 2.27
    printCaracs(resultsB, "b");  // expected 2.18

    // model.mymove->debug();
}
