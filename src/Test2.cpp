#include <cmath>
#include <cstdio>
#include <list>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
#include "utils/Random.hpp"
using namespace std;

double lambda = 4;
bool adaptive = true;

template <class T>
class MyDoubleMove : public MCUpdate {
    Rvar<T>& managedNode;

    // move memory
    vector<double> values;  // list of all values
    double mean;            // on-line mean (approximation)
    int nbVal;              // number of values computed so far (accepted and refused)
    double M2;              // variance * nbVals (approximation)
    int accept;             // number of accepted proposals

  public:
    MyDoubleMove(Rvar<T>& managedNode)
        : managedNode(managedNode), mean(0), nbVal(0), M2(0), accept(0) {}

    double Move(double) override {  // decided to ignore tuning modulator (ie, assume = 1)
        // if node is clamped print a warning message
        if (managedNode.isClamped()) {
            printf("WARNING: Trying to move a clamped node!\n");
            exit(1);
        } else {
            // It seems important to put this BEFORE proposemove. Corrupt sets value_updated to
            // false on the node and its direct children (in Rnode default implementation). The
            // parameter determines if a backup of logprob should be kept.
            managedNode.Corrupt(true);

            // double logHastings = managedNode.ProposeMove(1.0);  // ProposeMove modifies the
            // actual value of the node and returns the log of the Hastings ratio (proposal ratio)

            if (nbVal > 100 and adaptive) {
                (T&)managedNode = mean + Random::sNormal() * lambda * sqrt(M2 / nbVal);
            } else {
                double m = (Random::Uniform() - 0.5);
                managedNode += m;
            }
            while ((managedNode < 0) || (managedNode > 1)) {
                if (managedNode < 0) {
                    (T&)managedNode = -managedNode;
                }
                if (managedNode > 1) {
                    (T&)managedNode = 2 - managedNode;
                }
            }

            double logHastings = 0;

            double logMetropolis = managedNode.Update();

            bool accepted = log(Random::Uniform()) < logMetropolis + logHastings;
            if (!accepted) {
                managedNode.Corrupt(false);
                managedNode.Restore();
            } else {
                accept += 1;
            }
            values.push_back(managedNode);


            nbVal += 1;
            double delta = managedNode - mean;
            mean += delta / nbVal;
            double delta2 = managedNode - mean;
            M2 += delta * delta2;


            // return something
            return (double)accepted;  // for some reason Move seems to return (double)accepted where
                                      // accepted is a
                                      // bool that says if the move was accepted
        }
    }

    void debug() {
        printf("New value %f, mean=%f, variance=%f, acceptance=%f%%\n", double(managedNode), mean,
               M2 / nbVal, (accept * 100.0) / nbVal);
    }
};

class MyModel : public ProbModel {
  public:
    Const<PosReal>* posOne;
    Const<Real>* one;
    Normal* a;
    Exponential* b;
    list<Normal> leaves;
    MyDoubleMove<UnitReal>* mymove;

    MyModel()
        : posOne(new Const<PosReal>(1)),
          one(new Const<Real>(1)),
          a(new Normal(one, posOne)),
          b(new Exponential(posOne, Exponential::MEAN)) {
        for (int i = 0; i < 5; i++) {
            leaves.emplace_back(a, b);
            leaves.back().ClampAt(i < 3 ? 1 : 0);
        }
        RootRegister(one);
        RootRegister(posOne);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override {
        // mymove = new MyDoubleMove<UnitReal>(*p);
        // scheduler.Register(mymove, 1, "p");
        scheduler.Register(new SimpleMove(a, 1.0), 1, "a");
        scheduler.Register(new SimpleMove(b, 1.0), 1, "b");
    }

    double GetLogProb() override { return a->GetLogProb(); }

    void drawSample() override { a->drawSample(); }

    void ToStream(ostream& os) override { os << a->val() << endl; }

    void FromStream(istream&) override {}
};

int main() {
    MyModel model;
    vector<double> results;
    for (int i = 0; i < 100000; i++) {
        model.Move(1.0);
        results.push_back(model.a->val());
    }
    double mean = 0.0;
    for (auto i : results) {
        mean += i;
    }
    double variance = 0.0;
    for (auto i : results) {
        variance += i * i;
    }
    cout << "<DOUBLE_TEST> Mean: " << mean / results.size()
         << " ; variance: " << (variance - (mean * mean / results.size())) / results.size() << endl;
    // model.mymove->debug();
}
