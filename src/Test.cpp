#include <cmath>
#include <cstdio>
#include <list>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
#include "utils/Random.hpp"
using namespace std;

template <class T>
class MySimpleMove : public MCUpdate {
    Rvar<T>* managedNode;

  public:
    MySimpleMove(Rvar<T>* managedNode) : managedNode(managedNode) {}

    double Move(double) override {  // decided to ignore tuning modulator (ie, assume = 1)
        // if node is clamped print a warning message
        if (managedNode->isClamped()) {
            printf("WARNING: Trying to move a clamped node!\n");
            exit(1);
        } else {
            // It seems important to put this BEFORE proposemove. Corrupt sets value_updated to
            // false on the node and its direct children (in Rnode default implementation). The
            // parameter determines if a backup of logprob should be kept.
            managedNode->Corrupt(true);

            // double logHastings = managedNode->ProposeMove(1.0);  // ProposeMove modifies the actual
                                                                 // value of the node and returns
                                                                 // the log of the Hastings ratio
                                                                 // (proposal ratio)

            double m = (Random::Uniform() - 0.5);
            managedNode->value += m;
            while ((managedNode->value < 0) || (managedNode->value > 1)) {
                if (managedNode->value < 0) {
                    managedNode->value = -managedNode->value;
                }
                if (managedNode->value > 1) {
                    managedNode->value = 2 - managedNode->value;
                }
            }

            double logHastings = 0;

            double logMetropolis = managedNode->Update();

            bool accepted = log(Random::Uniform()) < logMetropolis + logHastings;
            if (!accepted) {
                managedNode->Corrupt(false);
                managedNode->Restore();
            }

            // return somehting
            return (double)accepted;  // for some reason Move seems to return (double)accepted where
                                      // accepted is a
                                      // bool that says if the move was accepted
        }
    }
};

class MyModel : public ProbModel {
  public:
    Const<PosReal>* one;
    Beta* p;
    list<Binomial> leaves;

    MyModel() : one(new Const<PosReal>(1)), p(new Beta(one, one)) {
        for (int i = 0; i < 5; i++) {
            leaves.emplace_back(1, p);
            leaves.back().ClampAt(i < 3 ? 1 : 0);
        }
        RootRegister(one);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override { scheduler.Register(new MySimpleMove<UnitReal>(p), 1, "p"); }

    double GetLogProb() override { return p->GetLogProb(); }

    void drawSample() override { p->drawSample(); }

    void ToStream(ostream& os) override { os << p->val() << endl; }

    void FromStream(istream&) override {}
};

int main() {
    MyModel model;
    vector<double> results;
    for (int i = 0; i < 100000; i++) {
        model.Move(1.0);
        results.push_back(model.p->val());
    }
    double mean = 0.0;
    for (auto i : results) {
        mean += i;
    }
    cout << mean / results.size() << endl;
}
