#include <cstdio>
#include <list>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
using namespace std;

class MySimpleMove : public MCUpdate {};

class MyModel : public ProbModel {
  public:
    Const<PosReal>* One;
    Beta* p;
    list<Binomial> leaves;

    MyModel() : One(new Const<PosReal>(1)), p(new Beta(One, One)) {
        for (int i=0; i<5; i++) {
            leaves.emplace_back(1, p);
            leaves.back().ClampAt(i<3?1:0);
        }
        RootRegister(One);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override { scheduler.Register(new SimpleMove(p, 1.0), 1, "p"); }

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
