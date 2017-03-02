#include <cstdio>
#include <vector>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
using namespace std;

class MySimpleMove : public MCUpdate {};

class MyModel : public ProbModel {
    Const<PosReal>* One;
    Beta* p;
    vector<Binomial> vec;

  public:
    MyModel() : One(new Const<PosReal>(1)), p(new Beta(One, One)), vec(5, Binomial(1, p)) {
        for (unsigned int i = 0; i < vec.size(); i++) {
            vec.at(i).ClampAt(i < 3 ? 1 : 0);
        }

        RootRegister(One);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override {
        scheduler.Register(new SimpleMove(p, 1.0), 1, "p");  // TODO what are those parameters
    }

    double GetLogProb() override { return p->GetLogProb(); }

    void drawSample() override { p->drawSample(); }

    void ToStream(ostream& os) override { os << p->val() << endl; }

    void FromStream(istream&) override {}
};

int main() {
    MyModel model;
    for (int i = 0; i < 100; i++) {
        model.Move(1.0);
        model.ToStream(cout);
    }
}
