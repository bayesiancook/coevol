#include <cstdio>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
using namespace std;

class MySimpleMove : public MCUpdate {};

class MyModel : public ProbModel {
    Const<PosReal>* One;
    Const<PosReal>* Two;
    Exponential* A;
    Exponential* B;
    Product* C;

  public:
    MyModel() {
        One = new Const<PosReal>(1);
        Two = new Const<PosReal>(2);
        A = new Exponential(One, Exponential::MEAN);
        B = new Exponential(Two, Exponential::MEAN);
        C = new Product(A, B);

        RootRegister(One);
        RootRegister(Two);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    void MakeScheduler() override {
        scheduler.Register(new SimpleMove(A, 0.5), 1, "A");
        scheduler.Register(new SimpleMove(B, 0.5), 1, "B");
    }

    double GetLogProb() override { return A->GetLogProb() + B->GetLogProb(); }

    void drawSample() override {
        A->drawSample();
        B->drawSample();
    }

    void ToStream(ostream&) override {
        cout << GetLogProb() << '\t' << A->val() << '\t' << B->val() << '\t' << C->val() << endl;
    }

    void FromStream(istream&) override {}
};

int main() {
    MyModel model;
    for (int i = 0; i < 100; i++) {
        model.Move(1.0);
        model.ToStream(cout);
    }
}
