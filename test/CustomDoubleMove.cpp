#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <list>
#include <random>
#include "TestUtils.hpp"
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
#include "utils/Chrono.hpp"
#include "utils/Random.hpp"
using namespace std;
using namespace Eigen;
using namespace TestUtils;


struct normal_random_variable {
    normal_random_variable(MatrixXd const& covar)
        : normal_random_variable(VectorXd::Zero(covar.rows()), covar) {}

    normal_random_variable(VectorXd const& mean, MatrixXd const& covar) : mean(mean) {
        SelfAdjointEigenSolver<MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    VectorXd mean;
    MatrixXd transform;

    VectorXd operator()() const {
        static std::mt19937 gen{std::random_device{}()};
        static std::normal_distribution<> dist;

        return mean +
               transform * VectorXd{mean.size()}.unaryExpr([&](double) { return dist(gen); });
    }
};


void test() {
    int size = 2;
    MatrixXd covar(size, size);
    covar << 1, .5, .5, 1;

    normal_random_variable sample{covar};

    std::cout << sample() << std::endl;
    std::cout << sample() << std::endl;
}


// =======================
//        CONSTANTS
// =======================
// #define REFERENCE_TEST2
double lambda = 0.5;
bool adaptive = false;
// =======================


// =======================
//       CUSTOM MOVE
// =======================
template <class T1, class T2>
class MyDoubleMove : public MCUpdate {
    Rvar<T1>& managedNode1;
    Rvar<T2>& managedNode2;

    Vector2d newValue;

    // move memory
    // vector<Vector2d> values;  // list of all values
    Vector2d mean;  // on-line mean
    Matrix2d covar;
    int t;       // number of values computed so far (accepted and refused)
    int accept;  // number of accepted proposals

  public:
    MyDoubleMove(Rvar<T1>& managedNode1, Rvar<T2>& managedNode2)
        : managedNode1(managedNode1), managedNode2(managedNode2), mean(0, 0), t(0), accept(0) {
        covar << 0, 0, 0, 0;
    }

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
            managedNode2.Corrupt(true);

            // double logHastings = managedNode1.ProposeMove(1.0);  // ProposeMove modifies the
            // actual value of the node and returns the log of the Hastings ratio (proposal ratio)


            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // UPDATE VALUE OF NODES
            if (t < 100 or !adaptive) {
                managedNode1 += Random::sNormal();
                managedNode2 += Random::sNormal();
                newValue = Vector2d(managedNode1, managedNode2);
            } else {
                printf("Adaptive move!\n");
                normal_random_variable sample{lambda * covar};
                newValue = mean + sample();
                (T1&)managedNode1 = newValue(0);
                (T2&)managedNode2 = newValue(1);
            }
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (managedNode2 < 0) {  // posReal specific :/
                (T2&)managedNode2 = -managedNode2;
            }

            double logHastings = 0;
            double logMetropolis = managedNode1.Update() + managedNode2.Update();


            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // DECIDE ACCEPTATION
            bool accepted = log(Random::Uniform()) < logMetropolis + logHastings;
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // ACT DEPENDING OF ACCEPTATION
            if (!accepted) {
                managedNode1.Corrupt(false);
                managedNode1.Restore();
                managedNode2.Corrupt(false);
                managedNode2.Restore();
            } else {
                accept += 1;
            }
            // values.push_back(Vector2d(managedNode1, managedNode2));
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // UPDATING THINGS IN AN ON-LINE FASHION
            t += 1;  // counting iterations (t=k)
            newValue = Vector2d(managedNode1, managedNode2);

            Matrix2d firstOuterProduct =
                mean * mean.transpose();  // is useful for the variance update below

            Vector2d tmp = ((t - 1.0) / t) * mean + (1.0 / t) * newValue;  // updating mean
            mean = tmp;

            Matrix2d tmp2 =
                ((t - 1.0) / t) * covar +
                (1.0 / t) * (t * firstOuterProduct - (t + 1.0) * (mean * mean.transpose()) +
                             newValue * newValue.transpose());
            covar = tmp2;
            // cout << "COVAR = \n" << covar << endl;
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            return (double)accepted;  // for some reason Move seems to return (double)accepted where
                                      // accepted is a
                                      // bool that says if the move was accepted
        }
    }

    void debug() {
        printf("New value %f/%f, first=%f, second=%f, acceptance=%f%%\n", double(managedNode1),
               double(managedNode2), mean(0), mean(1), (accept * 100.0) / t);
        cout << covar << endl;
        cassert(mean(0), 2.27, 0.1);
        cassert(mean(1), 2.18, 0.1);
        cassert(covar(0, 0), 0.32, 0.1);
        cassert(covar(1, 1), 1.0, 0.1);
    }
};


// ======================
//         MODEL
// ======================
class MyModel : public ProbModel {
  public:
    // graphical model
    Const<PosReal>* posOne;
    Const<Real>* one;
    Normal* a;
    Exponential* b;
    list<Normal> leaves;

    // moves
    MyDoubleMove<Real, PosReal>* myMove;

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
        myMove = new MyDoubleMove<Real, PosReal>(*a, *b);
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


// ======================
//      AUX FUNCTIONS
// ======================
void printCaracs(vector<double> data, string name) {
    double mean = 0.0;
    for (auto i : data) {
        mean += i;
    }
    double variance = 0.0;
    for (auto i : data) {
        variance += i * i;
    }
    double finalMean = mean / data.size();
    double finalVariance = (variance - (mean * mean / data.size())) / data.size();
    cout << "<" << name << "> Mean: " << finalMean << " ; variance: " << finalVariance << endl;
    double expectedMean = (name == "a") ? 2.27 : 2.18;
    double expectedVariance = name == "a" ? 0.32 : 1.0;
    cassert(expectedMean, finalMean, 0.1);
    cassert(expectedVariance, finalVariance, 0.1);
}


// ======================
//           MAIN
// ======================
int main() {
    MyModel model;
    vector<double> resultsA, resultsB;
    MeasureTime myTimer;
    for (int i = 0; i < 1000000; i++) {
        model.Move(1.0);
        resultsA.push_back(model.a->val());
        resultsB.push_back(model.b->val());
    }
    myTimer.print();
    printCaracs(resultsA, "a");  // expected 2.27
    printCaracs(resultsB, "b");  // expected 2.18
#ifndef REFERENCE_TEST2
    model.myMove->debug();
#endif
}
