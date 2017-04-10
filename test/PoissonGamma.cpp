#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <fstream>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"
#include "doctest.h"
using namespace std;
using namespace doctest;


class PoissonGammaModel : public ProbModel {
public:
    int N;
    vector<int> data;

    Const<PosReal>* One;

    Exponential* theta;
    Exponential* sigma;

    Gamma** omega;
    Product** rate;
    Poisson** X;

    PoissonGammaModel(int inN, vector<int> indata) {
        N = inN;
        data = indata;
        One = new Const<PosReal>(1);

        sigma = new Exponential(One, Exponential::MEAN);
        theta = new Exponential(One, Exponential::MEAN);

        omega = new Gamma*[N];
        rate = new Product*[N];
        X = new Poisson*[N];

        for (int i = 0; i < N; i++) {
            omega[i] = new Gamma(theta, theta);
            rate[i] = new Product(omega[i], sigma);
            X[i] = new Poisson(rate[i]);
            X[i]->ClampAt(data[i]);
        }

        RootRegister(One);
        Register();
        Update();
        MakeScheduler();
        getDot();
    }

    ~PoissonGammaModel() override = default;

    void MakeScheduler() override {
        scheduler.Register(new SimpleMove(sigma, 0.1), 1, "sigma");
        scheduler.Register(new SimpleMove(sigma, 0.01), 1, "sigma");

        scheduler.Register(new SimpleMove(theta, 0.1), 1, "theta");
        scheduler.Register(new SimpleMove(theta, 0.01), 1, "theta");

        cout << "SIZE = " << N << endl;

        for (int i = 0; i < N; i++) {
            scheduler.Register(new SimpleMove(omega[i], 1), 1, "omega");
            scheduler.Register(new SimpleMove(omega[i], 0.1), 1, "omega");
        }
    }

    double GetMeanOmega() {
        double mean = 0;
        for (int i = 0; i < N; i++) {
            mean += omega[i]->val();
        }

        mean /= N;
        return mean;
    }

    void TraceHeader(ostream& os) override { os << "logprob\ttheta\tsigma\tmeeanomega\n"; }

    void Trace(ostream& os) override {
        os << GetLogProb() << '\t' << theta->val() << '\t' << sigma->val() << '\t' << GetMeanOmega()
           << '\n';
    }

    double GetLogProb() override {
        double tot = 0;
        tot += sigma->GetLogProb();
        tot += theta->GetLogProb();
        return tot;
    }

    void drawSample() override {
        sigma->Sample();
        theta->Sample();
        for (int i = 0; i < N; i++) {
            omega[i]->Sample();
        }
    }

    void ToStream(ostream& /*os*/) override {}
    void FromStream(istream& /*is*/) override {}
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
    double expectedMean = 2.1;
    double expectedVariance = 1.45;

    CHECK(finalMean == Approx(expectedMean).epsilon(0.05));
    CHECK(finalVariance == Approx(expectedVariance).epsilon(0.15));
}


TEST_CASE("Testing PoissonGamma model against fixed values.") {
    string name{"youpi"};

    vector<int> indata{1, 2, 4, 3, 2};
    vector<double> data{};

    auto model = new PoissonGammaModel(5, indata);
    ofstream os((name + ".trace").c_str());
    model->TraceHeader(os);

    int i = 0;
    while (i < 250000) {
        model->Move(1.0);
        model->Trace(os);
        data.push_back(model->theta->val());
        i++;
    }

    printCaracs(data, "theta");
}
