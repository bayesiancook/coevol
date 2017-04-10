#include <fstream>
#include "core/ProbModel.hpp"
#include "core/RandomTypes.hpp"

using namespace std;

class PoissonGammaModel : public ProbModel {
    int N;
    vector<int> data;

    Const<PosReal>* One;

    Exponential* theta;
    Exponential* sigma;

    Gamma** omega;
    Product** rate;
    Poisson** X;

  public:
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

int main(int /*unused*/, char* /*unused*/ []) {
    cout << "Youpi" << std::endl;
    string name{"youpi"};

    vector<int> data{1,2,4,3,2};

    auto model = new PoissonGammaModel(5, data);
    cout << "youpla\n";
    ofstream os((name + ".trace").c_str());
    cout << "boum\n";
    model->TraceHeader(os);
    cout << "poin\n";

    int i = 0;
    while (i < 100000) {
        cout << "la\n";
        model->Move(1.0);
        model->Trace(os);
        i++;
    }
    cout << "tadaaaa\n";
}
