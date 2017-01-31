#include "ProbModel.hpp"
#include "RandomTypes.hpp"

class PoissonGammaModel : public ProbModel	{

  int N;
  int* data;

  Const<PosReal>* One;

  Exponential* theta;
  Exponential* sigma;


  Gamma** omega;
  Product** rate;
  Poisson** X;

public:

  PoissonGammaModel(int inN, int* indata)	{

    N = inN;
    data = indata;
    One = new Const<PosReal>(1);

    sigma = new Exponential(One,Exponential::MEAN);
    theta = new Exponential(One,Exponential::MEAN);

    omega = new Gamma*[N];
    rate = new Product*[N];
    X = new Poisson*[N];

    for (int i=0; i<N; i++)	{
      omega[i] = new Gamma(theta,theta);
      rate[i] = new Product(omega[i],sigma);
      X[i] = new Poisson(rate[i]);
      X[i]->ClampAt(data[i]);
    }

    RootRegister(One);
    Register();
    Update();
    MakeScheduler();
  }


  ~PoissonGammaModel() {}


  virtual void MakeScheduler()	{

    scheduler.Register(new SimpleMove(sigma,0.1),1,"sigma");
    scheduler.Register(new SimpleMove(sigma,0.01),1,"sigma");

    scheduler.Register(new SimpleMove(theta,0.1),1,"theta");
    scheduler.Register(new SimpleMove(theta,0.01),1,"theta");

    for (int i=0; i<N; i++)	{
      scheduler.Register(new SimpleMove(omega[i],1),1,"omega");
      scheduler.Register(new SimpleMove(omega[i],0.1),1,"omega");
    }
  }

  double GetMeanOmega()	{

    double mean = 0;
    for (int i=0; i<N; i++)	{
      mean += omega[i]->val();
    }

    mean /= N;
    return mean;
  }

  void TraceHeader(ostream& os)	{
    os << "logprob\ttheta\tsigma\tmeeanomega\n";
  }

  void Trace(ostream& os)	{
    os << GetLogProb() << '\t' << theta->val() << '\t' << sigma->val() << '\t' << GetMeanOmega() << '\n';
  }

  double GetLogProb()	{

    double tot = 0;
    tot += sigma->GetLogProb();
    tot += theta->GetLogProb();
    return tot;
  }

  void drawSample()	{

    sigma->Sample();
    theta->Sample();
    for (int i=0; i<N; i++)	{
      omega[i]->Sample();
    }
  }

  void ToStream(ostream&) {}
  void FromStream(istream&) {}

};


int main(int, char* argv[])	{

  ifstream is(argv[1]);
  string name = argv[2];

  int N;
  is >> N;
  int* data = new int[N];
  for (int i=0; i<N; i++)	{
    is >> data[i];
  }

  PoissonGammaModel* model = new PoissonGammaModel(N,data);
  ofstream os((name + ".trace").c_str());
  model->TraceHeader(os);

  while (1)	{
    model->Move(1.0);
    model->Trace(os);
  }
}
