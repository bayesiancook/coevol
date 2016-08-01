
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugatePoisson.h"

 
class Jeffreys : public Rvar<PosReal>	{

	public:

	Jeffreys() {}

	double logProb() {return -log(*this);}

	void drawSample() {value = 1;}

};

class PoissonBinModel : public ProbModel	{

	public:

	Jeffreys* mu;
	Dvar<PosReal>* alpha;
	ConjugateGamma* r;
	ConjugatePoisson* n;

	PoissonBinModel(int inn)	{

		mu = new Jeffreys;
		alpha = new Const<PosReal>(1);
		r = new ConjugateGamma(alpha,mu);
		n = new ConjugatePoisson(r);

		n->setval(inn);
		n->Clamp();

		RootRegister(mu);
		RootRegister(alpha);
		Register();

		Sample();
		Update();
		Trace(cerr);

	}

	void MakeScheduler() {}

	double Move(double tuning = 1)	{
		/*
		mu->Move(5);
		mu->Move(1);
		r->Move(5);
		r->Move(1);
		*/
		r->Integrate();
		mu->Move(5);
		mu->Move(1);
		r->Resample();
		return 1;
	}

	double GetLogProb()	{
		return mu->GetLogProb() + r->GetLogProb();
	}

	void drawSample()	{
		mu->Sample();
		r->Sample();
	}

	void ToStream(ostream& os)	{
		os << *mu << '\t' << *r << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu >> *r;
	}

	double GetMu()	{
		return mu->val();
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *mu << '\t' << *r << '\t' << *n << '\n';
	}
};

int main()	{

		
	int n = 9;
	PoissonBinModel* model = new PoissonBinModel(n);

	int burnin = 100;
	int nrep = 10000;
	ofstream os("out");

	for (int i=0; i<burnin; i++)	{
		for (int k=0; k<10; k++)	{
			model->Move(1);
		}
		model->Trace(os);
	}

	double empmean = 0;
	double empvar = 0;
	for (int i=0; i<nrep; i++)	{
		for (int k=0; k<10; k++)	{
			model->Move(1);
		}
		model->Trace(os);
		double tmp = model->GetMu();
		empmean += tmp;
		empvar += tmp * tmp;
	}
	empmean /= nrep;
	empvar /= nrep;
	empvar -= empmean * empmean;

	ofstream mos("monitor");
	model->Monitor(mos);

	cout << "true value : " << 1.0 / 8 << '\n';
	cout << "estimated  : " << empmean << '\t' << empvar << '\n';
}

