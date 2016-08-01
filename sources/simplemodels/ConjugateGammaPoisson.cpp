
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugatePoisson.h"


class binomial : public Rvar<Int>	{

	static const double InfProb = -1050;
	public:

	binomial(Rvar<Int>* inN, Var<PosReal>* inp)	{
		N = inN;
		p = inp;
		Register(N);
		Register(p);
	}
		
	virtual ~binomial() {}


	void 		drawSample()	{

		int val = 0;
		for (int k=0; k<N->val(); k++)	{
			if (Random::Uniform() < p->val())	{
				val++;
			}
		}
		cerr << val << '\n';
		setval(val);
	}

	private:

	virtual double 		logProb()	{
		if (((this->val()) < 0) || ((this->val()) > N->val()))	{
			return InfProb;
		}
		double ret = *this * log(p->val()) + (N->val() - *this) * log(1 - p->val());
		// ret += Random::logGamma(N->val()) - Random::logGamma(*this) - Random::logGamma(N->val() - this->val());
		for (int i=1; i<=*this; i++)	{
			ret -= log((double) i);
			ret += log((double) (N->val()+1-i));
		}
		return ret;
	}

	Var<Int>*	N;
	Var<PosReal>* 	p;
};
 
class PoissonBinModel : public ProbModel	{

	public:

	Dvar<PosReal>* mu;
	Dvar<PosReal>* alpha;
	ConjugateGamma* r;
	ConjugatePoisson* n;
	binomial* k;
	Dvar<PosReal>* p;

	PoissonBinModel(double inmu, double inp, int ink)	{

		mu = new Const<PosReal>(inmu);
		alpha = new Const<PosReal>(1);
		p = new Const<PosReal>(inp);
		r = new ConjugateGamma(alpha,mu);
		n = new ConjugatePoisson(r);
		k = new binomial(n,p);
		k->setval(ink);
		k->Clamp();

		/*
		n->setval(9);
		n->Clamp();
		*/

		RootRegister(mu);
		RootRegister(alpha);
		RootRegister(p);
		Register();

		Sample();
		Update();
		Trace(cerr);

		// MakeScheduler();
	}

	/*
	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(r,5),1,"r");
		scheduler.Register(new SimpleMove(r,1),1,"r");
		scheduler.Register(new SimpleMove(n,8),1,"n");
		scheduler.Register(new SimpleMove(n,2),1,"n");
	}
	*/

	void MakeScheduler() {}

	double Move(double tuning = 1)	{
		/*
		r->Move(5);
		r->Move(1);
		*/
		r->Integrate();
		n->Move(8);
		n->Move(2);
		r->Resample();
		cerr << r->val() << '\n';
		return 1;
	}

	double GetLogProb()	{
		return r->GetLogProb() + n->GetLogProb() + k->GetLogProb();
	}

	void drawSample()	{
		r->Sample();
		n->Sample();
	}

	void ToStream(ostream& os)	{
		os << *r << '\t' << *n << '\n';
	}

	void FromStream(istream& is)	{
		is >> *r >> *n;
	}

	double GetRate()	{
		return r->val();
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *r << '\t' << *n << '\n';
	}
};

int main()	{

	double mu = 0.05;
	double p = 0.1;
	int k = 10;
	double alpha = 1 + k;
	double beta = mu + p;
	double mean = alpha / beta;
	double var = alpha / beta / beta;

	PoissonBinModel* model = new PoissonBinModel(mu,p,k);

	int burnin = 100;
	int nrep = 1000;
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
		double tmp = model->GetRate();
		empmean += tmp;
		empvar += tmp * tmp;
	}
	empmean /= nrep;
	empvar /= nrep;
	empvar -= empmean * empmean;

	ofstream mos("monitor");
	model->Monitor(mos);

	cout << "true values : " << mean << '\t' << var << '\n';
	cout << "estimated   : " << empmean << '\t' << empvar << '\n';
}

