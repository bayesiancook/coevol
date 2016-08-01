
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugatePoisson.h"


class uniform : public Rvar<UnitReal>	{

	public:

	uniform() {}

	void drawSample()	{
		setval(Random::Uniform());
	}

	double logProb() {return 0;}

};

class binomial : public Rvar<Int>	{

	static const double InfProb = -1050;
	public:

	binomial(Rvar<Int>* inN, Var<UnitReal>* inp)	{
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
		// cerr << val << '\n';
		setval(val);
	}

	private:

	virtual double 		logProb()	{
		// cerr << "bin log prob\n";
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
	Var<UnitReal>* 	p;
};
 
class PoissonBinModel : public ProbModel	{

	public:

	Dvar<PosReal>* mu;
	Dvar<PosReal>* alpha;
	ConjugateGamma* r;
	ConjugatePoisson* n1;
	ConjugatePoisson* n2;
	binomial* k1;
	binomial* k2;
	uniform* p;

	PoissonBinModel(double inmu, int ink1, int ink2)	{

		mu = new Const<PosReal>(inmu);
		alpha = new Const<PosReal>(1);
		p = new uniform;
		r = new ConjugateGamma(alpha,mu);
		n1 = new ConjugatePoisson(r);
		k1 = new binomial(n1,p);
		n2 = new ConjugatePoisson(r);
		k2 = new binomial(n2,p);
		k1->setval(ink1);
		k1->Clamp();
		k2->setval(ink2);
		k2->Clamp();

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
		// r->Move(5);
		// r->Move(1);
		r->Integrate();
		p->Move(1);
		n1->Move(8);
		n1->Move(2);
		n2->Move(8);
		n2->Move(2);
		r->Resample();
		
		return 1;
	}

	double GetLogProb()	{
		return p->GetLogProb() + r->GetLogProb() + n1->GetLogProb() + k1->GetLogProb() + n2->GetLogProb() + k2->GetLogProb();
	}

	void drawSample()	{
		p->Sample();
		r->Sample();
		while (n1->val() < k1->val())	{
			n1->Sample();
		}
		while (n2->val() < k2->val())	{
			n2->Sample();
		}
	}

	void ToStream(ostream& os)	{
		os << *p << '\t' << *r << '\t' << *n1 << '\t' << *n2  << '\n';
	}

	void FromStream(istream& is)	{
		is >> *p >> *r >> *n1 >> *n2;
	}

	double GetRate()	{
		return r->val();
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *p << '\t' << *r << '\t' << *n1 << '\t' << *n2 << '\n';
	}
};

int main()	{

	double mu = 0.05;
	int k1 = 10;
	int k2 = 20;

	PoissonBinModel* model = new PoissonBinModel(mu,k1,k2);

	int burnin = 1000;
	int nrep = 100000;
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

	cout << "estimated   : " << empmean << '\t' << empvar << '\n';
}

