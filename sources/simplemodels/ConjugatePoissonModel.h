
#include "RandomTypes.h"
#include "ProbModel.h"
#include "Conjugate.h"
#include "ConjugatePoisson.h"

class ConjugatePoissonModel : public ProbModel	{

	public:

	int N;
	int* data;

	Dvar<PosReal>* Hyper;
	Exponential* alpha;
	Exponential* beta;
	ConjugateGamma** rates;
	ConjugatePoisson** counts;


	ConjugatePoissonModel(string datafile)	{
	
		cerr << "create conjugate poisson\n";
		// read data
		ifstream is(datafile.c_str());
		is >> N;
		data = new int[N];
		for (int i=0; i<N; i++)	{
			is >> data[i];
		}

		Hyper = new Const<PosReal>(PosReal(1));
		alpha = new Exponential(Hyper,Exponential::RATE);
		beta = new Exponential(Hyper,Exponential::RATE);
		
		rates = new ConjugateGamma*[N];
		counts = new ConjugatePoisson*[N];
		
		for (int i=0; i<N; i++)	{
			rates[i] = new ConjugateGamma(alpha,beta);
			counts[i] = new ConjugatePoisson(rates[i]);
		}

		RootRegister(Hyper);
		Register();

		ClampObservedNodes();

		cerr << "model created\n";
	}

	~ConjugatePoissonModel()	{
		for (int i=0; i<N; i++)	{
			delete rates[i];
			delete counts[i];
		}
		delete[] rates;
		delete[] counts;
		delete[] data;
	
		delete alpha;
		delete beta;
		delete Hyper;
	}

	void MakeScheduler() {}

	double GetLogProb()	{
		double total = 0;
		total += alpha->GetLogProb();
		total += beta->GetLogProb();
		for (int i=0; i<N; i++)	{
			total += rates[i]->GetLogProb();
			total += counts[i]->GetLogProb();
		}
		return total;
	}

	double GetMeanRate()	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += rates[i]->val();
		}
		return total / N;
	}
		
	double GetMeanCount()	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += counts[i]->val();
		}
		return total / N;
	}
		
	void ClampObservedNodes()	{
		for (int i=0; i<N; i++)	{
			counts[i]->ClampAt(data[i]);
		}
	}
	
	void TraceHeader(ostream& os)	{
		os << "#logprob\talpha\tbeta\tmeanrate\tmeancounts\n";
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *alpha << '\t' << *beta << '\t' << GetMeanRate() << '\t' << GetMeanCount();
		os << '\n';
	}

	double Move(double tuning)	{

		for (int i=0; i<N; i++)	{
			// rates[i]->ActivateSufficientStatistic();
			rates[i]->Integrate();
		}	
		
		alpha->Move(tuning);
		beta->Move(tuning);
		for (int i=0; i<N; i++)	{
			rates[i]->InactivateSufficientStatistic();
		}	
		
		/*
		for (int i=0; i<N; i++)	{
			rates[i]->Move(tuning);
		}	
		*/
		return 1;
	}

	void drawSample()	{
		alpha->Sample();
		beta->Sample();
		for (int i=0; i<N; i++)	{
			rates[i]->Sample();
		}	

		// will not do anything if counts nodes are clamped
		for (int i=0; i<N; i++)	{
			counts[i]->Sample();
		}
	}

	void ToStream(ostream& os)	{

		os << *alpha << '\n';
		os << *beta << '\n';
		for (int i=0; i<N; i++)	{
			os << *rates[i] << '\t';
		}
		os << '\n';
	}
	
	void FromStream(istream& is)	{

		is >> *alpha;
		is >> *beta;
		for (int i=0; i<N; i++)	{
			is >> *rates[i];
		}
	}

};

