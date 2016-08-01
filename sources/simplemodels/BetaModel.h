
#include "RandomTypes.h"
#include "ProbModel.h"

class BetaModel : public ProbModel	{

	public:

	int N;
	int* k;
	int* n;

	Dvar<PosReal>* Hyper;
	Dvar<PosReal>* Hyper2;
	Exponential* alpha;
	Exponential* beta;
	Beta** theta;
	Binomial** counts;
	double* meantheta;
	double* vartheta;
	int size;


	BetaModel(string datafile, double hyper = 1)	{
	
		Hyper2 = 0;
		// read data
		ifstream is(datafile.c_str());
		is >> N;
		k = new int[N];
		n = new int[N];

		meantheta = new double[N];
		vartheta = new double[N];
		for (int i=0; i<N; i++)	{
			meantheta[i] = 0;
			vartheta[i] = 0;
		}
		size = 0;

		for (int i=0; i<N; i++)	{
			is >> k[i] >> n[i];
		}

		Hyper = new Const<PosReal>(PosReal(hyper));
		alpha = new Exponential(Hyper,Exponential::RATE);
		beta = new Exponential(Hyper,Exponential::RATE);
		
		theta = new Beta*[N];
		counts = new Binomial*[N];
		
		for (int i=0; i<N; i++)	{
			theta[i] = new Beta(alpha,beta);
			counts[i] = new Binomial(n[i],theta[i]);
		}

		RootRegister(Hyper);
		Register();

		ClampObservedNodes();
	
		cerr << "model created\n";
	}

	BetaModel(int inN, int inn, double hyper = 1)	{
	
		N = inN;
		k = new int[N];
		n = new int[N];
		for (int i=0; i<N; i++)	{
			n[i] = inn;
		}

		Hyper = new Const<PosReal>(PosReal(hyper));
		Hyper2 = new Const<PosReal>(PosReal(hyper));
		alpha = new Exponential(Hyper,Exponential::RATE);
		beta = new Exponential(Hyper,Exponential::RATE);
		
		theta = new Beta*[N];
		counts = new Binomial*[N];
		
		for (int i=0; i<N; i++)	{
			theta[i] = new Beta(Hyper,Hyper2);
			counts[i] = new Binomial(n[i],theta[i]);
		}

		RootRegister(Hyper);
		Register();

		cerr << "model created\n";
		Sample();
	}

	~BetaModel()	{
		/*
		for (int i=0; i<N; i++)	{
			delete theta[i];
			delete counts[i];
		}
		delete[] theta;
		delete[] counts;
		delete[] n;
		delete[] k;
	
		delete alpha;
		delete beta;
		delete Hyper;
		delete Hyper2;
		*/
	}

	void MakeScheduler() {}

	double GetLogProb()	{
		double total = 0;
		total += alpha->GetLogProb();
		total += beta->GetLogProb();
		for (int i=0; i<N; i++)	{
			total += theta[i]->GetLogProb();
			total += counts[i]->GetLogProb();
		}
		return total;
	}

	double GetMeanTheta()	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += theta[i]->val();
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
		
	void AddUp()	{
		for (int i=0; i<N; i++)	{
			meantheta[i] += theta[i]->val();
			vartheta[i] += theta[i]->val() * theta[i]->val();
		}
		size++;
	}
	
	void Normalise()	{
		for (int i=0; i<N; i++)	{
			meantheta[i] /= size;
			vartheta[i] /= size;
			vartheta[i] -= meantheta[i] * meantheta[i];
		}
	}
		
	void PrintEstimates(ostream& os)	{
		for (int i=0; i<N; i++)	{
			os << meantheta[i] << '\t' << sqrt(vartheta[i]) << '\n';
		}
	}
		
	void ClampObservedNodes()	{
		for (int i=0; i<N; i++)	{
			counts[i]->ClampAt(k[i]);
		}
	}
	
	void TraceHeader(ostream& os)	{
		os << "#logprob\talpha\tbeta\tmeantheta\tmeancounts\n";
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *alpha << '\t' << *beta << '\t' << GetMeanTheta() << '\t' << GetMeanCount();
		os << '\n';
	}

	double Move(double tuning)	{
		alpha->Move(tuning);
		beta->Move(tuning);
		for (int i=0; i<N; i++)	{
			theta[i]->Move(tuning);
		}
		return 1;
	}

	
	void drawSample()	{
		alpha->Sample();
		beta->Sample();
		for (int i=0; i<N; i++)	{
			theta[i]->Sample();
		}	
		for (int i=0; i<N; i++)	{
			counts[i]->Sample();
		}
	}

	void ToStream(ostream& os)	{

		os << *alpha << '\t';
		os << *beta << '\t';
		for (int i=0; i<N; i++)	{
			os << *theta[i] << '\t';
		}
		os << '\n';
	}
	
	void FromStream(istream& is)	{

		is >> *alpha;
		is >> *beta;
		for (int i=0; i<N; i++)	{
			is >> *theta[i];
		}
	}

	void PrintSimu(ostream& os)	{
		for (int i=0; i<N; i++)	{
			os << *theta[i] << '\t' << *counts[i] << '\t' << n[i] << '\n';
		}
	}

	void PrintData(ostream& os)	{
		os << N << '\n';
		os << '\n';
		for (int i=0; i<N; i++)	{
			os << *counts[i] << '\t' << n[i] << '\n';
		}
	}

};

