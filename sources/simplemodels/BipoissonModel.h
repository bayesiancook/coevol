
#include "RandomTypes.h"
#include "ProbModel.h"
#include "Bipoisson.h"

class BipoissonModel : public ProbModel	{

	public:

	int N;
	int P;
	int** data;

	Dvar<PosReal>* Hyper;
	Dvar<PosReal>* Hyper2;
	Exponential* alpha;
	Exponential* lambda;
	Exponential* mu;
	Gamma** rates;
	Gamma** lengths;
	Bipoisson*** counts;


	BipoissonModel(string datafile)	{
	
		// read data
		ifstream is(datafile.c_str());
		int inN, inP;
		is >> inN >> inP;
		MakeNewModel(inN, inP);

		data = new int*[N];
		for (int i=0; i<N; i++)	{
			data[i] = new int[P];
			for (int j=0; j<P; j++)	{
				is >> data[i][j];
			}
		}
		ClampObservedNodes();
	}

	BipoissonModel(int inN, int inP)	{
		MakeNewModel(inN, inP);
	}

	void MakeNewModel(int inN, int inP)	{
		N = inN;
		P = inP;
		data = 0;

		Hyper = new Const<PosReal>(PosReal(0.02));
		Hyper2 = new Const<PosReal>(PosReal(0.01));
		alpha = new Exponential(Hyper,Exponential::RATE);
		lambda = new Exponential(Hyper,Exponential::RATE);
		mu = new Exponential(Hyper2,Exponential::RATE);
		
		rates = new Gamma*[N];
		for (int i=0; i<N; i++)	{
			// rates[i] = new Gamma(alpha,alpha);
			rates[i] = new Gamma(alpha,alpha);
		}

		lengths = new Gamma*[P];
		for (int j=0; j<P; j++)	{
			lengths[j] = new Gamma(lambda,mu);
		}

		counts = new Bipoisson**[N];
		for (int i=0; i<N; i++)	{
			counts[i] = new Bipoisson*[P];
			for (int j=0; j<P; j++)	{
				counts[i][j] = new Bipoisson(lengths[j],rates[i]);
			}
		}
		
		RootRegister(Hyper);
		RootRegister(Hyper2);
		Register();
		cerr << "model created\n";
	}

	~BipoissonModel()	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				delete counts[i][j];
			}
			delete[] counts[i];
			if (data)	{
				delete[] data[i];
			}
		}
		delete[] counts;
		delete[] data;

		for (int i=0; i<N; i++)	{
			delete rates[i];
		}
		delete[] rates;
	
		for (int j=0; j<P; j++)	{
			delete lengths[j];
		}
		delete[] lengths;
	
		delete alpha;
		delete lambda;
		delete mu;
		delete Hyper;
	}

	double GetLogProb()	{
		double total = 0;
		total += alpha->GetLogProb();
		total += lambda->GetLogProb();
		total += mu->GetLogProb();
		for (int i=0; i<N; i++)	{
			total += rates[i]->GetLogProb();
		}
		for (int j=0; j<P; j++)	{
			total += lengths[j]->GetLogProb();
		}

		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				total += counts[i][j]->GetLogProb();
			}
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
		
	double GetMeanLength()	{
		double total = 0;
		for (int j=0; j<P; j++)	{
			total += lengths[j]->val();
		}
		return total / P;
	}
		
	void ClampObservedNodes()	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				counts[i][j]->ClampAt(data[i][j]);
			}
		}
	}
	
	void TraceHeader(ostream& os)	{
		os << "#logprob\talpha\tlambda\tmu\tmeanrate\tmeanlength\n";
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *alpha << '\t' << *lambda << '\t' << *mu << '\t' << GetMeanRate() << '\t' << GetMeanLength();
		os << '\n';
	}

	void MakeScheduler() {}

	double Move(double tuning)	{
		int nrep = 10;
		for (int rep=0; rep<nrep; rep++)	{
			alpha->Move(tuning);
		}
		for (int rep=0; rep<nrep; rep++)	{
			lambda->Move(tuning);
			mu->Move(tuning);
			lambda->Move(5*tuning);
			mu->Move(5*tuning);
		}
		for (int i=0; i<N; i++)	{
			rates[i]->Move(tuning);
		}	
		for (int j=0; j<P; j++)	{
			lengths[j]->Move(tuning);
		}
		return 1;
	}

	void drawSample()	{
		alpha->Sample();
		lambda->Sample();
		mu->Sample();
		for (int i=0; i<N; i++)	{
			rates[i]->Sample();
		}	
		for (int j=0; j<P; j++)	{
			lengths[j]->Sample();
		}

		// this will usually not take place since they are clamped
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				counts[i][j]->Sample();
			}
		}
	}

	void Simulation(string name)	{
		Sample();
		ofstream os(name.c_str());
		os << N << '\t' << P << '\n';
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				os << *counts[i][j] << '\t';
			}
			os << '\n';
		}
		os.close();
		ofstream os2((name + ".param").c_str());
		ToStream(os2);
	}

	
	void ToStream(ostream& os)	{

		os << *alpha << '\n';
		os << *lambda << '\n';
		os << *mu << '\n';
		for (int i=0; i<N; i++)	{
			os << *rates[i] << '\t';
		}
		os << '\n';
		for (int j=0; j<P; j++)	{
			os << *lengths[j] << '\t';
		}
		os << '\n';
	}
	
	void FromStream(istream& is)	{

		is >> *alpha;
		is >> *lambda;
		is >> *mu;
		for (int i=0; i<N; i++)	{
			is >> *rates[i];
		}
		for (int j=0; j<P; j++)	{
			is >> *lengths[j];
		}
	}
};

