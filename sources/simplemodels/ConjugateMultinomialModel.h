
#include "RandomTypes.h"
#include "ProbModel.h"
#include "Conjugate.h"

class ConjugateMultinomialModel : public ProbModel	{

	public:

	int N;
	int dim;
	int** data;
	int* totcount;

	Dvar<PosReal>* Hyper;
	Exponential* concentration;
	ConjugateDirichlet* center;
	ConjugateDirichlet** rates;
	ConjugateMultinomial** counts;


	ConjugateMultinomialModel(string datafile)	{
	
		// read data
		ifstream is(datafile.c_str());
		is >> N;
		is >> dim;
		data = new int*[N];
		totcount = new int[N];
		for (int i=0; i<N; i++)	{
			data[i] = new int[dim];
			totcount[i] = 0;
		}

		for (int i=0; i<N; i++)	{
			for (int k=0; k<dim; k++)	{
				is >> data[i][k];
				totcount[i] += data[i][k];
			}
		}

		Hyper = new Const<PosReal>(PosReal(0.05));
		concentration = new Exponential(Hyper);
		// concentration->setval(20);
		center = new ConjugateDirichlet(dim);
		
		rates = new ConjugateDirichlet*[N];
		counts = new ConjugateMultinomial*[N];
		
		for (int i=0; i<N; i++)	{
			rates[i] = new ConjugateDirichlet(center,concentration);
			counts[i] = new ConjugateMultinomial(rates[i],totcount[i]);
		}

		RootRegister(Hyper);
		RootRegister(center);
		for (int i=0; i<N; i++)	{
			RootRegister(rates[i]);
		}
		Register();

		ClampObservedNodes();

		cerr << "model created\n";
	}

	~ConjugateMultinomialModel()	{
		for (int i=0; i<N; i++)	{
			delete[] data[i];
			delete rates[i];
			delete counts[i];
		}
		delete[] rates;
		delete[] counts;
		delete[] data;
		delete[] totcount;
	
		delete center;
		delete concentration;
		delete Hyper;
	}

	double GetLogProb()	{
		double total = 0;
		total += center->GetLogProb();
		total += concentration->GetLogProb();
		for (int i=0; i<N; i++)	{
			total += rates[i]->GetLogProb();
			total += counts[i]->GetLogProb();
		}
		return total;
	}

	double GetMeanEntropy()	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += rates[i]->val().GetEntropy();
		}
		total /= N;
		return total;
	}

	void ClampObservedNodes()	{
		for (int i=0; i<N; i++)	{
			counts[i]->ClampAt(data[i]);
		}
	}
	
	void TraceHeader(ostream& os)	{
		os << "#logprob\tconcentration\tcenter\tmeanentropy\n";
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << *concentration << '\t' << center->val().GetEntropy() << '\t' << GetMeanEntropy() << '\n';
	}

	double Move(double tuning)	{
		for (int i=0; i<N; i++)	{
			rates[i]->Integrate();
		}	
		
		concentration->Move(0.1*tuning);
		concentration->Move(0.01*tuning);
		concentration->Move(0.001*tuning);
		center->Move(1*tuning,1);
		center->Move(0.1*tuning,1);
		center->Move(0.1*tuning);
		center->Move(0.01*tuning);
		center->Move(0.001*tuning);

		for (int i=0; i<N; i++)	{
			rates[i]->InactivateSufficientStatistic();
		}	
		concentration->Move(0.1*tuning);
		concentration->Move(0.01*tuning);
		concentration->Move(0.001*tuning);
		center->Move(0.01*tuning);
		center->Move(0.001*tuning);
		for (int i=0; i<N; i++)	{
			rates[i]->Move(tuning);
		}
		return 1;
	}

	void drawSample()	{
		concentration->Sample();
		center->Sample();
		for (int i=0; i<N; i++)	{
			rates[i]->Sample();
		}	

		for (int i=0; i<N; i++)	{
			counts[i]->Sample();
		}
	}

	void ToStream(ostream& os)	{

		os << *concentration << '\n';
		os << *center << '\n';
		for (int i=0; i<N; i++)	{
			os << *rates[i] << '\t';
		}
		os << '\n';
	}
	
	void FromStream(istream& is)	{

		is >> *concentration;
		is >> *center;
		for (int i=0; i<N; i++)	{
			is >> *rates[i];
		}
	}

};

