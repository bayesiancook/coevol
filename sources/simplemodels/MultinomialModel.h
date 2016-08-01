
#include "RandomTypes.h"
#include "ProbModel.h"

class MultinomialModel : public ProbModel	{

	public:

	int N;
	int dim;
	int** data;
	int* totcount;

	Dvar<PosReal>* Hyper;
	Exponential* concentration;
	Dirichlet* center;
	Dirichlet** rates;
	Multinomial** counts;


	MultinomialModel(string datafile)	{
	
		// read data
		ifstream is(datafile.c_str());
		is >> N;
		is >> dim;
		data = new int*[N];
		totcount = new int[N];
		for (int i=0; i<N; i++)	{
			data[i] = new int[dim];
		}

		for (int i=0; i<N; i++)	{
			totcount[i] = 0;
			for (int k=0; k<dim; k++)	{
				is >> data[i][k];
				totcount[i] += data[i][k];
			}
		}

		Hyper = new Const<PosReal>(PosReal(0.05));
		concentration = new Exponential(Hyper);
		// concentration->setval(20);
		center = new Dirichlet(dim);
		
		rates = new Dirichlet*[N];
		counts = new Multinomial*[N];
		
		for (int i=0; i<N; i++)	{
			// rates[i] = new Dirichlet(dim);
			rates[i] = new Dirichlet(center,concentration);
			counts[i] = new Multinomial(rates[i],totcount[i]);
		}

		RootRegister(Hyper);
		RootRegister(center);
		/*
		for (int i=0; i<N; i++)	{
			RootRegister(rates[i]);
		}
		*/
		Register();

		ClampObservedNodes();

		cerr << "model created\n";
	}

	~MultinomialModel()	{
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
		if (isnan(total))	{
			cerr << "error\n";
			cerr << center->GetLogProb() << '\n';
			cerr << concentration->GetLogProb() << '\n';
			for (int i=0; i<N; i++)	{
				cerr << rates[i]->GetLogProb() << '\t';
				cerr << *rates[i] << '\t';
				cerr << counts[i]->GetLogProb() << '\n';
			}
			exit(1);
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
		concentration->Move(0.1*tuning);
		concentration->Move(0.01*tuning);
		concentration->Move(0.001*tuning);
		center->Move(0.01*tuning);
		center->Move(0.001*tuning);
		for (int i=0; i<N; i++)	{
			rates[i]->Move(tuning);
			rates[i]->Move(tuning,1);
		}	
		return 1;
	}

	void drawSample()	{
		concentration->Sample();
		center->Sample();
		for (int i=0; i<N; i++)	{
			rates[i]->Sample();
		}	

		// will not do anything if counts nodes clamped
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

