
#include "Chain.h"
#include "PoissonModel.h"
#include "ConjugatePoissonModel.h"

class PoissonChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	bool conjugate;

	public:

	string GetModelType() {return modeltype;}

	PoissonChain(string indata, string filename, int in_until, int in_every, bool inconjugate, int force = 0)	{
		conjugate = inconjugate; 
		if (conjugate)	{
			modeltype = "ConjugatePoisson";
		}
		else	{
			modeltype = "Poisson";
		}
		datafile = indata;
		name = filename;
		until = in_until;
		every = in_every;
		New(force);
	}

	PoissonChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "Poisson")	{
			model = new PoissonModel(datafile);
		}
		else if (modeltype == "ConjugatePoisson")	{
			model = new ConjugatePoissonModel(datafile);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "RESET\n";
		Reset(force);
		cerr << "START, initial lnL = " << model->GetLogProb() << "\n";
		Start();
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile;
		is >> every >> until >> size;

		if (modeltype == "Poisson")	{
			model = new PoissonModel(datafile);
		}
		else if (modeltype == "ConjugatePoisson")	{
			model = new ConjugatePoissonModel(datafile);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		cerr << "UPDATE\n";
		model->Update();
		cerr << "START : " << size << " points saved, current lnL = " << model->GetLogProb() << "\n";
		Start();
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << datafile << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}
	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
			model->Move(0.1);
		}
		SavePoint();
		Save();
		Monitor();
	}

};
