
#include "Chain.h"
#include "BetaModel.h"

class BetaChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	int burnin;

	public:

	string GetModelType() {return modeltype;}

	BetaChain(string indata, string filename, int in_burnin, int in_until, int force = 0)	{
		modeltype = "Beta";
		datafile = indata;
		name = filename;
		burnin = in_burnin;
		until = in_until;
		every = 10;
		New(force);
	}

	BetaChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "Beta")	{
			model = new BetaModel(datafile);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "RESET\n";
		Reset(force);
		cerr << "START, initial lnL = " << model->GetLogProb() << "\n";
		Start();
		cerr << "chain done\n";
		ofstream os((name + ".meantheta").c_str());
		((BetaModel*) model)->Normalise();
		((BetaModel*) model)->PrintEstimates(os);
		os.close();
		cerr << "estimates in " << name << ".meantheta\n";
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

		if (modeltype == "Beta")	{
			model = new BetaModel(datafile);
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
		((BetaModel*) model)->AddUp();
	}
};
