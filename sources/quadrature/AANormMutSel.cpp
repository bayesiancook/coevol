
#include "Chain.h"
#include "AANormMutSelModel.h"

class AANormMutSelChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	AANormMutSelChain(string indata, string intree, string filename, GeneticCodeType intype, int force = 1)	{
		modeltype = "AANORMMUTSEL";
		datafile = indata;
		treefile = intree;
		type = intype;
		name = filename;
		New(force);
	}

	AANormMutSelChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "AANORMMUTSEL")	{
			model = new AANormMutSelModel(datafile,treefile,true,type);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "RESET\n";
		Reset(force);
		cerr << "START\n";
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile;
		is >> every >> until >> size;

		if (modeltype == "AANORMMUTSEL")	{
			model = new AANormMutSelModel(datafile,treefile,false,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << type << '\n';
		param_os << datafile << '\t' << treefile << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}

	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
		}
		SavePoint();
		Save();
		Monitor();
	}
};

int main(int argc, char* argv[])	{

	// this is an already existing chain on the disk; reopen and restart
	if (argc == 2)	{
		string name = argv[1];
		AANormMutSelChain* chain = new AANormMutSelChain(name);
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
		
	}

	// this is a new chain
	else	{

		string datafile = argv[1];
		string treefile = argv[2];
		string name = argv[3];

		GeneticCodeType type = Universal;
		if (argc == 5)	{
			string t = argv[4];
			if (t == "mtmam")	{
				type = MtMam;
			}
			else if (t == "uni")	{
				type = Universal;
			}
		}
		AANormMutSelChain* chain = new AANormMutSelChain(datafile,treefile,name,type);
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
	}
}
