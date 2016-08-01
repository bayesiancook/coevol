
#include "Chain.h"
#include "PopSizeAAMutSelModel.h"

class PopSizeAAMutSelChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	PopSizeAAMutSelChain(string indata, string intree, int inP, string filename, GeneticCodeType intype, int force = 1)	{
		modeltype = "POPSIZEAAMUTSEL";
		datafile = indata;
		treefile = intree;
		type = intype;
		P = inP;
		name = filename;
		New(force);
	}

	PopSizeAAMutSelChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "POPSIZEAAMUTSEL")	{
			model = new PopSizeAAMutSelModel(datafile,treefile,P,true,type);
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
		is >> P;
		is >> every >> until >> size;

		if (modeltype == "POPSIZEAAMUTSEL")	{
			model = new PopSizeAAMutSelModel(datafile,treefile,P,false,type);
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
		param_os << P << '\n';
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
		PopSizeAAMutSelChain* chain = new PopSizeAAMutSelChain(name);
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
		
	}

	// this is a new chain
	else	{

		string datafile = argv[1];
		string treefile = argv[2];
		string name = argv[3];
		int P = 20;
		if (argc == 5)	{
			P = atoi(argv[4]);
		}

		GeneticCodeType type = Universal;
		if (argc == 6)	{
			string t = argv[5];
			if (t == "mtmam")	{
				type = MtMam;
			}
			else if (t == "uni")	{
				type = Universal;
			}
		}
		PopSizeAAMutSelChain* chain = new PopSizeAAMutSelChain(datafile,treefile,P,name,type);
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
	}
}
