
#include "Chain.h"
#include "BranchOmegaIIDModel.h"

class BranchOmegaIIDChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	BranchOmegaIIDChain(string indata, string intree, string filename, int force = 0)	{
		modeltype = "BRANCHOMEGAIID";
		datafile = indata;
		treefile = intree;
		name = filename;
		New(force);
	}

	BranchOmegaIIDChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BRANCHOMEGAIID")	{
			model = new BranchOmegaIIDModel(datafile,treefile);
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
		is >> datafile >> treefile;
		is >> every >> until >> size;

		if (modeltype == "BRANCHOMEGA")	{
			model = new BranchOmegaIIDModel(datafile,treefile);
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

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	BranchOmegaIIDChain* chain = new BranchOmegaIIDChain(datafile,treefile,name);
	cerr << "start\n";
	chain->Run();
}
