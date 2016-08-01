
#include "Chain.h"
#include "BranchOmegaGCModel.h"

class BranchOmegaGCChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	BranchOmegaGCChain(string indata, string intree, string incont, string filename, GeneticCodeType intype, int force = 1)	{
		modeltype = "BRANCHOMEGAGC";
		datafile = indata;
		treefile = intree;
		contdatafile = incont;
		type = intype;
		name = filename;
		New(force);
	}

	BranchOmegaGCChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BRANCHOMEGAGC")	{
			model = new BranchOmegaGCModel(datafile,treefile,contdatafile,true,type);
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
		is >> contdatafile;
		is >> every >> until >> size;

		if (modeltype == "BRANCHOMEGAGC")	{
			model = new BranchOmegaGCModel(datafile,treefile,contdatafile,false,type);
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
		param_os << contdatafile << '\n';
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
	string contdatafile = argv[3];
	string name = argv[4];

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
	BranchOmegaGCChain* chain = new BranchOmegaGCChain(datafile,treefile,contdatafile,name,type);
	cerr << "start\n";
	chain->Start();
	cerr << "exit\n";
}
