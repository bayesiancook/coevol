
#include "Chain.h"
#include "BranchGenOmegaModel.h"

class BranchGenOmegaChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string contdatafile;
	string treefile;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	BranchGenOmegaChain(string indata, string intree, string incontdata, string filename, GeneticCodeType intype, int force = 1)	{
		modeltype = "BRANCHGENOMEGA";
		type = intype;
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		name = filename;
		New(force);
	}

	BranchGenOmegaChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BRANCHGENOMEGA")	{
			model = new BranchGenOmegaModel(datafile,treefile,contdatafile);
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

		if (modeltype == "BRANCHGENOMEGA")	{
			model = new BranchGenOmegaModel(datafile,treefile,contdatafile);
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

	string datafile = argv[1];
	string treefile = argv[2];
	string contdatafile = argv[3];
	string name = argv[4];
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

	BranchGenOmegaChain* chain = new BranchGenOmegaChain(datafile,treefile,contdatafile,name,type);
	cerr << "start\n";
	chain->Start();
	cerr << "exit\n";
}
