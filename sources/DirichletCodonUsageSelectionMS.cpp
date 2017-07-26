
#include "Chain.h"
#include "DirichletCodonUsageSelectionModelMS.h"
#include <cmath>


class DirichletCodonUsageSelectionChainMS : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int category;
	int every;
	int burnin;
	int conjugate;
 	string type;
 	string mechanism;

	public:

	DirichletCodonUsageSelectionModelMS* GetModel() {return (DirichletCodonUsageSelectionModelMS*) model;}

	string GetModelType() {return modeltype;}

	DirichletCodonUsageSelectionChainMS(string indata, string intree, int incategory, int inburnin, int inevery, string filename, string intype, int inconjugate, string inmechanism, int force = 0 )	{
		modeltype = "SELECTIONGTR";
		datafile = indata;
		treefile = intree;
		category = incategory;
		burnin = inburnin;
		every = inevery;
		name = filename;
		conjugate = inconjugate;
		type = intype;
		mechanism = inmechanism;

		New(force);

	}


	DirichletCodonUsageSelectionChainMS(string indata, string intree, int incategory, int inevery, string filename, string intype, int inconjugate, string inmechanism, int force = 0 )	{
		modeltype = "SELECTIONGTR";
		datafile = indata;
		treefile = intree;
		category = incategory;
		every = inevery;
		name = filename;
		conjugate = inconjugate;
		type = intype;
		mechanism = inmechanism;

		New(force);

	}



	DirichletCodonUsageSelectionChainMS(string filename, string intype, string inmechanism)	{

		name = filename;
		type = intype;
		mechanism = inmechanism;

		conjugate = 1;
		Open();
		Save();
	}

	void New(int force)	{
		if (modeltype == "SELECTIONGTR")	{
			model = new DirichletCodonUsageSelectionModelMS(datafile,treefile,category,type,conjugate,mechanism);
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
		is >> datafile >> treefile >> category;
		is >> conjugate;
		is >> type >> mechanism;
		is >> every >> until >> size;

		if (modeltype == "SELECTIONGTR")	{
			model = new DirichletCodonUsageSelectionModelMS(datafile,treefile,category,type,conjugate,mechanism);
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
		param_os << datafile << '\t' << treefile << '\t' << category << '\n';
		param_os << conjugate << '\n';
		param_os << type << '\t' << mechanism << '\n';
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
	if (argc == 4)	{
		string name = argv[1];
		string type = argv[2];
		string mechanism = argv[3];

		DirichletCodonUsageSelectionChainMS* chain = new DirichletCodonUsageSelectionChainMS(name,type,mechanism);

		cerr << "start\n";
		chain->Start();
	}

	// this is a new chain
	else if (argc == 9)	{

		string datafile = argv[1];
		string treefile = argv[2];
		int category = atoi(argv[3]);
		int every = atoi(argv[4]);
		string name = argv[5];
		string type = argv[6];
		int conjugate = atoi(argv[7]);
		string mechanism = argv[8];

		cerr << "chain name : " << name << '\n';

		DirichletCodonUsageSelectionChainMS* chain = new DirichletCodonUsageSelectionChainMS(datafile,treefile,category,every,name,type,conjugate,mechanism,1);

		cerr << "start\n";
		chain->Start();
	}
}

