
#include "Chain.h"
#include "DolloModel.h"
#include "StringStreamUtils.h"

class DolloChain : public Chain	{

	private:
	string modeltype;
	string treefile;
	string datafile;
	string calibfile;
	double rootage;
	double rootstdev;

	public:

	string GetModelType() {return modeltype;}

	DolloChain(string intree, string indata, string incalibfile, double inrootage, double inrootstdev, string filename, int force = 1)	{
		modeltype = "DOLLO";
		treefile = intree;
		datafile = indata;
		calibfile = incalibfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		name = filename;
		New(force);
	}

	DolloChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "DOLLO")	{
			model = new DolloLogNormalModel(datafile,treefile,calibfile,rootage,rootstdev,true);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		Reset(force);
		// model->Update();
		cerr << "ln prob = " << model->GetLogProb() << "\n";
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> treefile >> datafile;
		is >> calibfile >> rootage >> rootstdev;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "LOGN")	{
			model = new DolloLogNormalModel(datafile,treefile,calibfile,rootage,rootstdev,true);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		cerr << "update\n";
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << treefile << '\t' << datafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << 0 << '\n';
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

	if (argc == 2)	{
		string name = argv[1];
		DolloChain* chain = new DolloChain(name);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		string treefile = "";
		string datafile = "";
		string name = "";
		string calibfile = "None";
		string rootfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int every = 1;
		int until = -1;

		int force = 0;

		try	{

			if (argc == 1)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

				if (s == "-d")	{
					i++;
					datafile = argv[i];
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-cal")	{
					i++;
					calibfile = argv[i];
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootage = atof(argv[i]);
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootstdev = atof(argv[i]);
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-x")	{
					i++;
					if (i == argc) throw(0);
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					until = atoi(argv[i]);
				}
				else	{
					if (i != (argc -1))	{
						throw(0);
					}
					name = argv[i];
				}
				i++;
			}
			if ((datafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
		}
		catch(...)	{
			cerr << "dollo -d <datafile> -t <tree> [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		DolloChain* chain = new DolloChain(treefile,datafile,calibfile,rootage,rootstdev,name,force);
		chain->SetEvery(every);
		chain->SetUntil(until);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
}

