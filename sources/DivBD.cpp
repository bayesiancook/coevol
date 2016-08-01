
#include "Chain.h"
#include "DivBDModel.h"
#include "StringStreamUtils.h"

class DivBDChain : public Chain	{

	private:
	string modeltype;
	string treefile;
	string datafile;
	string calibfile;
	int divmodel;
	int conjpath;
	int burnin;
	double initlambda;
	double finallambda;
	double initmu;
	double finalmu;

	public:

	string GetModelType() {return modeltype;}

	DivBDModel* GetDivModel() {return (DivBDModel*) model;}

	DivBDChain(string intree, string indata, string incalibfile, int indivmodel, double ininitlambda, double infinallambda, double ininitmu, double infinalmu, bool inconjpath, string filename, int force = 1)	{
		modeltype = "DIVBD";
		treefile = intree;
		datafile = indata;
		calibfile = incalibfile;
		divmodel = indivmodel;
		conjpath = inconjpath;
		name = filename;
		burnin = 0;
		initlambda = ininitlambda;
		finallambda = infinallambda;
		initmu = ininitmu;
		finalmu = infinalmu;
		New(force);
	}

	DivBDChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "DIVBD")	{
			model = new DivBDModel(datafile,treefile,calibfile,divmodel,initlambda,finallambda,initmu,finalmu,conjpath,true);
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
		is >> calibfile;
		is >> divmodel;
		is >> initlambda >> finallambda >> initmu >> finalmu;
		is >> conjpath;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> burnin >> every >> until >> size;

		if (modeltype == "DIVBD")	{
			model = new DivBDModel(datafile,treefile,calibfile,divmodel,initlambda,finallambda,initmu,finalmu,conjpath,true);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		cerr << "update\n";
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
		if (burnin)	{
			cerr << "restoring ais\n";
			ResetAIS();
		}
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << treefile << '\t' << datafile << '\n';
		param_os << calibfile << '\n';
		param_os << divmodel << '\n';
		param_os << initlambda << '\t' << finallambda << '\t' << initmu << '\t' << finalmu << '\n';
		param_os << conjpath << '\n';
		param_os << 0 << '\n';
		param_os << burnin << '\t' << every << '\t' << until << '\t' << size << '\n';

		model->ToStream(param_os);
	}

	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
		}
		SavePoint();
		Save();
		Monitor();
		if (burnin)	{
			if (size > burnin)	{
				double p = ((double) (size - burnin) / until);
				double ret = GetDivModel()->AISSet(p);
				ofstream os((name + ".ais").c_str(),ios_base::app);
				os << p << '\t' << GetDivModel()->GetLambda() << '\t' << GetDivModel()->GetMu() << '\t' << ret << '\n';
			}
		}
	}

	void ResetAIS()	{
		// GetDivModel()->MakeAISController();
		if (size > burnin)	{
			double p = ((double) (size - burnin) / until);
			cerr << "restoring ais at p\n";
			GetDivModel()->AISSet(p);
		}
	}

	void InitAIS(int inburnin)	{
		burnin = inburnin;
		ofstream os((name + ".ais").c_str());
		// GetDivModel()->MakeAISController();
		GetDivModel()->AISSet(0);
	}

};

int main(int argc, char* argv[])	{

	if (argc == 2)	{
		string name = argv[1];
		DivBDChain* chain = new DivBDChain(name);
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
		int divmodel = 1;
		bool conjpath = false;

		int burnin = 0;
		int every = 1;
		int until = -1;

		double initlambda = -1;
		double finallambda = -1;
		double initmu = -1;
		double finalmu = -1;

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
				}
				else if (s == "-yule")	{
					divmodel = 0;
				}
				else if (s == "-bd")	{
					divmodel = 1;
				}
				else if (s == "-ais")	{
					i++;
					burnin = atoi(argv[i]);
					i++;
					every = atoi(argv[i]);
					i++;
					until = atoi(argv[i]);
					i++;
					initlambda = atof(argv[i]);
					i++;
					finallambda = atof(argv[i]);
					i++;
					initmu = atof(argv[i]);
					i++;
					finalmu = atof(argv[i]);
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-x")	{
					if (burnin)	{
						throw(0);
					}
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
			if ((datafile == "") || (treefile == "") || (calibfile == "") || (name == ""))	{
				throw(0);
			}
		}
		catch(...)	{
			cerr << "divbd -d <datafile> -t <tree> -c <calibfile> [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		DivBDChain* chain = new DivBDChain(treefile,datafile,calibfile,divmodel,initlambda,finallambda,initmu,finalmu,conjpath,name,force);
		chain->SetEvery(every);
		chain->SetUntil(until);
		if (burnin)	{
			cerr << "init ais\n";
			chain->InitAIS(burnin);
		}
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
}

