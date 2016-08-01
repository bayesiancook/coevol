
#include "Chain.h"
#include "GTRModel.h"
#include "ConjugateGTRModel.h"

class GTRChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int mh;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	GTRChain(string indata, string intree, string filename, int inmh, int force = 0)	{
		modeltype = "GTR";
		datafile = indata;
		treefile = intree;
		mh = inmh;
		name = filename;
		New(force);
	}

	GTRChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "GTR")	{
			model = new GTRModel(datafile,treefile,mh);
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
		is >> mh;
		is >> every >> until >> size;

		if (modeltype == "GTR")	{
			model = new GTRModel(datafile,treefile,mh);
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
		param_os << mh << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}
	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
			// model->Move(0.1);
		}
		SavePoint();
		Save();
		Monitor();
	}
};


int main(int argc, char* argv[])	{

	if (argc == 2)	{
		string name = argv[1];
		GTRChain* chain = new GTRChain(name);
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
		bool mh = false;

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
				else if (s == "-mh")	{
					mh = 1;
				}
				else if (s == "-nonmh")	{
					mh = 0;
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
			cerr << "gtr -d <datafile> -t <tree> [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		GTRChain* chain = new GTRChain(datafile,treefile,name,mh,force);
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

