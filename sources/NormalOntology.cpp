
#include "Chain.h"
#include "NormalOntologyModel.h"
#include "StringStreamUtils.h"

class NormalOntologyChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string ontologyfile;
	int mingene;
	int minconcept;
	int withtoggle;
	double logvar;

	public:

	string GetModelType() {return modeltype;}

	NormalOntologyChain(string indatafile, string inontologyfile, int inmingene, int inminconcept, int inwithtoggle, double inlogvar, string inname, int force = 0)	{
		modeltype = "NORMALONTOLOGY";
		datafile = indatafile;
		ontologyfile = inontologyfile;
		mingene = inmingene;
		minconcept = inminconcept;
		withtoggle = inwithtoggle;
		logvar = inlogvar;
		name = inname;
		New(force);
	}

	NormalOntologyChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "NORMALONTOLOGY")	{
			model = new NormalOntologyModel(datafile,ontologyfile,mingene,minconcept,withtoggle,logvar);
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
		withtoggle = 0;
		logvar = -100;
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile >> ontologyfile;
		is >> mingene >> minconcept;
		int check;
		is >> check;
		if (check)	{
			is >> withtoggle;
			is >> check;
			if (check)	{
				is >> logvar;
				is >> check;
				if (check)	{
					cerr << "error when reading model\n";
					exit(1);
				}
			}
		}

		is >> every >> until >> size;

		if (modeltype == "NORMALONTOLOGY")	{
			model = new NormalOntologyModel(datafile,ontologyfile,mingene,minconcept,withtoggle,logvar);
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
		param_os << datafile << '\t' << ontologyfile << '\n';
		param_os << mingene << '\t' << minconcept << '\n';
		param_os << 1 << '\n';
		param_os << withtoggle << '\n';
		param_os << 1 << '\n';
		param_os << logvar << '\n';
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
		NormalOntologyChain* chain = new NormalOntologyChain(name);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		string datafile = "";
		string ontologyfile = "";
		int mingene = 1;
		int minconcept = 1;
		int force = 0;
		int withtoggle = 0;
		double logvar = -100;

		int every = 1;
		int until = -1;

		string name = "";
		// int nrep = 0;

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
				else if (s == "-o")	{
					i++;
					ontologyfile = argv[i];
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-mingene")	{
					i++;
					mingene = atoi(argv[i]);
				}
				else if (s == "-minconcept")	{
					i++;
					minconcept = atoi(argv[i]);
				}
				else if (s == "-withtoggle")	{
					withtoggle = 1;
				}
				else if (s == "-logvar")	{
					i++;
					logvar = atof(argv[i]);
				}
				else if ( (s == "-x") || (s == "-extract") )	{
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
			if ((datafile == "") || (ontologyfile == "") || (name == ""))	{
				throw(0);
			}
		}
		catch(...)	{
			cerr << "normalontology -d <data> -o <ontology> [-mingene <mingene> -minconcept <minconcept>]\n";
			cerr << '\n';
			exit(1);
		}

		NormalOntologyChain* chain = new NormalOntologyChain(datafile,ontologyfile,mingene,minconcept,withtoggle,logvar,name,force);
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

