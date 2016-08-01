
#include "Chain.h"
// #include "SeparateMultiVariateAAMutSelMatMixChronoModel.h"
// #include "SeparateMultiVariateAAMutSelMatMixModel.h"
#include "DirSeparateMultiVariateAAMutSelMatMixModel.h"

class SeparateMultiVariateAAMutSelChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	bool withgc;
	bool withpopsize;
	bool withaacomp;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	SeparateMultiVariateAAMutSelChain(string indata, string intree, int inP, bool inwithgc, bool inwithpopsize, bool inwithaacomp, string filename, GeneticCodeType intype, int force = 1)	{
		modeltype = "SEPARATEMULTIVARIATEAAMUTSEL";
		datafile = indata;
		treefile = intree;
		type = intype;
		P = inP;
		withgc = inwithgc;
		withpopsize = inwithpopsize;
		withaacomp = inwithaacomp;
		name = filename;
		New(force);
	}

	SeparateMultiVariateAAMutSelChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "SEPARATEMULTIVARIATEAAMUTSEL")	{
			model = new SeparateMultiVariateAAMutSelModel(datafile,treefile,P,withgc,withpopsize,withaacomp,true,type);
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
		is >> withgc >> withpopsize >> withaacomp;
		is >> every >> until >> size;

		if (modeltype == "SEPARATEMULTIVARIATEAAMUTSEL")	{
			model = new SeparateMultiVariateAAMutSelModel(datafile,treefile,P,withgc,withpopsize,withaacomp,false,type);
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
		param_os << withgc << '\t' << withpopsize << '\t' << withaacomp << '\n';
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

	string name = "";
	string datafile = "";
	string treefile = "";
	int P = 20;
	bool withgc = true;
	bool withpopsize = true;
	bool withaacomp = true;
	GeneticCodeType type = Universal;

	
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
			else if (s == "-ncat")	{
				i++;
				P = atoi(argv[i]);
			}
			else if ((s == "-acgt") || (s == "-gc"))	{
				withgc = true;
			}
			else if (s == "-popsize")	{
				withpopsize = true;
			}
			else if (s == "-aa")	{
				withaacomp = true;
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if (s == "-mtmam")	{
				type = MtMam;
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	SeparateMultiVariateAAMutSelChain* chain;
	if (datafile == "")	{
		chain = new SeparateMultiVariateAAMutSelChain(name);
	}
	else	{
		chain = new SeparateMultiVariateAAMutSelChain(datafile,treefile,P,withgc,withpopsize,withaacomp,name,type);
	}

	cerr << "start\n";
	chain->Start();
	cerr << "exit\n";
}
