
#include "Chain.h"
#include "CodonAAProfileMutSelMatInfMixModel.h"

class CodonAAProfileMutSelInfMixChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	//ProbModel* GetModel() {return (ProbModel*) model;}
	CodonAAProfileMutSelMatInfMixModel* GetModel() {return (CodonAAProfileMutSelMatInfMixModel*) model;}

	string GetModelType() {return modeltype;}

	CodonAAProfileMutSelInfMixChain(string indata, string intree, int inP, string filename, GeneticCodeType intype, int inevery = 1, int inuntil = -1, int force = 1)	{
		modeltype = "CODONAAPROFILEINFMIXMUTSEL";
		datafile = indata;
		treefile = intree;
		type = intype;
		P = inP;
		name = filename;
		every = inevery;
		until = inuntil;
		New(force);
	}

	CodonAAProfileMutSelInfMixChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "CODONAAPROFILEINFMIXMUTSEL")	{
			model = new CodonAAProfileMutSelMatInfMixModel(datafile,treefile,P,true,type);
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

		if (modeltype == "CODONAAPROFILEINFMIXMUTSEL")	{
			model = new CodonAAProfileMutSelMatInfMixModel(datafile,treefile,P,false,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		GetModel()->FromStreamWithUpdate(is);
		//model->Update();
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

	void Gnuplot()	{
		ofstream gp_os((name + ".gnuplot").c_str());
		gp_os << "p '" << (name + ".trace").c_str() << "' u 1 title 'lnPrior'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 2 title 'lnLike'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 3 title 'ncomp'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 4 title 'alpha'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 5 title 'centerent'\n" << "pause -1\n";
		/*
		gp_os << "p '" << (name + ".trace").c_str() << "' u 7 title 'centerA', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 8 title 'centerC', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 9 title 'centerD', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 10 title 'centerE', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 11 title 'centerF', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 12 title 'centerG', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 13 title 'centerH', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 14 title 'centerI', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 15 title 'centerK', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 16 title 'centerL', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 17 title 'centerM', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 18 title 'centerN', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 19 title 'centerP', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 20 title 'centerQ', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 21 title 'centerR', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 22 title 'centerS', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 23 title 'centerT', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 24 title 'centerV', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 25 title 'centerW', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 26 title 'centerY'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 27 title 'concen'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 28 title 'length'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 29 title 'statent'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 31 title 'statA', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 32 title 'statC', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 33 title 'statG', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 34 title 'statT'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 35 title 'rrent'\n" << "pause -1\n";
		gp_os << "p '" << (name + ".trace").c_str() << "' u 37 title 'rrAC', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 38 title 'rrAG', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 39 title 'rrAT', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 40 title 'rrCG', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 41 title 'rrCT', ";
		gp_os << "'" << (name + ".trace").c_str() << "' u 42 title 'rrGT'\n" << "pause -1\n";
		*/
	}
};

int main(int argc, char* argv[])	{

	string name = "";
	string datafile = "";
	string treefile = "";
	string profilefile = "none";
	int P = 1;
	GeneticCodeType type = Universal;
	int every = 1;
	int until = -1;

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
			else if (s == "-every")	{
				i++;
				every = atoi(argv[i]);
			}
			else if (s == "-until")	{
				i++;
				until = atoi(argv[i]);
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if (s == "-mtmam")	{
				type = MtMam;
			}
			else if (s == "-mtinv")	{
				type = MtInv;
			}
			else if (s == "-mtprot")	{
				type = MtProt;
			}
			else if (s == "-mtech")	{
				type = MtEch;
			}
			else if (s == "-gencode")	{
				i++;
				s = argv[i];
				if (s == "1")	{	// NCBI standard label
					type = Universal;
				}
				else if (s == "2")	{	// NCBI standard
					type = MtMam;
				}
				else if (s == "3")	{	// NCBI standard
					cout << "genetic code not yet implemented...\n";
					exit(1);
				}
				else if (s == "4")	{	// NCBI standard
					type = MtProt;
				}
				else if (s == "5")	{	// NCBI standard
					type = MtInv;
				}
				else if (s == "6")	{	// NCBI standard
					cout << "genetic code not yet implemented...\n";
					exit(1);
				}
				else if (s == "7")	{	// NCBI standard
					cout << "genetic code deleted...\n";
					exit(1);
				}
				else if (s == "8")	{	// NCBI standard
					cout << "genetic code deleted...\n";
					exit(1);
				}
				else if (s == "9")	{	// NCBI standard
					type = MtEch;
				}
				else	{
					cerr << "error with genetic code choice...\n";
				}

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


	CodonAAProfileMutSelInfMixChain* chain;
	if (datafile == "")	{
		// This is an already existing chain on the disk; reopen and restart
		chain = new CodonAAProfileMutSelInfMixChain(name);
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
	}
	else	{
		// This is a new chain.
		//
		chain = new CodonAAProfileMutSelInfMixChain(datafile,treefile,P,name,type,every,until);
		chain->Gnuplot();
		cerr << "start\n";
		chain->Start();
		cerr << "exit\n";
	}
}

