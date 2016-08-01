
#include "Chain.h"
// #include "DPFlexDiscreteNHLGAminoAcidModel.h.bk2"
// #include "DPFlexDiscreteNHLGAminoAcidModel.h"
// #include "DeltaNormDPFlexDiscreteNHLGAminoAcidModel.h"
#include "DirDPFlexDiscreteNHLGAminoAcidModel.h"
// #include "FlexDiscreteNHLGAminoAcidModel.h"

class FlexDiscreteNHAminoAcidChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	bool withsep;
	GeneticCodeType type;

	public:

	ProbModel* GetModel() {return (ProbModel*) model;}

	string GetModelType() {return modeltype;}

	FlexDiscreteNHAminoAcidChain(string indata, string intree, string incontdata, string filename, bool inwithsep, GeneticCodeType intype, int force = 1)	{
		modeltype = "FLEXDISCRETENHAMINOACID";
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		withsep = inwithsep;
		type = intype;
		name = filename;
		New(force);
	}

	FlexDiscreteNHAminoAcidChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "FLEXDISCRETENHAMINOACID")	{
			model = new FlexDiscreteNHAminoAcidModel(datafile,treefile,contdatafile,withsep,true,type);
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
		is >> datafile >> treefile >> contdatafile;
		is >> withsep;
		is >> every >> until >> size;

		if (modeltype == "FLEXDISCRETENHAMINOACID")	{
			model = new FlexDiscreteNHAminoAcidModel(datafile,treefile,contdatafile,withsep,false,type);
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
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		param_os << withsep << '\n';
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
	string contdatafile = "None";
	GeneticCodeType type = Universal;

	bool withsep = false;
	
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
			else if (s == "-c")	{
				i++;
				contdatafile = argv[i];
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if (s == "-mtmam")	{
				type = MtMam;
			}
			else if (s == "-sep")	{
				withsep = true;
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

	FlexDiscreteNHAminoAcidChain* chain;
	if (datafile == "")	{
		chain = new FlexDiscreteNHAminoAcidChain(name);
	}
	else	{
		chain = new FlexDiscreteNHAminoAcidChain(datafile,treefile,contdatafile,name,withsep,type);
	}

	cerr << "start\n";
	chain->Start();
	cerr << "exit\n";
}

