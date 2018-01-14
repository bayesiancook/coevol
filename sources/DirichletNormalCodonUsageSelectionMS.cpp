
#include "Chain.h"
#include "DirichletNormalCodonUsageSelectionModelMS.h"
#include <cmath>


class DirichletNormalCodonUsageSelectionChainMS : public Chain	{

	private:
	string modeltype;
	string datafile;
    string contdatafile;
	string treefile;
    int ncond;
    int fixglob;
    int fixvar;
    int codonmodel;

	public:

	DirichletNormalCodonUsageSelectionModelMS* GetModel() {return (DirichletNormalCodonUsageSelectionModelMS*) model;}

	string GetModelType() {return modeltype;}

	DirichletNormalCodonUsageSelectionChainMS(string indata, string incontdata, string intree, int inncond, int infixglob, int infixvar, int incodonmodel, int inevery, int inuntil, string filename, int force)	{
		modeltype = "DIFFSELNORM";
		datafile = indata;
        contdatafile = incontdata;
		treefile = intree;
        ncond = inncond;
        fixglob = infixglob;
        fixvar = infixvar;
        codonmodel = incodonmodel;
		every = inevery;
		name = filename;
		New(force);
	}

	DirichletNormalCodonUsageSelectionChainMS(string filename)	{
		name = filename;
		Open();
		Save();
	}

	void New(int force)	{
		if (modeltype == "DIFFSELNORM")	{
			model = new DirichletNormalCodonUsageSelectionModelMS(datafile,contdatafile,treefile,ncond,fixglob,fixvar,codonmodel,true);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "RESET\n";
		Reset(force);
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile >> contdatafile >> treefile;
        is >> ncond;
        is >> fixglob >> fixvar >> codonmodel;
		is >> every >> until >> size;

		if (modeltype == "DIFFSELNORM")	{
			model = new DirichletNormalCodonUsageSelectionModelMS(datafile,contdatafile,treefile,ncond,fixglob,fixvar,codonmodel,true);
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
		param_os << datafile << '\t' << contdatafile << '\t' << treefile << '\n';
        param_os << ncond << '\n';
        param_os << fixglob << '\t' << fixvar << '\t' << codonmodel << '\n';
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

    DirichletNormalCodonUsageSelectionChainMS* chain = 0;
	// this is an already existing chain on the disk; reopen and restart
	if (argc == 2)	{
		string name = argv[1];
		chain = new DirichletNormalCodonUsageSelectionChainMS(name);
	}

	// this is a new chain
	else if (argc == 9)	{

        string datafile = "";
        string contdatafile = "";
        string treefile = "";
        int ncond = 2;
        int fixglob = 1;
        int fixvar = 1;
        int codonmodel = 1;

        string name = "";
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
                else if (s == "-c") {
                    i++;
                    contdatafile = argv[i];
                }
                else if ((s == "-t") || (s == "-T"))	{
                    i++;
                    treefile = argv[i];
                }
                else if (s == "-f") {
                    force = 1;
                }
                else if (s == "-ncond")	{
                    i++;
                    ncond = atoi(argv[i]);
                }
                else if (s == "-ms")    {
                    codonmodel = 1;
                }
                else if (s == "-sr")    {
                    codonmodel = 0;
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
        }
        catch(...)	{
            cerr << "error in command\n";
            exit(1);
        }

        chain = new DirichletNormalCodonUsageSelectionChainMS(datafile,contdatafile,treefile,ncond,fixglob,fixvar,codonmodel,every,until,name,force);

	}
    cerr << "START\n";
	chain->Start();
}

