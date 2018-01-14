
#include "Chain.h"
#include "OmegaSelectionModel.h"
#include <cmath>


class OmegaSelectionChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string contdatafile;
	string treefile;
	int category;

	public:

	OmegaSelectionModel* GetModel() {return (OmegaSelectionModel*) model;}

	string GetModelType() {return modeltype;}

	OmegaSelectionChain(string indata, string incontdata, string intree, int incategory, int inevery, int inuntil, string filename, int force = 0 )	{
		modeltype = "OMEGADIFF";
		datafile = indata;
		contdatafile = incontdata;
		treefile = intree;
		category = incategory;
		every = inevery;
		until = inuntil;
		name = filename;
		New(force);

	}

	OmegaSelectionChain(string filename)	{
		name = filename;
		Open();
		Save();
	}

	void New(int force)	{
		if (modeltype == "OMEGADIFF")	{
			model = new OmegaSelectionModel(datafile,contdatafile,treefile,category,true);
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
		is >> datafile >> contdatafile >> treefile >> category;
		is >> every >> until >> size;

		if (modeltype == "OMEGADIFF")	{
			model = new OmegaSelectionModel(datafile,contdatafile,treefile,category,true);
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
		param_os << datafile << '\t' << contdatafile << '\t' << treefile << '\t' << category << '\n';
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

    OmegaSelectionChain* chain = 0;
	// this is an already existing chain on the disk; reopen and restart
	if (argc == 2)	{
		string name = argv[1];
		chain = new OmegaSelectionChain(name);
	}

	// this is a new chain
	else    {

        string datafile = "";
        string contdatafile = "";
        string treefile = "";
        int ncond = 2;

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

		chain = new OmegaSelectionChain(datafile,contdatafile,treefile,ncond,every,until,name,force);

	}
    cerr << "START\n";
	chain->Start();
}


