#include <cmath>
#include "DirichletNormalCodonUsageSelectionModelMS.hpp"
#include "core/Chain.hpp"
using namespace std;

class DirichletNormalCodonUsageSelectionChainMS : public Chain {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int category;
    int level;
    int fixglob;
    int fixvar;
    string codonmodel;
    int conjugate;

  public:
    DirichletNormalCodonUsageSelectionModelMS* GetModel() {
        return (DirichletNormalCodonUsageSelectionModelMS*)model;
    }

    string GetModelType() override { return modeltype; }

    DirichletNormalCodonUsageSelectionChainMS(string indata, string intree, int incategory, int inlevel,
                                              int inevery, int inuntil,
					      int infixglob, int infixvar, string incodonmodel, int inconjugate,
					      string inname,
                                              int force) {
        modeltype = "DIFFSEL";
        datafile = indata;
        treefile = intree;
        category = incategory;
	level = inlevel;
        every = inevery;
	until = inuntil;
	fixglob = infixglob;
	fixvar = infixvar;
	codonmodel = incodonmodel;
        conjugate = inconjugate;
        name = inname;

        New(force);
    }

    DirichletNormalCodonUsageSelectionChainMS(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        if (modeltype == "DIFFSEL") {
            model = new DirichletNormalCodonUsageSelectionModelMS(datafile, treefile, category, level,
                                                                  fixglob, fixvar, conjugate, codonmodel);
        } else {
            cerr << "error, does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        cerr << "-- Reset" << endl;
        Reset(force);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> category >> level;
        is >> fixglob >> fixvar >> codonmodel;
        is >> conjugate;
	int tmp;
	is >> tmp;
	if (tmp)	{
		cerr << "error when reading model\n";
		exit(1);
	}
        is >> every >> until >> size;

        if (modeltype == "DIFFSEL") {
            model = new DirichletNormalCodonUsageSelectionModelMS(datafile, treefile, category, level,
                                                                  fixglob, fixvar, conjugate, codonmodel);
        } else {
            cerr << "error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        model->FromStream(is);

        model->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << category << '\t' << level << '\n';
	param_os << fixglob << '\t' << fixvar << '\t' << codonmodel << '\n';
        param_os << conjugate << '\n';
	param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';

        model->ToStream(param_os);
    }


    void Move() override {
#if DEBUG > 1
        MeasureTime myTimer;
#endif
        for (int i = 0; i < every; i++) {
            model->Move(1);
        }

        SavePoint();
        Save();
        Monitor();

#if DEBUG > 1
        myTimer << "Performed " << every << " iterations. ";
        myTimer.print<0>();
#endif
    }
};


int main(int argc, char* argv[]) {
    // this is an already existing chain on the disk; reopen and restart
    DirichletNormalCodonUsageSelectionChainMS* chain;

    if (argc == 2) {
        string name = argv[1];
        chain = new DirichletNormalCodonUsageSelectionChainMS(name);
    }

    // this is a new chain
    else	{

        string datafile = "";
        string treefile = "";
        int category = 1;
	int level = 2;
        int every = 1;
	int until = -1;
        string name = "";
	int conjugate = 1;
	int fixglob = 1;
	int fixvar = 1;
	string codonmodel = "MS";
	int force = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}
		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-f")	{
				force = 1;
			}
			else if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-t")	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-ncond")	{
				i++;
				category = atoi(argv[i]);
			}
			else if (s == "-levels")	{
				i++;
				level = atoi(argv[i]);
			}
			else if (s == "-conj")	{
				conjugate = 1;
			}
			else if (s == "-nonconj")	{
				conjugate = 0;
			}
			else if (s == "-fixglob")	{
				fixglob = 1;
			}
			else if (s == "-freeglob")	{
				fixglob = 0;
			}
			else if (s == "-fixvar")	{
				fixvar = 1;
			}
			else if (s == "-freevar")	{
				fixvar = 0;
			}
			else if ((s == "-ms") || (s == "-mutsel"))	{
				codonmodel = "MS";
			}
			else if ((s == "-sr")  || (s == "-squareroot"))	{
				codonmodel = "SR";
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
	}
	catch(...)	{
		cerr << "-- Error: incorrect number of command line parameters!" << endl;
		cerr << "diffsel -d datafile ...\n";
		exit(1);
	}
        cerr << "-- Chain name: " << name << endl;

        chain = new DirichletNormalCodonUsageSelectionChainMS(datafile, treefile, category, level, every, until, fixglob, fixvar, codonmodel, conjugate, name, force);
    }
    cerr << "-- Starting the chain!" << endl;
    // chain->SetUntil(20);
    chain->Start();
}
