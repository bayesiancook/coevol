#include <tclap/CmdLine.h>
#include <cmath>
#include "DirichletNormalCodonUsageSelectionModelMS.hpp"
#include "core/Chain.hpp"
using namespace std;
using namespace TCLAP;

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

    DirichletNormalCodonUsageSelectionChainMS(string indata, string intree, int incategory,
                                              int inlevel, int inevery, int inuntil, int infixglob,
                                              int infixvar, string incodonmodel, int inconjugate,
                                              string inname, int force) {
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
            model = new DirichletNormalCodonUsageSelectionModelMS(
                datafile, treefile, category, level, fixglob, fixvar, conjugate, codonmodel);
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
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DIFFSEL") {
            model = new DirichletNormalCodonUsageSelectionModelMS(
                datafile, treefile, category, level, fixglob, fixvar, conjugate, codonmodel);
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
    cerr << "-- Parsing command line arguments\n";

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        cerr << "-- Trying to reopen existing chain named " << name << " on disk\n";
        DirichletNormalCodonUsageSelectionChainMS chain{name};
    }

    // this is a new chain
    else {
        try {
            CmdLine cmd("The diffsel application.", ' ', "mini-coevol");

            // Positional arguments
            UnlabeledValueArg<string> name_arg(
                "name",
                "The name of the run (used to name output file and when resuming computing later).",
                false, "tmp_diffsel", "name", cmd);

            // Switches
            SwitchArg force_arg(
                "f", "force", "Forces diffsel to erase existing chain with same name", cmd, false);
            SwitchArg conjugate_arg("c", "noconj", "Disable conjugate sampling", cmd, false);
            SwitchArg freeglob_arg("g", "freeglob", "Frees glob", cmd, false);
            SwitchArg freevar_arg("v", "freevar", "Frees var", cmd, false);
            SwitchArg squareroot_arg("s", "squareroot", "Uses the square root method.", cmd, false);


            // Value args
            ValueArg<string> datafile_arg("d", "data", "Alignment file in phylip format", true, "",
                                          "filename", cmd);
            ValueArg<string> treefile_arg("t", "tree", "Tree file in newick format", true, "",
                                          "filename", cmd);
            ValueArg<int> ncond_arg("n", "ncond", "Number of conditions", false, 1, "integer", cmd);
            ValueArg<int> levels_arg("l", "levels", "Number of levels", false, 2, "integer", cmd);
            ValueArg<int> until_arg("u", "until", "Number of iteration blocks to perform.", false,
                                    -1, "integer", cmd);
            ValueArg<int> every_arg(
                "e", "every",
                "Number of iteration per iteration block, ie, number of iterations "
                "run before saving data to disk.",
                false, 2, "integer", cmd);

            cmd.parse(argc, argv);

            string datafile = datafile_arg.getValue();
            string treefile = treefile_arg.getValue();
            int category = ncond_arg.getValue();
            int level = levels_arg.getValue();
            int every = every_arg.getValue();
            int until = until_arg.getValue();
            string name = name_arg.getValue();
            int conjugate = not static_cast<int>(conjugate_arg.getValue());
            int fixglob = not static_cast<int>(freeglob_arg.getValue());
            int fixvar = not static_cast<int>(freevar_arg.getValue());
            string codonmodel = squareroot_arg.getValue() ? "SR" : "MS";
            int force = static_cast<int>(force_arg.getValue());

            cerr << "-- Chain name: " << name << endl;

            DirichletNormalCodonUsageSelectionChainMS chain{datafile,   treefile,  category, level,
                                                            every,      until,     fixglob,  fixvar,
                                                            codonmodel, conjugate, name,     force};
            cerr << "-- Starting the chain!" << endl;
            // chain->SetUntil(20);
            chain.Start();


        } catch (ArgException& e) {
            cerr << "-- Error: " << e.error() << " for arg " << e.argId() << endl;
        } catch (...) {
            cerr << "-- Error: incorrect number of command line parameters!" << endl;
            cerr << "diffsel -d datafile ...\n";
            exit(1);
        }
    }
}
