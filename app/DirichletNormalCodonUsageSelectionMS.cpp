#include <tclap/CmdLine.h>
#include <cmath>
#include "DirichletNormalCodonUsageSelectionModelMS.hpp"
#include "core/Chain.hpp"
using namespace std;
using namespace TCLAP;

class DirichletNormalCodonUsageSelectionChainMS : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile, codonmodel;
    int category, level, fixglob, fixvar, conjugate;

  public:
    DirichletNormalCodonUsageSelectionModelMS* GetModel() {
        return static_cast<DirichletNormalCodonUsageSelectionModelMS*>(model);
    }

    string GetModelType() override { return modeltype; }

    DirichletNormalCodonUsageSelectionChainMS(string indata, string intree, int incategory,
                                              int inlevel, int inevery, int inuntil, int infixglob,
                                              int infixvar, string incodonmodel, int inconjugate,
                                              string inname, int force)
        : modeltype("DIFFSEL"),
          datafile(indata),
          treefile(intree),
          codonmodel(incodonmodel),
          category(incategory),
          level(inlevel),
          fixglob(infixglob),
          fixvar(infixvar),
          conjugate(inconjugate) {
        SetEvery(inevery);
        SetUntil(inuntil);
        SetName(inname);
        New(force);
    }

    DirichletNormalCodonUsageSelectionChainMS(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DirichletNormalCodonUsageSelectionModelMS(
            datafile, treefile, category, level, fixglob, fixvar, conjugate, codonmodel);
        cerr << "-- Reset" << endl;
        Reset(force);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> category >> level;
        is >> fixglob >> fixvar >> codonmodel;
        is >> conjugate;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DIFFSEL") {
            model = new DirichletNormalCodonUsageSelectionModelMS(
                datafile, treefile, category, level, fixglob, fixvar, conjugate, codonmodel);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        model->FromStream(is);

        model->Update();
        cerr << size << "-- Points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
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
        for (int i = 0; i < every; i++) {
            model->Move(1);
        }

        SavePoint();
        Save();
        Monitor();
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
            CmdLine cmd("To reopen an existing chain, use:\t\t diffsel <name>", ' ', "mini-coevol");

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
                false, 1, "integer", cmd);

            // Performing the actual parsing
            cmd.parse(argc, argv);

            // Extracting the values
            string name = name_arg.getValue();
            int conjugate = not static_cast<int>(conjugate_arg.getValue());
            int fixglob = not static_cast<int>(freeglob_arg.getValue());
            int fixvar = not static_cast<int>(freevar_arg.getValue());
            string codonmodel = squareroot_arg.getValue() ? "SR" : "MS";
            int force = static_cast<int>(force_arg.getValue());

            cerr << "-- Chain name: " << name << endl;

            DirichletNormalCodonUsageSelectionChainMS chain{datafile_arg.getValue(),
                                                            treefile_arg.getValue(),
                                                            ncond_arg.getValue(),
                                                            levels_arg.getValue(),
                                                            every_arg.getValue(),
                                                            until_arg.getValue(),
                                                            fixglob,
                                                            fixvar,
                                                            codonmodel,
                                                            conjugate,
                                                            name,
                                                            force};
            cerr << "-- Starting the chain!" << endl;
            chain.Start();

        } catch (ArgException& e) {
            cerr << "-- Error: " << e.error() << " for arg " << e.argId() << endl;
            exit(1);
        }
    }
}
