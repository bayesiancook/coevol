#define DEBUG 1


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
    int every;
    int burnin;
    string type;
    int conjugate;
    string mechanism;

  public:
    DirichletNormalCodonUsageSelectionModelMS* GetModel() {
        return (DirichletNormalCodonUsageSelectionModelMS*)model;
    }

    string GetModelType() override { return modeltype; }

    DirichletNormalCodonUsageSelectionChainMS(string indata, string intree, int incategory,
                                              int inburnin, int inevery, string filename,
                                              string intype, int inconjugate, string inmechanism,
                                              int force = 0) {
        modeltype = "SELECTIONGTR";
        datafile = indata;
        treefile = intree;
        category = incategory;
        burnin = inburnin;
        every = inevery;
        name = filename;
        type = intype;
        conjugate = inconjugate;
        mechanism = inmechanism;

        New(force);
    }


    DirichletNormalCodonUsageSelectionChainMS(string indata, string intree, int incategory,
                                              int inevery, string filename, string intype,
                                              int inconjugate, string inmechanism, int force = 0) {
        modeltype = "SELECTIONGTR";
        datafile = indata;
        treefile = intree;
        category = incategory;
        every = inevery;
        name = filename;
        type = intype;
        conjugate = inconjugate;
        mechanism = inmechanism;

        New(force);
    }


    DirichletNormalCodonUsageSelectionChainMS(string filename) {
        name = filename;

        conjugate = 1;
        Open();
        Save();
    }

    void New(int force) override {
        if (modeltype == "SELECTIONGTR") {
            model = new DirichletNormalCodonUsageSelectionModelMS(datafile, treefile, category,
                                                                  type, conjugate, mechanism);
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
        is >> datafile >> treefile >> category;
        is >> conjugate;
        is >> type >> mechanism;
        is >> every >> until >> size;

        if (modeltype == "SELECTIONGTR") {
            model = new DirichletNormalCodonUsageSelectionModelMS(datafile, treefile, category,
                                                                  type, conjugate, mechanism);
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
        param_os << datafile << '\t' << treefile << '\t' << category << '\n';
        param_os << conjugate << '\n';
        param_os << type << '\t' << mechanism << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';

        model->ToStream(param_os);
    }


    void Move() override {
#if DEBUG > 0
        MeasureTime myTimer;
#endif
        for (int i = 0; i < every; i++) {
            model->Move(1);
        }

        SavePoint();
        Save();
        Monitor();

#if DEBUG > 0
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
    else if (argc == 9) {
        string datafile = argv[1];
        string treefile = argv[2];
        int category = atoi(argv[3]);
        int every = atoi(argv[4]);
        string name = argv[5];
        string type = argv[6];
        int conjugate = atoi(argv[7]);
        string mechanism = argv[8];

        cerr << "-- Chain name: " << name << endl;

        chain = new DirichletNormalCodonUsageSelectionChainMS(datafile, treefile, category, every,
                                                              name, type, conjugate, mechanism, 1);
    }

    else if (argc == 10) {
        string datafile = argv[1];
        string treefile = argv[2];
        int category = atoi(argv[3]);
        int burnin = atoi(argv[4]);
        int every = atoi(argv[5]);
        string name = argv[6];
        string type = argv[7];
        int conjugate = atoi(argv[8]);
        string mechanism = argv[9];

        cerr << "-- Chain name: " << name << endl;

        chain = new DirichletNormalCodonUsageSelectionChainMS(
            datafile, treefile, category, burnin, every, name, type, conjugate, mechanism, 1);
    } else {
        cerr << "-- Error: incorrect number of command line parameters!" << endl;
        exit(1);
    }
    cerr << "-- Starting the chain!" << endl;
    // chain->SetUntil(20);
    chain->Start();
}
