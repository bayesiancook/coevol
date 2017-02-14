#ifndef SELECTIONGTR_H
#define SELECTIONGTR_H

#include <stdio.h>

#include "BranchProcess.hpp"
#include "CodonSequenceAlignment.hpp"
#include "GTRSubMatrix.hpp"
#include "IIDNormalIIDArray.hpp"
#include "MSCodonSubMatrix.hpp"
#include "ProfileConjugatePath.hpp"
#include "SelectionPhyloProcess.hpp"
#include "core/ProbModel.hpp"


class DirichletNormalCompMove : public MCUpdate {
  public:
    DirichletNormalCompMove(DirichletIIDArray* inglobal, IIDNormalIIDArray* indiff, double intuning,
                            int inn, int innrep) {
        global = inglobal;
        diff = indiff;
        tuning = intuning;
        n = inn;
        nrep = innrep;
        mnodearray = new Mnode*[GetSize()];
        for (int i = 0; i < GetSize(); i++) {
            mnodearray[i] = new Mnode;
            global->GetDirichletVal(i)->Register(mnodearray[i]);
            diff->GetNormalVal(i)->Register(mnodearray[i]);
        }
    }

    int GetSize() { return global->GetSize(); }

    double Move(double /*tuning_modulator*/) override {
        int Naccepted = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < GetSize(); i++) {
                Mnode* mnode = mnodearray[i];
                Dirichlet* profile = global->GetDirichletVal(i);
                IIDNormal* normal = diff->GetNormalVal(i);

                mnode->Corrupt(true);

                double logHastings = 0;

                int dim = normal->GetDim();
                if (2 * n > dim) {
                    n = dim / 2;
                }

                auto bkprofile = new double[dim];
                for (int k = 0; k < dim; k++) {
                    bkprofile[k] = (*profile)[k];
                }
                auto indices = new int[2 * n];
                Random::DrawFromUrn(indices, 2 * n, dim);
                for (int i = 0; i < n; i++) {
                    int i1 = indices[2 * i];
                    int i2 = indices[2 * i + 1];
                    double tot = (*profile)[i1] + (*profile)[i2];
                    double x = (*profile)[i1];

                    double h = tot * tuning * (Random::Uniform() - 0.5);
                    x += h;
                    while ((x < 0) || (x > tot)) {
                        if (x < 0) {
                            x = -x;
                        }
                        if (x > tot) {
                            x = 2 * tot - x;
                        }
                    }

                    (*profile)[i1] = x;
                    (*profile)[i2] = tot - x;
                }

                for (int k = 0; k < dim; k++) {
                    (*normal)[k] += log(bkprofile[k] / (*profile)[k]);
                }
                delete[] indices;
                delete[] bkprofile;

                double logratio = mnode->Update() + logHastings;
                bool accepted = (log(Random::Uniform()) < logratio);

                if (!accepted) {
                    mnode->Corrupt(false);
                    mnode->Restore();
                    Naccepted++;
                }
            }
        }
        return (double)Naccepted / GetSize();
    }

  private:
    Mnode** mnodearray;
    DirichletIIDArray* global;
    IIDNormalIIDArray* diff;
    double tuning;
    int n;
    int nrep;
};

class NormalNormalCompMove : public MCUpdate {
  public:
    NormalNormalCompMove(IIDNormalIIDArray* indiff1, IIDNormalIIDArray* indiff2, double intuning,
                         int inn, int innrep) {
        diff1 = indiff1;
        diff2 = indiff2;
        tuning = intuning;
        n = inn;
        nrep = innrep;
        mnodearray = new Mnode*[GetSize()];
        for (int i = 0; i < GetSize(); i++) {
            mnodearray[i] = new Mnode;
            diff1->GetNormalVal(i)->Register(mnodearray[i]);
            diff2->GetNormalVal(i)->Register(mnodearray[i]);
        }
    }

    int GetSize() { return diff1->GetSize(); }

    double Move(double /*tuning_modulator*/) override {
        int Naccepted = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < GetSize(); i++) {
                Mnode* mnode = mnodearray[i];
                IIDNormal* normal1 = diff1->GetNormalVal(i);
                IIDNormal* normal2 = diff2->GetNormalVal(i);

                mnode->Corrupt(true);

                double logHastings = 0;

                int dim = normal1->GetDim();
                if (n > dim) {
                    n = dim;
                }

                auto indices = new int[n];
                Random::DrawFromUrn(indices, n, dim);
                for (int i = 0; i < n; i++) {
                    int j = indices[i];
                    double h = tuning * (Random::Uniform() - 0.5);
                    (*normal1)[j] += h;
                    (*normal2)[j] -= h;
                }
                delete[] indices;

                double logratio = mnode->Update() + logHastings;
                bool accepted = (log(Random::Uniform()) < logratio);

                if (!accepted) {
                    mnode->Corrupt(false);
                    mnode->Restore();
                    Naccepted++;
                }
            }
        }
        return (double)Naccepted / GetSize();
    }

  private:
    Mnode** mnodearray;
    IIDNormalIIDArray* diff1;
    IIDNormalIIDArray* diff2;
    double tuning;
    int n;
    int nrep;
};

class ComplexDirichletIIDArrayMove : public MCUpdate {
  public:
    ComplexDirichletIIDArrayMove(DirichletIIDArray* inselectarray, double intuning, int innrep)
        : selectarray(inselectarray), tuning(intuning), nrep(innrep) {}

    double Move(double /*tuning_modulator*/) override {
        auto tot = new double[selectarray->GetSize()];
        for (int i = 0; i < selectarray->GetSize(); i++) {
            tot[i] = 0;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif

        for (int i = 0; i < selectarray->GetSize(); i++) {
            for (int rep = 0; rep < nrep; rep++) {
                tot[i] += selectarray->GetDirichletVal(i)->Move(tuning * 1, 2);
                tot[i] += selectarray->GetDirichletVal(i)->Move(tuning * 0.3, 5);
                tot[i] += selectarray->GetDirichletVal(i)->Move(tuning * 0.1, 10);
                // tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.3,2);
                // tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*1,5);
                // tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.03,10);
                // tot[i] += selectarray->GetDirichletVal(i)->Move(tuning*0.03,10);
            }
        }

        double total = 0;
        for (int i = 0; i < selectarray->GetSize(); i++) {
            total += tot[i];
        }
        delete[] tot;

        return total / selectarray->GetSize() / nrep / 6;
    }

  private:
    DirichletIIDArray* selectarray;
    double tuning;
    int nrep;
};

class RenormalizedIIDStat : public Dvar<Profile> {
  public:
    // the constructor takes 4 arguments: pointers to the 3 variables a and b and c and the variable
    // beta
    // which we want to multiply and renormalized, all this in a componentwise fashion.
    //
    RenormalizedIIDStat(Var<Profile>* ina, Var<RealVector>* inb, Var<RealVector>* inc,
                        Var<PosReal>* inbeta) {
        if (ina->GetDim() != inb->GetDim()) {
            cerr << "error in RenormalizedIIDStat : non matching dimension (" << ina->GetDim()
                 << " vs " << inb->GetDim() << ")\n";
            throw;
        }
        setval(Profile(ina->GetDim()));
        bkvalue = Profile(ina->GetDim());
        a = ina;
        b = inb;
        c = inc;
        beta = inbeta;
        Register(a);
        Register(b);
        Register(beta);
        Register(c);
        specialUpdate();
    }

    // the only function that we need to define is the one that makes multiplication and
    // renormalization
    void specialUpdate() override {
        // divide
        double total = 0;

        for (int i = 0; i < GetDim(); i++) {
            (*this)[i] = exp(beta->val() * (log((*a)[i]) + (*b)[i] + (*c)[i]));
            total += (*this)[i];
        }

        // renormalize
        for (int i = 0; i < GetDim(); i++) {
            (*this)[i] /= total;
        }
    }

  private:
    Var<Profile>* a;
    Var<RealVector>* b;
    Var<PosReal>* beta;
    Var<RealVector>* c;
};

/*
  class SumConstrainedRealVector : public Dvar<RealVector>	{

  public:

  SumConstrainedRealVector(Var<RealVector>* insource)	{
  setval(RealVector(insource->GetDim()));
  bkvalue = RealVector(insource->GetDim());
  source = insource;
  Register(source);
  specialUpdate();
  }

  ~SumConstrainedRealVector()	{
  }

  protected:

  void specialUpdate()	{
  double total = 0;
  for (int i=0; i<GetDim(); i++)	{
  total += (*source)[i];
  }
  total /= GetDim();
  for (int i=0; i<GetDim(); i++)	{
  (*this)[i] = (*source)[i] - total;
  }
  }

  Var<RealVector>* source;
  };
*/


class DirichletNormalCodonUsageSelectionModelMS : public ProbModel {
    // data fields

    // ---------
    // the fixed parameters of the model
    // ---------

    // a fixed tree (read from file)
    Tree* tree;
    FileSequenceAlignment* data;
    TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    // number of sites
    int Nsite;

    // number of states (4 for nucleic acids, 20 for amino-acids, 61 for codons)
    int Nstate;

    // number of categories
    int K;

    // ---------
    // the random variables of the model
    // ---------

    Const<Real>* Zero;
    Const<PosReal>* One;
    Const<PosReal>* Ten;
    Const<PosReal>* beta;

    Exponential* lambda;
    GammaTree* gamtree;

    Dirichlet* relrate;
    Dirichlet* nucstationary;
    GTRRandomSubMatrixWithNormRates* nucmatrix;

    Const<RealVector>* zero;

    Dirichlet* center;
    Const<PosReal>* PriorConcentration;
    Exponential* concentration;
    DirichletIIDArray* globalselectionprofile;

    Exponential** var;
    IIDNormalIIDArray** selectionnormal;
    // SumConstrainedRealVector*** selectionprofile;

    RenormalizedIIDStat*** statarray;
    RandomSubMatrix*** submatrix;

    AllocationTree* allocatetree;

    ProfilePathConjugateArray* patharray;
    SelectionMatrixTree* matrixtree;
    PhyloProcess* phyloprocess;

    double** selectionmean;
    double accrate;

    Dirichlet* codonusageselection;

    int conjugate;
    string type;
    string mechanism;

  public:
    // constructor
    // this is where the entire graph structure of the model is created

    DirichletNormalCodonUsageSelectionModelMS(string datafile, string treefile, int inK,
                                              string intype, int inconjugate, string inmechanism,
                                              bool sample = true) {
        conjugate = inconjugate;
        type = intype;
        mechanism = inmechanism;
        K = inK;


        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();    // # columns
        Nstate = codondata->GetNstate();  // # states (61 for codons)

        int Nnuc = 4;
        int Naa = 20;

        cerr << Nsite << '\t' << Nstate << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        cerr << "tree and data ok\n";

        // ----------
        // construction of the graph
        // ----------

        Zero = new Const<Real>(0);
        One = new Const<PosReal>(1);
        Ten = new Const<PosReal>(10);
        beta = new Const<PosReal>(1);

        // a tree with all branch lengths iid from an exponential distribution of rate lambda
        // this is a gamma distribution of shape mu=1 and scale lambda
        // lambda is itself endowed with an exponential prior of mean 10
        lambda = new Exponential(Ten, Exponential::MEAN);
        lambda->setval(50);
        gamtree = new GammaTree(tree, One, lambda);

        allocatetree = new AllocationTree(tree, K);

        relrate = new Dirichlet(Nnuc * (Nnuc - 1) / 2);
        nucstationary = new Dirichlet(Nnuc);
        nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate, nucstationary, true);
        // true because we want it to be normalised
        // here the branch length is equal to the number of proposed mutations for that branch

        codonusageselection = new Dirichlet(Nstate);
        codonusageselection->setuniform();
        codonusageselection->Clamp();

        // prepare a temporary real vector with entries all = 1/20
        Profile tmp(Naa);
        for (int k = 0; k < Naa; k++) {
            tmp[k] = 1.0 / Naa;
        }

        zero = new Const<RealVector>(Naa);
        for (int i = 0; i < Naa; i++) {
            (*zero)[i] = 0;
        }

        PriorConcentration = new Const<PosReal>(Naa);
        concentration = new Exponential(PriorConcentration, Exponential::MEAN);
        center = new Dirichlet(Naa);

        center->setuniform();
        concentration->setval(Naa);

        /*
          double min = 0.001;
          double max = 1000;
        */

        var = new Exponential*[K];
        for (int k = 1; k < K; k++) {
            var[k] = new Exponential(One, Exponential::MEAN);
            var[k]->setval(1.0);
        }


        if (type == "clamp_MCMC") {
            center->Clamp();
            concentration->Clamp();
        }

        if (type == "clamp_MCMC_var") {
            center->Clamp();
            concentration->Clamp();
            for (int k = 1; k < K; k++) {
                var[k]->Clamp();
            }
        }


        globalselectionprofile = new DirichletIIDArray(Nsite, center, concentration);

        selectionnormal = new IIDNormalIIDArray*[K];
        // selectionprofile = new SumConstrainedRealVector**[K];
        statarray = new RenormalizedIIDStat**[K];


        cerr << "selection profiles\n";
        for (int k = 1; k < K; k++) {
            selectionnormal[k] = new IIDNormalIIDArray(Nsite, Naa, Zero, var[k]);
        }

        /*
          cerr << "stat arrays\n";
          for (int k=1; k<K; k++)	{
          selectionprofile[k] = new SumConstrainedRealVector*[Nsite];
          for (int i=0; i<Nsite;i++){
          selectionprofile[k][i] = new SumConstrainedRealVector(selectionnormal[k]->GetVal(i));
          }
          }
        */

        for (int k = 0; k < K; k++) {
            statarray[k] = new RenormalizedIIDStat*[Nsite];
            for (int i = 0; i < Nsite; i++) {
                if (K == 1) {
                    statarray[k][i] = new RenormalizedIIDStat(globalselectionprofile->GetVal(i),
                                                              zero, zero, beta);
                } else {
                    if (k == 0) {
                        statarray[k][i] = new RenormalizedIIDStat(globalselectionprofile->GetVal(i),
                                                                  zero, zero, beta);
                    }

                    else if (k == 1) {
                        statarray[k][i] =
                            new RenormalizedIIDStat(globalselectionprofile->GetVal(i),
                                                    selectionnormal[1]->GetVal(i), zero, beta);
                    }

                    else if (k == 2) {
                        statarray[k][i] = new RenormalizedIIDStat(
                            globalselectionprofile->GetVal(i), selectionnormal[1]->GetVal(i),
                            selectionnormal[2]->GetVal(i), beta);
                    }

                    else if (k == 3) {
                        statarray[k][i] = new RenormalizedIIDStat(
                            globalselectionprofile->GetVal(i), selectionnormal[1]->GetVal(i),
                            selectionnormal[3]->GetVal(i), beta);
                    }
                }
            }
        }

        cerr << "submatrices\n";

        // Square Root //phenimenological
        if (mechanism == "SR") {
            submatrix = new RandomSubMatrix**[K];
            for (int k = 0; k < K; k++) {
                submatrix[k] = new RandomSubMatrix*[Nsite];
                for (int i = 0; i < Nsite; i++) {
                    submatrix[k][i] = new RandomMGSRFitnessCodonUsageSubMatrix(
                        (CodonStateSpace*)codondata->GetStateSpace(), nucmatrix, statarray[k][i],
                        codonusageselection, false);
                }
            }
        }

        // Mutation Selection // mechanistic
        if (mechanism == "MS") {
            submatrix = new RandomSubMatrix**[K];
            for (int k = 0; k < K; k++) {
                submatrix[k] = new RandomSubMatrix*[Nsite];
                for (int i = 0; i < Nsite; i++) {
                    submatrix[k][i] = new RandomMGMSFitnessCodonUsageSubMatrix(
                        (CodonStateSpace*)codondata->GetStateSpace(), nucmatrix, statarray[k][i],
                        codonusageselection, false);
                }
            }
        }


        if (conjugate != 0) {
            cerr << "conjugate\n";
            matrixtree = nullptr;
            patharray = new ProfilePathConjugateArray(Nsite, K, submatrix);
            patharray->InactivateSufficientStatistic();
            phyloprocess = new ProfileConjugateSelectionPhyloProcess(allocatetree, gamtree,
                                                                     codondata, patharray);
        } else {
            matrixtree = new SelectionMatrixTree(allocatetree, submatrix);
            cerr << "create phylo process\n";
            phyloprocess = new SelectionPhyloProcess(gamtree, nullptr, matrixtree, codondata);
        }

        cerr << "unfold\n";
        phyloprocess->Unfold();

        cerr << "register\n";
        RootRegister(Zero);
        RootRegister(One);
        RootRegister(Ten);
        RootRegister(beta);
        RootRegister(relrate);
        RootRegister(nucstationary);
        RootRegister(zero);
        RootRegister(codonusageselection);
        RootRegister(PriorConcentration);
        RootRegister(center);

        Register();

        if (sample) {
            cerr << "initialise\n";
            // Sample();
            cerr << "sample completed\n";
            Update();
            cerr << "update completed\n";
            cerr << "ln L " << GetLogLikelihood() << '\n';
            // cerr << "random calls " << Random::GetCount() << '\n';
            // Trace(cerr);
            Trace(cerr);
        }

        cerr << "trace completed\n";

        cerr << "scheduler\n";
        MakeScheduler();

        cerr << "model created\n";
    }

    // destructor
    // deallocations should normally be done here
    // but in general, the model is deleted just before the program exits, so we can dispense with
    // it for the moment

    ~DirichletNormalCodonUsageSelectionModelMS() override = default;

    Tree* GetTree() { return tree; }


    double Update() {
        cerr << "update with phyloprocess\n";
        double ret = ProbModel::Update();
        phyloprocess->Move(1);
        ret = ProbModel::Update();
        return ret;
    }

    // log probability of the model is the sum of the log prior and the log likelihood

    double GetLogProb() override { return GetLogPrior() + GetLogLikelihood(); }

    double GetLogPrior() {
        double total = 0;
        total += lambda->GetLogProb();
        total += gamtree->GetLogProb();
        total += relrate->GetLogProb();
        total += nucstationary->GetLogProb();

        total += center->GetLogProb();
        total += concentration->GetLogProb();
        for (int j = 0; j < Nsite; j++) {
            total += globalselectionprofile->GetVal(j)->GetLogProb();
        }

        for (int i = 1; i < K; i++) {
            total += var[i]->GetLogProb();
            for (int j = 0; j < Nsite; j++) {
                total += selectionnormal[i]->GetVal(j)->GetLogProb();
            }
        }

        return total;
    }


    double GetLogLikelihood() {
        // return 0;
        return phyloprocess->GetLogProb();
    }


    void MakeScheduler() override {
        scheduler.Register(new SimpleMove(phyloprocess, 1), 1, "mapping");

        for (int n = 0; n < 3; n++) {
            // first register the moves on the global variables
            scheduler.Register(new SimpleMove(lambda, 1), 10, "lambda");
            scheduler.Register(new SimpleMove(lambda, 0.1), 10, "lambda");

            scheduler.Register(new SimpleMove(gamtree, 3), 4, "branch lengths");
            scheduler.Register(new SimpleMove(gamtree, 2.3), 4, "branch lengths");

            if (conjugate != 0) {
                scheduler.Register(new ProfileConjugateMove(patharray, true), 1,
                                   "activate suff stat");
            }

            for (int m = 0; m < 20; m++) {
                scheduler.Register(new ProfileMove(relrate, 0.1, 1), 1, "relrates");
                scheduler.Register(new ProfileMove(relrate, 0.03, 3), 1, "relrates");
                scheduler.Register(new ProfileMove(relrate, 0.01, 3), 1, "relrates");

                scheduler.Register(new ProfileMove(nucstationary, 0.1, 1), 1, "nucstationary");
                scheduler.Register(new ProfileMove(nucstationary, 0.01, 1), 1, "nucstationary");

                scheduler.Register(new ProfileMove(codonusageselection, 0.3, 1), 1,
                                   "codonusageselection");
                scheduler.Register(new ProfileMove(codonusageselection, 0.1, 3), 1,
                                   "codonusageselection");
                scheduler.Register(new SimpleMove(codonusageselection, 0.01), 1,
                                   "codonusageselection");

                scheduler.Register(
                    new ComplexDirichletIIDArrayMove(globalselectionprofile, 0.15, 10), 1,
                    "global stat");

                /*
                  scheduler.Register(new
                  DirichletNormalCompMove(globalselectionprofile,selectionnormal[1],0.03,1,5),1,"global
                  diff comp");
                  scheduler.Register(new
                  DirichletNormalCompMove(globalselectionprofile,selectionnormal[1],0.01,3,5),1,"global
                  diff comp");
                  scheduler.Register(new
                  DirichletNormalCompMove(globalselectionprofile,selectionnormal[1],0.001,10,5),1,"global
                  diff comp");
                */

                scheduler.Register(new ProfileMove(center, 0.1, 1), 10, "center");
                scheduler.Register(new ProfileMove(center, 0.01, 5), 10, "center");

                scheduler.Register(new SimpleMove(concentration, 0.1), 10, "conc");
                scheduler.Register(new SimpleMove(concentration, 0.03), 10, "conc");

                for (int i = 1; i < K; i++) {
                    scheduler.Register(new SimpleMove(var[i], 1), 10, "var");
                    scheduler.Register(new SimpleMove(var[i], 0.1), 10, "var");

                    std::stringstream temp;
                    std::string indice;
                    temp << i;
                    temp >> indice;

                    scheduler.Register(new IIDNormalIIDArrayMove(selectionnormal[i], 5, 1), 10,
                                       "stat" + indice);
                    scheduler.Register(new IIDNormalIIDArrayMove(selectionnormal[i], 3, 5), 10,
                                       "stat" + indice);
                    scheduler.Register(new IIDNormalIIDArrayMove(selectionnormal[i], 1, 10), 10,
                                       "stat" + indice);
                    scheduler.Register(new IIDNormalIIDArrayMove(selectionnormal[i], 1, 20), 10,
                                       "stat" + indice);
                }

                /*
                  for(int i=2;i<K;i++){
                  scheduler.Register(new
                  NormalNormalCompMove(selectionnormal[1],selectionnormal[i],1,10,5),1,"diff diff
                  comp 1/10");
                  scheduler.Register(new
                  NormalNormalCompMove(selectionnormal[1],selectionnormal[i],3,10,5),1,"diff diff
                  comp 3/10");
                  scheduler.Register(new
                  NormalNormalCompMove(selectionnormal[1],selectionnormal[i],5,3,5),1,"diff diff
                  comp 5/3");
                  }
                */
            }

            if (conjugate != 0) {
                scheduler.Register(new ProfileConjugateMove(patharray, false), 1,
                                   "inactivate suff stat");
            }
        }
    }


    double Move(double /*tuning_modulator*/) override {
        scheduler.Cycle(1, 1, false, false);
        return 1;
    }


    void drawSample() override {
        cerr << "in sample\n";
        exit(1);
    }

    int GetSite() { return Nsite; }

    int GetaaState() { return Naa; }

    int GetCodonState() { return Nstate; }

    int GetCategory() { return K; }

    double GetBeta() { return beta->val(); }

    void OutputSelectionProfile(string basename, int sitemin, int sitemax) {
        if (sitemin == -1) {
            sitemin = 0;
            sitemax = Nsite;
        }
        ofstream os((basename + ".trueprofiles").c_str());
        for (int k = 1; k < K; k++) {
            os << k << '\n';
            for (int i = sitemin; i < sitemax; i++) {
                os << i << '\t';
                for (int j = 0; j < Naa; j++) {
                    os << (*selectionnormal[k]->GetVal(i))[j] << '\t';
                }
                os << '\n';
            }
            os << '\n';
        }
    }

    void GetSelectionProfile(string basename) {
        ifstream is((basename + ".trueprofiles").c_str());

        for (int k = 1; k < K; k++) {
            int tmp;
            is >> tmp;
            if (tmp != k) {
                cerr << "error when reading true selection profiles\n";
                exit(1);
            }
            for (int i = 0; i < Nsite; i++) {
                int tmp;
                is >> tmp;
                /*
                  if (tmp != i)	{
                  cerr << "error when reading true selection profiles\n";
                  exit(1);
                  }
                */

                for (int j = 0; j < Naa; j++) {
                    is >> (*selectionnormal[k]->GetVal(i))[j];
                }
            }
        }
    }

    double GetSelectionProfile(int cat, int site, int state) {
        int k = cat;
        int i = site;
        int j = state;
        if (cat != 0) {
            return exp((*selectionnormal[k]->GetVal(i))[j]);
        }
        return (*globalselectionprofile->GetVal(i))[j];
    }

    double GetLength() { return gamtree->GetTotalLength(); }

    double GetCenterEntropy() {
        double tot = 0;
        for (int i = 0; i < center->GetDim(); i++) {
            double tmp = (*center)[i];
            tot -= tmp * log(tmp);
        }
        return tot;
    }


    double GetMeanVar(int k) {
        double totvar = 0;
        for (int i = 0; i < Nsite; i++) {
            double var = 0;
            for (int j = 0; j < Naa; j++) {
                double tmp = (*selectionnormal[k]->GetVal(i))[j];
                var += tmp * tmp;
            }
            var /= Naa;
            totvar += var;
        }
        return totvar / Nsite;
    }

    double GetStationary(int cat, int site, int state) {
        int k = cat;
        int i = site;
        int j = state;

        return (*statarray[k][i])[j];
    }

    GTRRandomSubMatrixWithNormRates* Getnucmatrix() { return nucmatrix; }

    RandomSubMatrix* Getsubmatrix(int cat, int site) {
        int k = cat;
        int i = site;
        return (submatrix[k][i]);
    }

    double GetCodonUsage(int state) {
        int i = state;
        return (*codonusageselection)[i];
    }


    CodonStateSpace* GetCodonStateSpace() { return codondata->GetCodonStateSpace(); }

    // creates the header of the <model_name>.trace file
    void TraceHeader(ostream& os) override {
        os << "#logprior\tlnL\tlength\t";
        os << "globent\tcenter\tconc\t";
        for (int i = 1; i < K; i++) {
            os << "selvar" << i << '\t';
            os << "var" << i << '\t';
        }
        os << "rrent\n";
    }


    // writes all summary statistics on one single line
    // in the same order as that provided by the header
    void Trace(ostream& os) override {
        os << GetLogPrior();
        os << '\t' << GetLogLikelihood();
        os << '\t' << GetLength();
        os << '\t' << globalselectionprofile->GetMeanEntropy();
        os << '\t' << GetCenterEntropy();
        os << '\t' << concentration->val();
        for (int i = 1; i < K; i++) {
            os << '\t' << GetMeanVar(i);
            os << '\t' << var[i]->val();
        }
        os << '\t' << relrate->val().GetEntropy();
        os << '\n';
        os.flush();
    }

    void ToStream(ostream& os) override {
        os.precision(7);
        os << *lambda << '\n';
        os << *gamtree << '\n';
        os << *relrate << '\n';
        os << *nucstationary << '\n';
        os << *codonusageselection << '\n';
        os << *center << '\n';
        os << *concentration << '\n';
        os << *globalselectionprofile << '\n';
        for (int i = 1; i < K; i++) {
            os << *var[i] << '\n';
            os << *selectionnormal[i] << '\n';
        }
    }

    void FromStream(istream& is) override {
        is >> *lambda;
        is >> *gamtree;
        is >> *relrate;
        is >> *nucstationary;
        is >> *codonusageselection;
        is >> *center;
        is >> *concentration;
        is >> *globalselectionprofile;
        for (int i = 1; i < K; i++) {
            is >> *var[i];
            // cerr << "var : " << *var[i]<< '\n';
            is >> *selectionnormal[i];
            // cerr << " selectionnormal : " << *selectionnormal[1]<< '\n';
        }
    }

    void PostPredAli(string name, int sitemin, int sitemax) {
        auto datacopy = new CodonSequenceAlignment(codondata);
        phyloprocess->PostPredSample(false);
        phyloprocess->GetLeafData(datacopy);
        ostringstream s;
        ofstream os((name + ".ali").c_str());
        if (sitemax != -1) {
            auto datacopy2 = new CodonSequenceAlignment(datacopy, sitemin, sitemax);
            datacopy2->ToStream(os);
        } else {
            datacopy->ToStream(os);
        }
    }
};

#endif
