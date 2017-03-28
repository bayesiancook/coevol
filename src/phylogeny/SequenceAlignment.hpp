#ifndef SEQUENCEALIGNMENT_H
#define SEQUENCEALIGNMENT_H

#include <vector>
#include "StateSpace.hpp"
#include "TaxonSet.hpp"

// this class works like an interface
// it does not do any job
class SequenceAlignment {
  public:
    SequenceAlignment() : Data(nullptr) {}

    SequenceAlignment(SequenceAlignment *from) {
        Ntaxa = from->Ntaxa;
        Nsite = from->Nsite;
        taxset = from->taxset;
        statespace = from->statespace;

        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = from->Data[i][j];
            }
        }
    }

    SequenceAlignment(SequenceAlignment *from, int start, int length) {
        Ntaxa = from->Ntaxa;
        if ((start < 0) || (length < 0) || (start + length > from->Nsite)) {
            std::cerr << "error in sequence alignment: overflow\n";
            exit(1);
        }

        Nsite = length;
        taxset = from->taxset;
        statespace = from->statespace;
        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = from->Data[i][j + start];
            }
        }
    }

    void SubSelect(int sitemin, int sitemax) {
        auto Data2 = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data2[i] = new int[sitemax - sitemin];
            for (int j = sitemin; j < sitemax; j++) {
                Data2[i][j - sitemin] = Data[i][j];
            }
        }
        for (int i = 0; i < Ntaxa; i++) {
            delete[] Data[i];
        }
        delete[] Data;
        Data = Data2;
        Nsite = sitemax - sitemin;
    }

    SequenceAlignment(SequenceAlignment *from, TaxonSet *subset) {
        Ntaxa = subset->GetNtaxa();
        Nsite = from->GetNsite();
        statespace = from->statespace;
        taxset = subset;
        Data = new int *[Ntaxa];
        for (int k = 0; k < Ntaxa; k++) {
            Data[k] = new int[Nsite];
        }
        for (int i = 0; i < from->GetNtaxa(); i++) {
            int k = subset->GetTaxonIndex(from->GetTaxonSet()->GetTaxon(i));
            if (k != -1) {
                for (int j = 0; j < Nsite; j++) {
                    Data[k][j] = from->Data[i][j];
                }
                if (k >= Ntaxa) {
                    std::cerr << "error in sequence alignment subset\n";
                    std::cerr << Ntaxa << '\n';
                    exit(1);
                }
            }
        }
    }

    SequenceAlignment(SequenceAlignment *from, std::map<std::string, int> &group, int ngroup) {
        int n = 0;
        for (int i = 0; i < from->GetNsite(); i++) {
            if (from->GroupCompatible(i, group, ngroup)) {
                n++;
            }
        }
        statespace = from->statespace;
        taxset = from->taxset;
        Ntaxa = from->GetNtaxa();
        Nsite = n;
        Data = new int *[Ntaxa];
        for (int k = 0; k < Ntaxa; k++) {
            Data[k] = new int[Nsite];
        }
        int j = 0;
        for (int i = 0; i < from->GetNsite(); i++) {
            if (from->GroupCompatible(i, group, ngroup)) {
                for (int k = 0; k < Ntaxa; k++) {
                    Data[k][j] = from->Data[k][i];
                }
                j++;
            }
        }
    }

    void MaskOutgroup(std::string *group1, std::string *group2, int n1, int n2) {
        int nrem = 0;
        int npos = 0;
        for (int i = 0; i < GetNsite(); i++) {
            int ok = 0;
            for (int j = 0; j < n1; j++) {
                if (Data[GetTaxonSet()->GetTaxonIndex(group1[j])][i] != unknown) {
                    ok = 1;
                }
            }
            if (ok == 0) {
                npos++;
                for (int j = 0; j < n2; j++) {
                    if (Data[GetTaxonSet()->GetTaxonIndex(group2[j])][i] != unknown) {
                        Data[GetTaxonSet()->GetTaxonIndex(group2[j])][i] = unknown;
                        nrem++;
                    }
                }
            }
        }
        std::cerr << "number of positions masked   : " << npos << '\n';
        std::cerr << "total number of cells masked : " << nrem << '\n';
    }

    void MaskDiversity(std::string *group, int n) {
        int nrem = 0;
        int npos = 0;
        std::vector<int> present(GetStateSpace()->GetNstate());
        for (int i = 0; i < GetNsite(); i++) {
            for (int k = 0; k < GetStateSpace()->GetNstate(); k++) {
                present[k] = 0;
            }
            for (int j = 0; j < GetNtaxa(); j++) {
                int fromgroup = 0;
                for (int k = 0; k < n; k++) {
                    if (GetTaxonSet()->GetTaxon(j) == group[k]) {
                        fromgroup = 1;
                    }
                }
                if (fromgroup == 0) {
                    if (Data[j][i] != unknown) {
                        present[Data[j][i]] = 1;
                    }
                }
            }

            int pos = 0;
            for (int j = 0; j < n; j++) {
                if (Data[GetTaxonSet()->GetTaxonIndex(group[j])][i] != unknown) {
                    if (present[Data[GetTaxonSet()->GetTaxonIndex(group[j])][i]] == 0) {
                        Data[GetTaxonSet()->GetTaxonIndex(group[j])][i] = unknown;
                        nrem++;
                        pos = 1;
                    }
                }
            }
            if (pos != 0) {
                npos++;
            }
        }
        std::cerr << "number of positions masked   : " << npos << '\n';
        std::cerr << "total number of cells masked : " << nrem << '\n';
    }

    SequenceAlignment(SequenceAlignment *from, double missingfrac) {
        int n = 0;
        for (int i = 0; i < from->GetNsite(); i++) {
            if (from->MissingFrac(i) < missingfrac) {
                n++;
            }
        }
        statespace = from->statespace;
        taxset = from->taxset;
        Ntaxa = from->GetNtaxa();
        Nsite = n;
        Data = new int *[Ntaxa];
        for (int k = 0; k < Ntaxa; k++) {
            Data[k] = new int[Nsite];
        }
        int j = 0;
        for (int i = 0; i < from->GetNsite(); i++) {
            if (from->MissingFrac(i) < missingfrac) {
                for (int k = 0; k < Ntaxa; k++) {
                    Data[k][j] = from->Data[k][i];
                }
                j++;
            }
        }
    }

    void Mask(SequenceAlignment *from) {
        for (int i = 0; i < from->GetNsite(); i++) {
            for (int j = 0; j < Ntaxa; j++) {
                if (from->Data[j][i] == unknown) {
                    Data[j][i] = unknown;
                }
            }
        }
    }

    bool GroupCompatible(int i, std::map<std::string, int> &group, int ngroup) {
        std::vector<int> presence(ngroup);
        for (int n = 0; n < ngroup; n++) {
            presence[n] = 0;
        }
        for (int k = 0; k < GetNtaxa(); k++) {
            if (Data[k][i] != unknown) {
                //             std::cerr << GetTaxonSet()->GetTaxon(k) << '\t' <<
                //             group[GetTaxonSet()->GetTaxon(k)]
                // << '\n';
                presence[group[GetTaxonSet()->GetTaxon(k)]] = 1;
            }
        }
        int ok = 1;
        for (int n = 1; n < ngroup; n++) {
            ok *= presence[n];
        }
        return ok != 0;
    }

    double MissingFrac(int i) {
        double n = 0;
        for (int k = 0; k < GetNtaxa(); k++) {
            if (Data[k][i] != unknown) {
                n++;
            }
        }
        return n / GetNtaxa();
    }

    SequenceAlignment(int **inData, std::string *names, int inNsite, StateSpace *instatespace,
                      TaxonSet *intaxset) {
        Nsite = inNsite;
        taxset = intaxset;
        Ntaxa = taxset->GetNtaxa();
        statespace = instatespace;

        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }
        for (int i = 0; i < Ntaxa; i++) {
            int mapi = taxset->GetTaxonIndex(names[i]);
            for (int j = 0; j < Nsite; j++) {
                Data[mapi][j] = inData[i][j];
            }
        }
    }

    virtual ~SequenceAlignment() = default;

    // the set of characters (A,C,G,T for nucleotides, etc..)
    StateSpace *GetStateSpace() { return statespace; }

    void RegisterWith(TaxonSet *intaxset) {
        std::cerr << "register\n";
        if (taxset->GetNtaxa() != intaxset->GetNtaxa()) {
            std::cerr << "error in register seq\n";
            exit(1);
        }
        auto tmp = new int *[GetNtaxa()];
        for (int i = 0; i < GetNtaxa(); i++) {
            int k = 0;
            while ((k < GetNtaxa()) && (taxset->GetTaxon(k) != intaxset->GetTaxon(i))) {
                k++;
            }
            if (k == GetNtaxa()) {
                std::cerr << "error in register seq : overflow\n";
                exit(1);
            }
            tmp[i] = Data[k];
        }
        for (int i = 0; i < GetNtaxa(); i++) {
            Data[i] = tmp[i];
        }
        delete taxset;
        taxset = intaxset;
        delete[] tmp;
        std::cerr << "register ok\n";
    }

    void Unclamp() {
        for (int i = 0; i < Ntaxa; i++) {
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = unknown;
            }
        }
    }

    int GetNstate() { return statespace->GetNstate(); }

    // the list of taxa
    TaxonSet *GetTaxonSet() { return taxset; }

    int GetNsite() { return Nsite; }

    int GetNtaxa() { return taxset->GetNtaxa(); }

    bool isMissing(int taxon, int site) { return Data[taxon][site] == -1; }

    bool AllMissingTaxon(int tax) {
        bool ret = true;
        int site = 0;
        while ((site < GetNsite()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            site++;
        }
        return ret;
    }

    bool AllMissingTaxon(std::string taxname) {
        int index = taxset->GetTaxonIndex(taxname);
        if (index == -1) {
            std::cerr << "error in all missing taxon: did not recognize " << taxname << '\n';
            exit(1);
        }
        return AllMissingTaxon(index);
    }

    bool AllMissingColumn(int site) {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            tax++;
        }
        return ret;
    }

    bool NoMissingColumn(int site) {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] != unknown);
            tax++;
        }
        return ret;
    }

    bool ConstantColumn(int site) {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && (Data[tax][site] == unknown)) {
            tax++;
        }

        if (tax < GetNtaxa()) {
            int refstate = Data[tax][site];

            while ((tax < GetNtaxa()) && ret) {
                if (Data[tax][site] != -1) {
                    ret &= static_cast<int>(Data[tax][site] == refstate);
                }
                tax++;
            }
        }
        return ret;
    }

    void SetState(int taxon, int site, int state) { Data[taxon][site] = state; }

    int GetState(int taxon, int site) { return Data[taxon][site]; }

    void GetEmpiricalFreq(double *in);

    void GetSiteEmpiricalFreq(double **in, double pseudocount = 0);

    void ToStream(std::ostream &os);
    void ToStreamTriplet(std::ostream &os);
    int GetNonMissingTriplet();
    void ToFasta(std::ostream &os);

    void DeleteConstantSites() {
        int i = 0;
        int j = 0;
        int Eliminated = 0;
        while (i < Nsite) {
            int k = 0;
            while ((k < Ntaxa) && (Data[k][i] == unknown)) {
                k++;
            }
            if (k < Ntaxa) {
                int a = Data[k][i];
                k++;
                while ((k < Ntaxa) && ((Data[k][i] == unknown) || (Data[k][i] == a))) {
                    k++;
                }
                if (k == Ntaxa) {
                    Eliminated++;
                } else {
                    for (int k = 0; k < Ntaxa; k++) {
                        Data[k][j] = Data[k][i];
                    }
                    j++;
                }
            }
            i++;
        }
        Nsite -= Eliminated;
        std::cout << "number of positions eliminated : " << Eliminated << '\n';
    }

    void GetMissingCellsFromTemplate(SequenceAlignment *from) {
        int fromnsite = from->GetNsite();

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int ii = i % fromnsite;
                int state = from->GetState(j, ii);
                if (state == unknown) {
                    Data[j][i] = unknown;
                }
            }
        }
    }

    double GetMeanDiversity() {
        int Nstate = GetNstate();
        std::vector<int> found(Nstate);

        double mean = 0;
        for (int i = 0; i < Nsite; i++) {
            for (int k = 0; k < Nstate; k++) {
                found[k] = 0;
            }
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    found[state] = 1;
                }
            }
            double div = 0;
            for (int k = 0; k < Nstate; k++) {
                div += found[k];
            }
            mean += div;
        }
        mean /= Nsite;
        return mean;
    }

    void GetTaxEmpiricalFreq(double **taxfreq) {
        int Nstate = GetNstate();
        for (int j = 0; j < Ntaxa; j++) {
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] = 0;
            }
        }

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    taxfreq[j][state]++;
                }
            }
        }

        for (int j = 0; j < Ntaxa; j++) {
            double total = 0;
            for (int k = 0; k < Nstate; k++) {
                total += taxfreq[j][k];
            }
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] /= total;
            }
        }
    }

    double CompositionalHeterogeneity(std::ostream *os) {
        int Nstate = GetNstate();
        auto taxfreq = new double *[Ntaxa];
        for (int j = 0; j < Ntaxa; j++) {
            taxfreq[j] = new double[Nstate];
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] = 0;
            }
        }

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    taxfreq[j][state]++;
                }
            }
        }

        // make global freqs out of tax-specific freqs
        auto globalfreq = new double[Nstate];
        for (int k = 0; k < Nstate; k++) {
            globalfreq[k] = 0;
            for (int j = 0; j < Ntaxa; j++) {
                globalfreq[k] += taxfreq[j][k];
            }
        }

        // normalise
        double total = 0;
        for (int k = 0; k < Nstate; k++) {
            total += globalfreq[k];
        }
        for (int k = 0; k < Nstate; k++) {
            globalfreq[k] /= total;
        }
        for (int j = 0; j < Ntaxa; j++) {
            double total = 0;
            for (int k = 0; k < Nstate; k++) {
                total += taxfreq[j][k];
            }
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] /= total;
                if (os != nullptr) {
                    (*os) << taxfreq[j][k] << '\t';
                }
            }
            if (os != nullptr) {
                (*os) << '\n';
            }
        }
        if (os != nullptr) {
            (*os) << '\n';
        }

        // compute max distance
        double maxdist = 0;
        for (int j = 0; j < Ntaxa; j++) {
            double dist = 0;
            for (int k = 0; k < Nstate; k++) {
                double tmp = (taxfreq[j][k] - globalfreq[k]);
                dist += tmp * tmp;
            }
            if (maxdist < dist) {
                maxdist = dist;
            }
        }

        delete[] globalfreq;
        for (int j = 0; j < Ntaxa; j++) {
            delete[] taxfreq[j];
        }
        delete[] taxfreq;

        return maxdist;
    }

    // data fields

    int Ntaxa;
    int Nsite;
    TaxonSet *taxset;
    StateSpace *statespace;
    int **Data;
};

class FileSequenceAlignment : public SequenceAlignment {
  public:
    FileSequenceAlignment(std::istream &is);
    FileSequenceAlignment(std::string filename, int fullline = 0);

  private:
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int ReadNexus(std::string filespec);
    int ReadSpecial(std::string filename);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);

    std::string *SpeciesNames;
};

#endif  // SEQUENCEALIGNMENT_H
