#ifndef CONTINUOUSDATA_H
#define CONTINUOUSDATA_H

#include <cmath>
#include <fstream>
#include <sstream>
#include "TaxonSet.hpp"

// this class works like an interface
// it does not do any job
class ContinuousData {
  public:
    ContinuousData() = default;

    ContinuousData(const TaxonSet* intaxset, int inNsite) {
        taxset = intaxset;
        Nsite = inNsite;
        Data = new double*[GetNtaxa()];
        for (int i = 0; i < GetNtaxa(); i++) {
            Data[i] = new double[Nsite];
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = 0;
            }
        }
        charname = new std::string[GetNsite()];
        for (int j = 0; j < Nsite; j++) {
            charname[j] = "none";
        }
    }

    ContinuousData(ContinuousData* from) {
        taxset = from->GetTaxonSet();
        Nsite = from->GetNsite();
        Data = new double*[GetNtaxa()];
        for (int i = 0; i < GetNtaxa(); i++) {
            Data[i] = new double[Nsite];
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = from->Data[i][j];
            }
        }
        charname = new std::string[GetNsite()];
        for (int j = 0; j < Nsite; j++) {
            charname[j] = from->charname[j];
        }
    }

    void ToStream(std::ostream& os, TaxonSet* taxset = nullptr) {
        if (taxset == nullptr) {
            //             std::cerr << "??? in taxon set\n";
            os << GetNtaxa() << '\t' << GetNsite() << '\n';
            for (int i = 0; i < GetNtaxa(); i++) {
                os << GetTaxonSet()->GetTaxon(i);
                /*
                  std::string s = taxset->GetTaxon(i);
                  unsigned int l = s.length();
                  unsigned int k = 0;
                  while ((k < l) && (s[k] != '_')) k++;
                  if (k == l)	{
                              std::cerr << "error in get name\n";
                  exit(1);
                  }
                  k++;
                  os << s.substr(k,l-k);
                */
                for (int j = 0; j < GetNsite(); j++) {
                    os << '\t' << Data[i][j];
                    /*
                      if (Data[i][j] == -1)	{
                      os << '\t' << -1;
                      }
                      else	{
                      os << '\t' << log(Data[i][j]);
                      // os << '\t' << exp(Data[i][j]);
                      }
                    */
                }
                os << '\n';
            }

        } else {
            int ntaxa = 0;
            for (int i = 0; i < GetNtaxa(); i++) {
                if (taxset->GetTaxonIndex(GetTaxonSet()->GetTaxon(i)) != -1) {
                    ntaxa++;
                }
            }

            os << ntaxa << '\t' << GetNsite() << '\n';
            for (int i = 0; i < GetNtaxa(); i++) {
                if (taxset->GetTaxonIndex(GetTaxonSet()->GetTaxon(i)) != -1) {
                    os << GetTaxonSet()->GetTaxon(i);
                    for (int j = 0; j < GetNsite(); j++) {
                        os << '\t' << Data[i][j];
                    }
                    os << '\n';
                }
            }
        }
    }

    void ToStreamLog(std::ostream& os) {
        os << GetNtaxa() << '\t' << GetNsite() << '\n';
        for (int i = 0; i < GetNtaxa(); i++) {
            os << GetTaxonSet()->GetTaxon(i);
            /*
              std::string s = taxset->GetTaxon(i);
              unsigned int l = s.length();
              unsigned int k = 0;
              while ((k < l) && (s[k] != '_')) k++;
              if (k == l)	{
                          std::cerr << "error in get name\n";
              exit(1);
              }
              k++;
              os << s.substr(k,l-k);
            */
            for (int j = 0; j < GetNsite(); j++) {
                if (Data[i][j] == -1) {
                    os << '\t' << -1;
                } else {
                    os << '\t' << log(Data[i][j]);
                }
            }
            os << '\n';
        }
    }

    // the list of taxa
    const TaxonSet* GetTaxonSet() { return taxset; }

    int GetNtaxa() { return GetTaxonSet()->GetNtaxa(); }

    int GetNsite() { return Nsite; }

    bool isMissing(int taxon, int site) { return Data[taxon][site] == -1; }

    bool isMissing(int taxon) {
        bool mis = true;
        for (int i = 0; i < Nsite; i++) {
            mis &= static_cast<int>(Data[taxon][i] == -1);
        }
        return mis;
    }

    bool isMissing(std::string taxname) {
        int index = GetTaxonSet()->GetTaxonIndex(taxname);
        return (index == -1);
    }
    /*
      if (index == -1)	{
                  std::cerr << "error : taxon not found : " << taxname << '\n';
      GetTaxonSet()->ToStream(            std::cerr);
      exit(1);
      }
      return isMissing(index);
      }
    */

    double GetMeanLog(int site) {
        double total = 0;
        int n = 0;
        for (int j = 0; j < GetNtaxa(); j++) {
            if (Data[j][site] != -1) {
                total += log(Data[j][site]);
                n++;
            }
        }
        return total / n;
    }

    std::string GetCharacterName(int site) { return charname[site]; }

    double GetState(std::string taxon, int site) {
        return Data[taxset->GetTaxonIndex(taxon)][site];
    }

    // char       GetCharState(int taxon, int state);
    double GetState(int taxon, int site) { return Data[taxon][site]; }
    // char       GetCharState(int taxon, int state);

    int Nsite;
    const TaxonSet* taxset;
    double** Data;
    std::string* charname;
};

class FileContinuousData : public ContinuousData {
  public:
    FileContinuousData(std::istream& is) { ReadDataFromFile(is); }

    FileContinuousData(std::string filename) {
        std::ifstream is(filename.c_str());
        if (!is) {
            std::cerr << "error when opening file : " << filename << '\n';
            exit(1);
        }
        ReadDataFromFile(is);
    }

  private:
    int ReadDataFromFile(std::istream& is) {
        std::string temp;
        is >> temp;
        int Ntaxa;
        if (temp == "#TRAITS") {
            is >> Ntaxa;
            is >> Nsite;
            charname = new std::string[Nsite];
            for (int j = 0; j < Nsite; j++) {
                is >> charname[j];
            }
        } else {
            Ntaxa = atoi(temp.c_str());
            is >> Nsite;
            charname = new std::string[Nsite];
            for (int j = 0; j < Nsite; j++) {
                std::ostringstream s;
                s << "character" << j + 1;
                charname[j] = s.str();
            }
        }


        auto name = new std::string[Ntaxa];
        Data = new double*[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            is >> name[i];
            Data[i] = new double[Nsite];
            for (int j = 0; j < Nsite; j++) {
                is >> Data[i][j];
            }
        }
        taxset = new TaxonSet(name, Ntaxa);
        delete[] name;
        return 1;
    }
};

#endif