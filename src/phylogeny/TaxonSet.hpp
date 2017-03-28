#ifndef TAXONSET_H
#define TAXONSET_H

#include <iostream>
#include <map>
#include <string>

class Tree;
class Link;

class TaxonSet {
  public:
    TaxonSet(const std::string *names, int ntaxa);
    TaxonSet(const Tree *tree, const Link *subgroup = nullptr);
    ~TaxonSet();

    int GetNtaxa() const;
    std::string GetTaxon(int index) const;
    int GetTaxonIndex(std::string intaxon) const;
    int GetTaxonIndexWithIncompleteName(std::string taxname) const;

    void ToStream(std::ostream &os) const;

  private:
    void RecursiveGetSubSet(const Link *from, int &i);

    int Ntaxa;
    mutable std::map<std::string, int> taxmap;
    std::string *taxlist;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline int TaxonSet::GetNtaxa() const { return Ntaxa; }
inline std::string TaxonSet::GetTaxon(int index) const { return taxlist[index]; }
inline int TaxonSet::GetTaxonIndex(std::string intaxon) const { return taxmap[intaxon] - 1; }

#endif  // TAXONSET_H
