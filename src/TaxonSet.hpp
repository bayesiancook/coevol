#ifndef TAXONSET_H
#define TAXONSET_H

#include <map>
#include <string>
#include <iostream>

using namespace std;


class Tree;
class Link;

class TaxonSet	{

	public:
				TaxonSet(const string* names, int ntaxa);
				TaxonSet(const Tree* tree, const Link* subgroup = 0);
				~TaxonSet();

	int 			GetNtaxa() const;
	string 			GetTaxon(int index) const;
	int 			GetTaxonIndex(string intaxon) const;
	int 			GetTaxonIndexWithIncompleteName(string intaxon) const;

	void			ToStream(ostream& os) const;
	private:

	void			RecursiveGetSubSet(const Link* from, int& index);

	int Ntaxa;
	mutable map<string,int> taxmap;
	string*	taxlist;

};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


inline int TaxonSet::GetNtaxa() const {return Ntaxa;}
inline string TaxonSet::GetTaxon(int index) const {return taxlist[index];}
inline int TaxonSet::GetTaxonIndex(string intaxon) const {return taxmap[intaxon]-1;}

#endif // TAXONSET_H
