
#include <iostream>
#include <cstdlib>
using namespace std;

#include "BiologicalSequences.h"
#include "TaxonSet.h"
#include "Tree.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 TaxonSet
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TaxonSet::TaxonSet(const string* names, int ntaxa)	{

	Ntaxa = ntaxa;
	taxlist = new string[ntaxa];
	for (int i=0; i<ntaxa; i++)	{
		if (taxmap[names[i]])	{
			cerr << "found several taxa with same name : " << names[i] << '\n';
			exit(1);
		}
		taxlist[i] = names[i];
		taxmap[names[i]] = i+1;
	}
}

TaxonSet::~TaxonSet()	{
	delete[] taxlist;
}

TaxonSet::TaxonSet(const Tree* tree, const Link* subgroup)	{

	Ntaxa = tree->GetSize(subgroup);
	taxlist = new string[Ntaxa];
	if (!subgroup)	{
		subgroup = tree->GetRoot();
	}
	int i = 0;
	RecursiveGetSubSet(subgroup,i);
}

void TaxonSet::RecursiveGetSubSet(const Link* from, int& i)	{

	if (from->isLeaf())	{
		taxlist[i] = from->GetNode()->GetName();
		taxmap[from->GetNode()->GetName()] = i+1;
		i++;
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetSubSet(link->Out(),i);
		}
	}
}

void TaxonSet::ToStream(ostream& os)	const {
	os << Ntaxa << '\n';
	for (int i=0; i<Ntaxa; i++)	{
		os << taxlist[i] << '\n';
	}
}

int TaxonSet::GetTaxonIndexWithIncompleteName(string taxname) const {

	int found = -1;
	for (int i=0; i<Ntaxa; i++)	{
		if (taxlist[i].substr(0,taxname.length()) == taxname)	{
			if (found != -1)	{
				cerr << "error : taxon found twice : " << taxname << '\n';
				exit(1);
			}
			found = i;
		}
	}
	if (found == -1)	{
		for (int i=0; i<Ntaxa; i++)	{
			if (taxname.substr(0,taxlist[i].length()) == taxlist[i])	{
				if (found != -1)	{
					cerr << "error : taxon found twice : " << taxname << '\n';
					exit(1);
				}
				found = i;
			}
		}
	}
	return found;
}

