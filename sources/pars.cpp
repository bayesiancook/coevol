
#include "Random.h"
#include "Tree.h"
#include "SequenceAlignment.h"

SequenceAlignment* data;
Tree* tree;

int pars(const Link* from, int* obs, int site)	{

	int ret = 0;
	if (from->isLeaf())	{
		int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
		int state = data->GetState(tax,site);
		if (state == -1)	{
			for (int k=0; k<data->GetNstate(); k++)	{
				obs[k] = 1;
			}
		}
		else	{
			for (int k=0; k<data->GetNstate(); k++)	{
				obs[k] = 0;
			}
			obs[state] = 1;
		}
	}
	else	{
		for (int k=0; k<data->GetNstate(); k++)	{
			obs[k] = 1;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		}
	}

}

int main(int argc, char* argv[])	{

	SequenceAlignment* data = new FileSequenceAlignment(argv[1]);
	const TaxonSet* taxonset = data->GetTaxonSet();
	Tree* tree = new Tree(argv[2]);
	tree->RegisterWith(taxonset);

	cerr << data->GetNsite() << '\t' << data->GetNtaxa() << '\n';
	for (int i=0; i<data->GetNsite(); i++)	{
	}
}




