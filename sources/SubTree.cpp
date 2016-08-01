

#include "SequenceAlignment.h"
#include "ContinuousData.h"
#include "Tree.h"
#include <fstream>
#include <cstdlib>
using namespace std;


int main(int argc, char* argv[])	{

	/*
	string datafile1 = argv[1];
	string datafile2 = argv[2];

	FileSequenceAlignment* data1 = new FileSequenceAlignment(datafile1);
	FileSequenceAlignment* data2 = new FileSequenceAlignment(datafile2);

	TaxonSet* tax1 = data1->GetTaxonSet();
	TaxonSet* tax2 = data2->GetTaxonSet();

	int tot = 0;
	for (int i=0; i<tax1->GetNtaxa(); i++)	{
		if (tax2->GetTaxonIndex(tax1->GetTaxon(i)) != -1)	{
			cerr << tax1->GetTaxon(i) << '\n';
			tot++;
		}
	}

	cerr << '\n';
	cerr << "total : " << tot << '\n';
	cerr << '\n';

	*/



	string taxfile = argv[1];
	string datafile = argv[2];

	int Ntaxa = 0;
	ifstream is(taxfile.c_str());
	is >> Ntaxa;

	string* names = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		is >> names[i];
	}

	FileSequenceAlignment* ali = new FileSequenceAlignment(datafile);
	for (int j=0; j<ali->GetNtaxa(); j++)	{
		int k = 0;
		while ((k < Ntaxa) && (names[k] != ali->GetTaxonSet()->GetTaxon(j))) k++;
		if (k != Ntaxa)	{
			cout << ali->GetTaxonSet()->GetTaxon(j) << '\t';
			for (int i=0; i<ali->GetNsite(); i++)	{
				if (ali->GetState(j,i) != unknown)	{
					cout << ali->GetStateSpace()->GetState(ali->GetState(j,i));
				}
				else	{
					cout << '?';
				}
			}
			cout << '\n';
		}
	}


	/*
	string datafile = argv[1];
	string treefile = argv[2];
	string tax1 = argv[3];
	string tax2 = argv[4];

	FileSequenceAlignment* data = new FileSequenceAlignment(datafile);
	Tree * tree = new Tree(treefile);

	cerr << "taxonset\n";
	TaxonSet* taxset = new TaxonSet(tree, tree->GetLCA(tax1,tax2));

	taxset->ToStream(cout);

	cerr << "ok\n";

	SequenceAlignment* newdata = new SequenceAlignment(data,taxset);
	string name = argv[5];
	ofstream os((name + ".ali").c_str());
	newdata->ToStream(os);

	cerr << "tree\n";

	ofstream tos((name + ".tree").c_str());
	tree->ToStream(tos,tree->GetLCA(tax1,tax2));
	tos.close();
	*/

	/*
	// get data from file
	string datafile = argv[1];
	string contdatafile = argv[2];
	ofstream os(argv[3]);

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	TaxonSet* taxonset = data->GetTaxonSet();

	ContinuousData* contdata = new FileContinuousData(contdatafile);
	contdata->ToStream(os,taxonset);
	*/

	/*
	string contdatafile = argv[1];
	ofstream os(argv[2]);
	ContinuousData* contdata = new FileContinuousData(contdatafile);
	contdata->ToStreamLog(os);
	*/
}

