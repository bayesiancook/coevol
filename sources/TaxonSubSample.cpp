

#include "Random.h"
#include "SequenceAlignment.h"
#include "ContinuousData.h"
#include "Tree.h"
#include <fstream>
#include <cstdlib>
using namespace std;


int main(int argc, char* argv[])	{


	string datafile = argv[1];
	double p = atof(argv[2]);
	string name = argv[3];

	FileSequenceAlignment* data = new FileSequenceAlignment(datafile);
	int member[data->GetNtaxa()];
	int ntaxa = 0;
	for (int i=0; i<data->GetNtaxa(); i++)	{
		if (Random::Uniform() < p)	{
			member[i] = 1;
			ntaxa++;
		}
		else	{
			member[i] = 0;
		}
	}
	string* names = new string[ntaxa];
	int k= 0;
	for (int i=0; i<data->GetNtaxa(); i++)	{
		if (member[i])	{
			names[k] = data->GetTaxonSet()->GetTaxon(i);
			k++;
		}
	}

	TaxonSet* taxset = new TaxonSet(names,ntaxa);
	SequenceAlignment* newdata = new SequenceAlignment(data,taxset);
	ofstream os((name + ".ali").c_str());
	newdata->ToStream(os);

}

