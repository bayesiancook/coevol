#include "Random.h"
#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	ifstream is(argv[2]);
	ofstream os(argv[3]);

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	int N;
	is >> N;
	string* taxa = new string[N];
	
	for (int i=0; i<N; i++)	{
		is >> taxa[i];
	}
	TaxonSet* taxset = new TaxonSet(taxa,N);

	SequenceAlignment* subdata = new SequenceAlignment(data,taxset);
	cerr << data->GetNtaxa() << '\t' << subdata->GetNtaxa() << '\n';
	subdata->ToStream(os);
}

