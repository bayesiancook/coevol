#include "Random.h"
#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string name = argv[1];
	SequenceAlignment* data = new FileSequenceAlignment(name);
	int Ntaxa = data->GetNtaxa();
	double** freq = new double*[Ntaxa];
	if (data->GetNstate() != Nnuc)	{
		cerr << "error: only nuc ali\n";
		exit(1);
	}
	for (int i=0; i<Ntaxa; i++)	{
		freq[i] = new double[Nnuc];
	}
	data->GetTaxEmpiricalFreq(freq);
	
	ofstream os(argv[2]);
	os << Ntaxa << '\t' << 1 << '\n';
	for (int i=0; i<Ntaxa; i++)	{
		double gc = freq[i][1] + freq[i][2];
		os << data->GetTaxonSet()->GetTaxon(i) << '\t' << gc / (1-gc) << '\n';
	}
}
