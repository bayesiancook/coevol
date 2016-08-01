#include "Random.h"
#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	SequenceAlignment* data = new FileSequenceAlignment(datafile);

	int dim = data->GetNstate();
	int N = data->GetNsite();

	double** freq = new double*[N];
	for (int i=0; i<N; i++)	{
		freq[i] = new double[dim];
	}
	data->GetSiteEmpiricalFreq(freq,-1.0/dim);

	ofstream os(argv[2]);
	os << N << '\t' << dim << '\n';
	for (int i=0; i<N; i++)	{
		for (int j=0; j<dim; j++)	{
			os << freq[i][j] << '\t';
		}
		os << '\n';
	}
}

