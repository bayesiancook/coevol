#include "CodonSequenceAlignment.h"


int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N;
	is >> N;
	string* gene = new string[N];
	CodonSequenceAlignment** ali = new CodonSequenceAlignment*[N];
	FileSequenceAlignment** nucali = new FileSequenceAlignment*[N];
	double* freq = new double[Nnuc];
	double* comp = new double[N];
	double* gc = new double[N];
	for (int i=0; i<N; i++)	{
		is >> gene[i];
		cerr << gene[i] << '\n';
		nucali[i] = new FileSequenceAlignment(gene[i]);
		ali[i] = new CodonSequenceAlignment(nucali[i],true);
		double tmp = ali[i]->NucleotideCompositionalHeterogeneity(0);
		nucali[i]->GetEmpiricalFreq(freq);
		gc[i] = freq[1] + freq[2];
		comp[i] = tmp;
	}

	for (int i=0; i<N; i++)	{
		cerr << gene[i] << '\t' << gc[i] << '\t' << comp[i] << '\n';
	}

	int* permut = new int[N];
	for (int i=0; i<N; i++)	{
		permut[i] = i;
	}
	for (int i=0; i<N; i++)	{
		for (int j=N-1; j>i; j--)	{
			if (gc[permut[i]] > gc[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}
	ofstream os("sort.gc");
	for (int i=0; i<N; i++)	{
		os << gene[permut[i]] << '\t' << gc[permut[i]] << '\t' << comp[permut[i]] << '\n';
	}
	cerr << '\n';

	for (int i=0; i<N; i++)	{
		permut[i] = i;
	}
	for (int i=0; i<N; i++)	{
		for (int j=N-1; j>i; j--)	{
			if (comp[permut[i]] > comp[permut[j]])	{
				int tmp = permut[i];
				permut[i] = permut[j];
				permut[j] = tmp;
			}
		}
	}
	ofstream pos("sort.comp");
	for (int i=0; i<N; i++)	{
		pos << gene[permut[i]] << '\t' << comp[permut[i]] << '\t' << gc[permut[i]] << '\n';
	}
}

