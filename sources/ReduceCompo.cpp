
#include "SumConstrained.h"

int main(int argc, char* argv[])	{

	string infile = argv[1];
	string outfile = argv[2];

	int anc = 1;
	// int anc = atoi(argv[3]);

	ifstream is(infile.c_str());
	int Ntaxa, Nsite;
	is >> Ntaxa >> Nsite;

	double* x = new double[Nsite];
	double* y = new double[Nsite];

	SumConstrainedMapping* basis = new SumConstrainedMapping(Nsite);

	ofstream bos((outfile + ".basis").c_str());
	basis->ToStream(bos);
	bos.close();

	ofstream os(outfile.c_str());
	os << Ntaxa << '\t' << Nsite-1 << '\n';

	double** b = basis->base;
	for (int k=0; k<Ntaxa; k++)	{
		string name1, name2;
		is >> name1;
		if (anc)	{
			is >> name2;
		}
		double total = 0;
		for (int i=0; i<Nsite; i++)	{
			double tmp;
			is >> tmp;
			x[i] = log(tmp);
			total += x[i];
		}
		total /= Nsite;
		for (int i=0; i<Nsite; i++)	{
			x[i] -= total;
		}
		for (int i=0; i<Nsite; i++)	{
			double tmp = 0;
			for (int j=0; j<Nsite; j++)	{
				tmp += b[i][j] * x[j];
			}
			y[i] = tmp;
			if (!i)	{
				if (fabs(y[i]) > 1e-6)	{
					cerr << "error : " << y[i] << '\n';
				}
			}
		}
		os << name1;
		if (anc)	{
			os << '\t' << name2;
		}
		for (int i=1; i<Nsite; i++)	{
			os << '\t' << exp(y[i]);
		}
		os << '\n';
	}
	os.close();
}





