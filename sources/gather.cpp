
#include "Random.h"
#include "BiologicalSequences.h"

int main(int argc, char* argv[])	{

	ifstream is1(argv[1]);
	ifstream is2(argv[2]);
	ifstream is3(argv[3]);
	int Nsite = atoi(argv[4]);
	double cutoff1 = atof(argv[5]);
	double cutoff2 = atof(argv[6]);
	double cutoff3 = atof(argv[7]);
	double cutoff4 = atof(argv[8]);
	
	int offset = 22;

	double om1[Nsite];
	double ppom1[Nsite];
	double om3[Nsite][3];
	double ppom3[Nsite][3];
	double sel[Nsite][20];
	double ppsel[Nsite][20];
	
	for (int i=0; i<Nsite; i++)	{
		int tmp;
		is1 >> tmp;
		if (tmp != i)	{
			cerr << "error when reading om1\n";
			exit(1);
		}
		is1 >> om1[i] >> ppom1[i];
	}

	for (int i=0; i<Nsite; i++)	{
		int tmp;
		is2 >> tmp;
		if (tmp != i)	{
			cerr << "error when reading om1\n";
			exit(1);
		}
		for (int k=0; k<3; k++)	{
			is2 >> om3[i][k] >> ppom3[i][k];
		}
	}

	for (int i=0; i<Nsite; i++)	{
		int tmp;
		is3 >> tmp;
		if (tmp != i)	{
			cerr << "error when reading om1\n";
			exit(1);
		}
		for (int a=0; a<20; a++)	{
			is3 >> sel[i][a] >> ppsel[i][a];
		}
	}

	for (int i=0; i<Nsite; i++)	{
		int included = 0;
		if (ppom1[i] > cutoff1)	{
			included = 1;
		}
		for (int k=0; k<3; k++)	{
			if (ppom3[i][k] > cutoff2)	{
				included = 1;
			}
		}
		for (int a=0; a<20; a++)	{
			if ((ppsel[i][a] > cutoff3) || (ppsel[i][a] < (1-cutoff3)))	{
				included = 1;
			}
		}

		if (included)	{
			cout << i + offset << '\t';
			cout << om1[i] << '\t' << "(" << ppom1[i] << ")" << '\t';
			for (int k=2; k>=0; k--)	{
				cout << om3[i][k] << '\t' << "(" << ppom3[i][k] << ")" << '\t';
			}
			for (int a=0; a<20; a++)	{
				if (ppsel[i][a] > cutoff4)	{
					cout << "+" << AminoAcids[a] << " (" << ppsel[i][a] << ") " << '\t';
				}
				if (ppsel[i][a] < (1-cutoff4))	{
					cout << "-" << AminoAcids[a] << " (" << 1-ppsel[i][a] << ") " << '\t';
				}
			}
			cout << '\n';
		}
	}
}

