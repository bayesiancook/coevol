
#include "MG3OmegaCodonSubMatrix.h"


void MG3CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				Q[i][j] = (*GetNucMatrix(pos))(a,b);
				total += Q[i][j];
			}
			else	{
				Q[i][j] = 0;
			}
		}
	}
	Q[i][i] = -total;
}

void MG3CodonSubMatrix::ComputeStationary()	{

	// compute stationary probabilities
	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] = GetNucMatrix(0)->Stationary(GetCodonPosition(0,i)) * GetNucMatrix(1)->Stationary(GetCodonPosition(1,i)) * GetNucMatrix(2)->Stationary(GetCodonPosition(2,i));
		total += mStationary[i];
	}

	// renormalize stationary probabilities
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i]  /= total;
	}
}

void MG3OmegaCodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = (*GetNucMatrix(pos))(a,b);
				if (! Synonymous(i,j))	{
					Q[i][j] *= GetOmega();
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

