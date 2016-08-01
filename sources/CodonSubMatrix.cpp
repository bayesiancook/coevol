#include "CodonSubMatrix.h"

void MGCodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				Q[i][j] = (*NucMatrix)(a,b);
				total += Q[i][j];
			}
			else	{
				Q[i][j] = 0;
			}
		}
	}
	Q[i][i] = -total;
}

void MGCodonSubMatrix::ComputeStationary()	{

	// compute stationary probabilities
	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
		total += mStationary[i];
	}

	// renormalize stationary probabilities
	for (int i=0; i<GetNstate(); i++)	{
		mStationary[i]  /= total;
	}
}

void MGCodonSubMatrix::ComputeNucArrays()	{

	cerr << "error: in mg codon sub matrix compute nuc array\n";
	exit(1);
	for (int i=0; i<Nnuc; i++)	{
		for (int j=0; j<Nnuc; j++)	{
			if (i != j)	{
				synnucarray[i][j] = (*NucMatrix)(i,j);
				nonsynnucarray[i][j] = (*NucMatrix)(i,j);
			}
		}
	}
	const double* stat = NucMatrix->GetStationary();
	const int* stoppos1 = GetCodonStateSpace()->GetStopPos1();
	const int* stoppos2 = GetCodonStateSpace()->GetStopPos2();
	const int* stoppos3 = GetCodonStateSpace()->GetStopPos3();
	stopstat = 0;
	for (int i=0; i<GetCodonStateSpace()->GetNstop(); i++)	{
		stopstat += stat[stoppos1[i]] * stat[stoppos2[i]] * stat[stoppos3[i]];
	}
}

void MGOmegaCodonSubMatrix::ComputeArray(int i)	{

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
				Q[i][j] = (*NucMatrix)(a,b);
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


void MGOmegaCodonSubMatrix::ComputeNucArrays()	{

	for (int i=0; i<Nnuc; i++)	{
		for (int j=0; j<Nnuc; j++)	{
			if (i != j)	{
				synnucarray[i][j] = (*NucMatrix)(i,j);
				nonsynnucarray[i][j] = omega * (*NucMatrix)(i,j);
			}
		}
	}
	const double* stat = NucMatrix->GetStationary();
	const int* stoppos1 = GetCodonStateSpace()->GetStopPos1();
	const int* stoppos2 = GetCodonStateSpace()->GetStopPos2();
	const int* stoppos3 = GetCodonStateSpace()->GetStopPos3();
	stopstat = 0;
	for (int i=0; i<GetCodonStateSpace()->GetNstop(); i++)	{
		stopstat += stat[stoppos1[i]] * stat[stoppos2[i]] * stat[stoppos3[i]];
	}
}


void AminoAcidReducedCodonSubMatrix::ComputeStationary()	{

	// Stat[a] = \sum_{i|a} CodonStat[i]

	for (int a=0; a<GetNstate(); a++)	{
		mStationary[a] = 0;
	}

	for (int i=0; i<GetCodonStateSpace()->GetNstate(); i++)	{
		mStationary[GetCodonStateSpace()->Translation(i)] += GetCodonSubMatrix()->Stationary(i);
	}

}

void AminoAcidReducedCodonSubMatrix::ComputeArray(int a)	{

	// Q[a][[b] = [ \sum _{i|a, j|b} CodonQ[i][j] ] / [ \sum_{i|a} CodonStat[i] ]

	for (int b=0; b<GetNstate(); b++)	{
		Q[a][b] = 0;
	}

	for (int i=0; i<GetCodonStateSpace()->GetNstate(); i++)	{
		if (GetCodonStateSpace()->Translation(i) == a)	{
			for (int j=0; j<GetCodonStateSpace()->GetNstate(); j++)	{
				int b = GetCodonStateSpace()->Translation(j);
				if (b != a)	{
					Q[a][b] += GetCodonSubMatrix()->Stationary(i) * (*GetCodonSubMatrix())(i,j);
				}
			}
		}
	}

	for (int b=0; b<GetNstate(); b++)	{
		if (b != a)	{
			Q[a][b] /= mStationary[a];
		}
	}

	double total = 0;
	for (int b=0; b<GetNstate(); b++)	{
		if (b != a)	{
			total += Q[a][b];
		}
	}
	Q[a][a] = -total;
}
