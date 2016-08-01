#include "MGCodonAAProfileMutSelCodonSubMatrix.h"


void MGCodonAAProfileMutSelCodonSubMatrix::ComputeStationary()    {

	// compute stationary probabilities
	double total = 0;
	for (int i=0; i<GetNstate(); i++)       {
		mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) *
				NucMatrix->Stationary(GetCodonPosition(1,i)) *
                        	NucMatrix->Stationary(GetCodonPosition(2,i)) *
                        	exp(Neff->val() * log( (*CodonProfile)[i] )) *
                        	exp(Neff->val() * log( (*AAProfile)[GetCodonStateSpace()->Translation(i)] ));
		total += mStationary[i];
	}
	// renormalize stationary probabilities
        for (int i=0; i<GetNstate(); i++)       {
		mStationary[i]  /= total;
	}
}


void MGCodonAAProfileMutSelCodonSubMatrix::ComputeArray(int i)    {

        double total = 0;
        for (int j=0; j<GetNstate(); j++)       {
                if (i!=j)       {
                        int pos = GetDifferingPosition(i,j);
                        if ((pos != -1) && (pos != 3))  {
                                int a = GetCodonPosition(pos,i);
                                int b = GetCodonPosition(pos,j);
                                if (a == b)     {
                                        cerr << "identical states\n";
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
                                        cerr << pos << '\n';
                                        exit(1);
                                }
                                Q[i][j] = (*NucMatrix)(a,b);
                                if (! Synonymous(i,j))  {//When event is nonsynonymous, NeffDelta is a function of CodonProfile and AAProfile
                                        double NeffDeltaF = Neff->val() *
					(	log((*CodonProfile)[j] / (*CodonProfile)[i]) +
						log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)])
					);
                                        if (fabs(NeffDeltaF) < TOOSMALL)        {
                                                Q[i][j] /= ( 1.0 - (NeffDeltaF / 2) );
						//cerr << "TOOSMALL\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
                                        }
					else if (NeffDeltaF > TOOLARGE)	{
						Q[i][j] *= NeffDeltaF;
						//cerr << "TOOLARGE\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
					}
					else if (NeffDeltaF < TOOLARGENEGATIVE)	{
						Q[i][j] = 0;
						//cerr << "TOOLARGENEGATIVE\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
	
					}
                                        else    {
                                                Q[i][j] *=  (NeffDeltaF)/(1.0 - exp(-NeffDeltaF));
                                        }
                                }
				else	{// When event is synonymous, NeffDelta is only a function of CodonProfile
                                        double NeffDeltaF = Neff->val() * log( (*CodonProfile)[j] / (*CodonProfile)[i] );
                                        if (fabs(NeffDeltaF) < TOOSMALL)        {
                                                Q[i][j] /= ( 1.0 - (NeffDeltaF / 2) );
						//cerr << "TOOSMALL\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
                                        }
					else if (NeffDeltaF > TOOLARGE)	{
						Q[i][j] *= NeffDeltaF;
						//cerr << "TOOLARGE\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
					}
					else if (NeffDeltaF < TOOLARGENEGATIVE)	{
						Q[i][j] = 0;
						//cerr << "TOOLARGENEGATIVE\tNeff: " << Neff->val() << ", DeltaF: " << log((*AAProfile)[GetCodonStateSpace()->Translation(j)] / (*AAProfile)[GetCodonStateSpace()->Translation(i)]) << ", ";
						//cerr << "NeffDeltaF: " << NeffDeltaF << "\n";
						//cerr.flush();
	
					}
                                        else    {
                                                Q[i][j] *=  (NeffDeltaF)/(1.0 - exp(-NeffDeltaF));
                                        }

				}
                        }
                        else    {
                                Q[i][j] = 0;
                        }
                        total += Q[i][j];
                       
			if (Q[i][j] < 0)        {
                                cerr << "negative entry in matrix\n";
                                exit(1);
                        }
			if (isinf(Q[i][j]))	{
				cerr << "inf Q[i][j]\n";
				exit(1);
			}
			if (isnan(Q[i][j]))	{
				cerr << "nan Q[i][j]\n";
				exit(1);
			}

                }
        }
        Q[i][i] = -total;
        if (total <0)   {
                cerr << "negative rate away\n";
                exit(1);
        }
}

                                                                                                                                                                              
