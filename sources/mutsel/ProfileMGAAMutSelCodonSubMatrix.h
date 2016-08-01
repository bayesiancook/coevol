
#ifndef PROFILEMGAAMUTSELCODONSUBMATRIX_H
#define PROFILEMGAAMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
#include "CovMatrix.h"
#include "Var.h"

inline double mutselomega(double s, double es = -1)	{
	double omega = 1;
	if (s < -100)	{
		omega = 0;
	}
	else if (s>100)	{
		omega = s;
	}
	else if (fabs(s) > 1e-12)	{
		if (es == -1)	{
			es = exp(s);
		}
		omega = s / (1 - 1.0 / es);
	}
	return omega;
}

class MGAAMutSelCodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix) : 
	// MGAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* inaalogfitness, double* inmu = 0, double in popsize = 1) : 
			SubMatrix(instatespace->GetNstate()),
			MGCodonSubMatrix(instatespace, inNucMatrix) {

		mu = new double[Naa];
		emu = new double[Naa];
		for (int k=0; k<Naa; k++)	{
			mu[k] = 0;
			emu[k] = 1;
		}
		aascaledlogfitness = new double[Naa];
		expaascaledlogfitness = new double[Naa];
		for (int k=0; k<Naa; k++)	{
			aascaledlogfitness[k] = 0;
			expaascaledlogfitness[k] = 1;
		}
		mMutStationary = new double[GetNstate()];
		popsize = 1;
		statnorm = 0;
	}

	~MGAAMutSelCodonSubMatrix()	{
		delete[] mu;
		delete[] emu;
		delete[] mMutStationary;
		delete[] aascaledlogfitness;
		delete[] expaascaledlogfitness;
	}

	protected:

	void SetMu(double* inmu)	{
		for (int i=0; i<Naa; i++)	{
			mu[i] = popsize * inmu[i];
		}
		double max = mu[0];
		for (int i=1; i<Naa; i++)	{
			if (max < mu[i])	{
				max = mu[i];
			}
		}
		for (int i=0; i<Naa; i++)	{
			emu[i] = exp(mu[i] - max);
		}
	}

	/*
	void SetMu(double* inemu)	{
		for (int i=0; i<Naa; i++)	{
			emu[i] = inemu[i];
			if (emu[i] > 1e-10)	{
				mu[i] = log(emu[i]);
			}
			else	{
				mu[i] = -30;
			}
			mu[i] *= popsize;
			emu[i] = exp(mu[i]);
		}
	}
	*/

	void SetAAScaledLogFitness(double* inaa)	{
		for (int j=0; j<Naa; j++)	{
			expaascaledlogfitness[j] = inaa[j];
			if (expaascaledlogfitness[j] < 1e-10)	{
				aascaledlogfitness[j] = -30;
			}
			else	{
				aascaledlogfitness[j] = log(expaascaledlogfitness[j]);
			}
			if (popsize != 1)	{
				aascaledlogfitness[j] *= popsize;
				expaascaledlogfitness[j] = exp(aascaledlogfitness[j]);
			}
		}
	}

	void SetPopSize(double inpopsize)	{
		popsize = inpopsize;
	}

	void ComputeStationaryNormalisationFactor()	{
		statnorm = 0;
		for (int k=0; k<GetNstate(); k++) {
			statnorm += mMutStationary[k] * expaascaledlogfitness[GetCodonStateSpace()->Translation(k)] * emu[GetCodonStateSpace()->Translation(k)];
		}
	}

	double FStationary(int refcodon)	{
		return mMutStationary[refcodon] * expaascaledlogfitness[GetCodonStateSpace()->Translation(refcodon)] * emu[GetCodonStateSpace()->Translation(refcodon)] / statnorm;
	}

	double FTransition(int refcodon1, int refcodon2)	{

		double* y = aascaledlogfitness;
		double* ey = expaascaledlogfitness;

		int aa1 = GetCodonStateSpace()->Translation(refcodon1);
		int aa2 = GetCodonStateSpace()->Translation(refcodon2);

		double s = y[aa2] + mu[aa2] - y[aa1] - mu[aa1];
		double es = ey[aa2] * emu[aa2] / ey[aa1] / emu[aa1];
		return mutselomega(s,es);
	}
		
	void ComputeStationary()	{

		for (int i=0; i<GetNstate(); i++)	{
			mMutStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
		}

		ComputeStationaryNormalisationFactor();

		double total = 0;
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] = FStationary(i);
			total += mStationary[i];
		}
		if (fabs(total - 1)>1e-8)	{
			cerr << "error in MGMix: total stationary not equal to unity: " << total << '\n';
		}
	}

	void  ComputeArray(int i)	{

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
						double tmp = FTransition(i,j);
						Q[i][j] *= tmp;
						if (tmp < 0)	{
							cerr << "negative transition\n";
							cerr << tmp << '\n';
							exit(1);
						}
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
			cerr << total << '\n';
			exit(1);
		}
	}

	// data members
	protected:

	double* mMutStationary;
	double* mu;
	double* emu;
	double* aascaledlogfitness;
	double* expaascaledlogfitness;
	double statnorm;
	double popsize;
};


class RandomMGAAMutSelCodonSubMatrix : public MGAAMutSelCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inAAScaledLogFitness, Var<RealVector>* inMu = 0, Var<Real>* inLogPopSize = 0, Var<PosReal>* inPopSize = 0) :
			SubMatrix(instatespace->GetNstate()),
			MGAAMutSelCodonSubMatrix(instatespace,inmatrix) ,
			RandomCodonSubMatrix(instatespace),
			matrix(inmatrix), AAScaledLogFitness(inAAScaledLogFitness), Mu(inMu), LogPopSize(inLogPopSize), PopSize(inPopSize) {

		Register(matrix);
		if (Mu)	{
			Register(Mu);
		}
		Register(AAScaledLogFitness);
		if (LogPopSize && PopSize)	{
			cerr << "error in RandomMGMix: cannot define both PopSize and LogPopSize\n";
			exit(1);
		}
		if (PopSize)	{
			Register(PopSize);
		}
		if (LogPopSize)	{
			Register(LogPopSize);
		}
	} 

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		if (PopSize)	{
			SetPopSize(PopSize->val());
		}
		if (LogPopSize)	{
			SetPopSize(exp(LogPopSize->val()));
		}
		if (Mu)	{
			SetMu(Mu->GetArray());
		}
		SetAAScaledLogFitness(AAScaledLogFitness->GetArray());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<Profile>* AAScaledLogFitness;
	Var<RealVector>* Mu;
	Var<Real>* LogPopSize;
	Var<PosReal>* PopSize;
};


/*
class MGBGCAAMutSelCodonSubMatrix : public MGAAMutSelCodonSubMatrix	{

	public:

	MGBGCAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix) : 
			SubMatrix(instatespace->GetNstate()),
			MGAAMutSelCodonSubMatrix(instatespace, inNucMatrix)	{
		bgc = 0;
		ebgc = 1;
	}

	protected:

	void SetBGC(double inbgc)	{
		bgc = popsize * inbgc;
		ebgc = exp(bgc);
	}

	void ComputeStationaryNormalisationFactor()	{
		statnorm = 0;
		for (int k=0; k<GetNstate(); k++) {
			double bgcfactor = 1;
			for (int pos=0; pos<3; pos++)	{
				int a = GetCodonPosition(pos,k);
				if ((a == 1) || (a == 2))	{
					bgcfactor *= ebgc;
				}
			}
			statnorm += mMutStationary[k] * expaascaledlogfitness[GetCodonStateSpace()->Translation(k)] * emu[GetCodonStateSpace()->Translation(k)] * bgcfactor;
		}
	}

	double FStationary(int refcodon)	{
		double bgcfactor = 1;
		for (int pos=0; pos<3; pos++)	{
			int a = GetCodonPosition(pos,refcodon);
			if ((a == 1) || (a == 2))	{
				bgcfactor *= ebgc;
			}
		}
		return mMutStationary[refcodon] * expaascaledlogfitness[GetCodonStateSpace()->Translation(refcodon)] * emu[GetCodonStateSpace()->Translation(refcodon)] * bgcfactor / statnorm;
	}

	double FTransition(int refcodon1, int refcodon2, double bgcdiff, double ebgcdiff)	{

		double* y = aascaledlogfitness;
		double* ey = expaascaledlogfitness;

		int aa1 = GetCodonStateSpace()->Translation(refcodon1);
		int aa2 = GetCodonStateSpace()->Translation(refcodon2);

		double s = y[aa2] + mu[aa2] - y[aa1] - mu[aa1] + bgcdiff;
		double es = ey[aa2] * emu[aa2] / ey[aa1] / emu[aa1] * ebgcdiff;
		return mutselomega(s,es);
	}
		
	void  ComputeArray(int i)	{

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
					bool initialgc = ((a == 1) || (a == 2));
					bool finalgc = ((b == 1) || (b == 2));
					double bgcdiff = 0;
					if (initialgc)	{
						bgcdiff -= bgc;
					}
					if (finalgc)	{
						bgcdiff += bgc;
					}
					Q[i][j] = (*NucMatrix)(a,b);
					if (! Synonymous(i,j))	{
						double tmp = FTransition(i,j,bgcdiff,exp(bgcdiff));
						Q[i][j] *= tmp;
						if (tmp < 0)	{
							cerr << "negative transition\n";
							cerr << tmp << '\n';
							exit(1);
						}
					}
					else	{
						Q[i][j] *= mutselomega(bgcdiff,exp(bgcdiff));
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
			cerr << total << '\n';
			exit(1);
		}
	}

	private:
	double bgc;
	double ebgc;

};

class RandomMGBGCAAMutSelCodonSubMatrix : public MGBGCAAMutSelCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGBGCAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inAAScaledLogFitness, Var<Real>* inBGC, Var<RealVector>* inMu = 0, Var<Real>* inLogPopSize = 0, Var<PosReal>* inPopSize = 0) :
			SubMatrix(instatespace->GetNstate()),
			MGBGCAAMutSelCodonSubMatrix(instatespace,inmatrix) ,
			RandomCodonSubMatrix(instatespace),
			matrix(inmatrix), AAScaledLogFitness(inAAScaledLogFitness), BGC(inBGC), Mu(inMu), LogPopSize(inLogPopSize), PopSize(inPopSize) {

		Register(matrix);
		if (Mu)	{
			Register(Mu);
		}
		Register(AAScaledLogFitness);
		Register(BGC);
		if (LogPopSize && PopSize)	{
			cerr << "error in RandomMGMix: cannot define both PopSize and LogPopSize\n";
			exit(1);
		}
		if (PopSize)	{
			Register(PopSize);
		}
		if (LogPopSize)	{
			Register(LogPopSize);
		}
	} 

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		if (PopSize)	{
			SetPopSize(PopSize->val());
		}
		if (LogPopSize)	{
			SetPopSize(exp(LogPopSize->val()));
		}
		SetBGC(BGC->val());
		if (Mu)	{
			SetMu(Mu->GetArray());
		}
		SetAAScaledLogFitness(AAScaledLogFitness->GetArray());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<RealVector>* AAScaledLogFitness;
	Var<Real>* BGC;
	Var<RealVector>* Mu;
	Var<Real>* LogPopSize;
	Var<PosReal>* PopSize;
};
*/
#endif
