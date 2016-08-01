
#ifndef MGMIXAAMUTSELCODONSUBMATRIX_H
#define MGMIXAAMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
#include "CovMatrix.h"
#include "Var.h"

class MGMixAAMutSelCodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGMixAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, int inN) : 
			SubMatrix(instatespace->GetNstate()),
			MGCodonSubMatrix(instatespace, inNucMatrix) {

		N = inN;
		mu = new double[Naa];
		emu = new double[Naa];
		for (int k=0; k<Naa; k++)	{
			mu[k] = 0;
			emu[k] = 1;
		}
		weight = new double[N];
		statnorm = new double[N];
		for (int i=0; i<N; i++)	{
			weight[i] = 1.0 / N;
		}
		mixture = new double*[N];
		emixture = new double*[N];
		for (int i=0; i<N; i++)	{
			mixture[i] = new double[Naa];
			emixture[i] = new double[Naa];
			for (int k=0; k<Naa; k++)	{
				mixture[i][k] = 0;
				emixture[i][k] = 1;
			}
		}
		mMutStationary = new double[GetNstate()];
		popsize = 1;
	}

	~MGMixAAMutSelCodonSubMatrix()	{
		delete[] mu;
		delete[] emu;
		delete[] mMutStationary;
		delete[] weight;
		delete[] statnorm;
		delete[] mixture;
		delete[] emixture;
	}

	protected:

	void SetMu(double* inmu)	{
		for (int i=0; i<Naa; i++)	{
			mu[i] = inmu[i];
		}
		double max = mu[0];
		for (int i=1; i<Naa; i++)	{
			if (max < mu[i])	{
				max = mu[i];
			}
		}
		for (int i=0; i<Naa; i++)	{
			emu[i] = exp(popsize * (mu[i] - max));
		}
	}

	void	SetWeights(const double* inw)	{
		for (int i=0; i<N; i++)	{
			weight[i] = inw[i];
		}
	}

	void SetMixture(double* inmixture, int i)	{
		for (int j=0; j<Naa; j++)	{
			mixture[i][j] = inmixture[j];
		}	
		double max = mixture[i][0];
		for (int j=1; j<Naa; j++)	{
			if (max < mixture[i][j])	{
				max = mixture[i][j];
			}
		}
		for (int j=0; j<Naa; j++)	{
			emixture[i][j] = exp(popsize * (mixture[i][j] - max));
		}
	}

	void SetPopSize(double inpopsize)	{
		popsize = inpopsize;
	}

	void ComputeStationaryNormalisationFactor()	{
		for (int i=0; i<N; i++)	{
			statnorm[i] = 0;
			for (int k=0; k<GetNstate(); k++) {
				statnorm[i] += mMutStationary[k] * emixture[i][GetCodonStateSpace()->Translation(k)] * emu[GetCodonStateSpace()->Translation(k)];
			}
		}
	}

	double FStationary(int i, int refcodon)	{
		return mMutStationary[refcodon] * emixture[i][GetCodonStateSpace()->Translation(refcodon)] * emu[GetCodonStateSpace()->Translation(refcodon)] / statnorm[i];
	}

	double GetMeanStationary(int codon)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += weight[i] * FStationary(i,codon);
		}
		return total;
	}

	double FTransition(int i, int refcodon1, int refcodon2)	{

		double* y = mixture[i];
		double* ey = emixture[i];

		int aa1 = GetCodonStateSpace()->Translation(refcodon1);
		int aa2 = GetCodonStateSpace()->Translation(refcodon2);

		double pi = FStationary(i,refcodon1);
		double s = popsize * (y[aa2] + mu[aa2] - y[aa1] - mu[aa1]);
		double es = ey[aa2] * emu[aa2] / ey[aa1] / emu[aa1];
		double omega = 1;
		if (s < -100)	{
			// cerr << "<";
			omega = 0;
		}
		/*
		else if (s < -100)	{
			omega = -s * es;
		}
		*/
		else if (s>100)	{
			// cerr << ">";
			omega = s;
		}
		else if (fabs(s) > 1e-12)	{
			omega = s * es / (es - 1);
		}
		if (omega < 0)	{
			cerr << "negative omega : " << omega << '\n';
			exit(1);
		}
		return pi * omega;
	}
		
	double GetMeanTransition(int codon1, int codon2)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += weight[i] * FTransition(i,codon1, codon2);
		}
		return total;
	}

	void ComputeStationary()	{

		for (int i=0; i<GetNstate(); i++)	{
			mMutStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
		}

		ComputeStationaryNormalisationFactor();

		double total = 0;
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] = GetMeanStationary(i);
			total += mStationary[i];
		}
		/*
		if (fabs(total - 1)>1e-8)	{
			cerr << "error in MGMix: total stationary not equal to unity: " << total << '\n';
		}
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] /= total;
		}
		*/
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
						double tmp = GetMeanTransition(i,j);
						Q[i][j] *= tmp / mStationary[i];
						if (tmp < 0)	{
							cerr << "negative transition\n";
							cerr << tmp << '\n';
							exit(1);
						}
						if (mStationary[i] < 0)	{
							cerr << "negative stat\n";
							cerr << mStationary[i] << '\n';
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

	int N;
	double* mMutStationary;
	double* mu;
	double* emu;
	double* weight;
	double** mixture;
	double** emixture;
	double* statnorm;
	
	double popsize;
};


class RandomMGMixAAMutSelCodonSubMatrix : public MGMixAAMutSelCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGMixAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inMu, VarArray<RealVector>* inMixture, Var<Profile>* inWeight, Var<Real>* inLogPopSize = 0, Var<PosReal>* inPopSize = 0) :
			SubMatrix(instatespace->GetNstate()),
			MGMixAAMutSelCodonSubMatrix(instatespace,inmatrix,inMixture->GetSize()) ,
			RandomCodonSubMatrix(instatespace),
			matrix(inmatrix), Mu(inMu), Mixture(inMixture), Weight(inWeight), LogPopSize(inLogPopSize), PopSize(inPopSize) {

		Register(matrix);
		Register(Mu);
		Register(Weight);
		for (int i=0; i<Mixture->GetSize(); i++)	{
			Register(Mixture->GetVal(i));
		}
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
		SetMu(Mu->GetArray());
		SetWeights(Weight->GetArray());
		for (int i=0; i<Mixture->GetSize(); i++)	{
			SetMixture(Mixture->GetVal(i)->GetArray(), i);
		}
	}

	protected:

	RandomSubMatrix* matrix;
	Var<RealVector>* Mu;
	VarArray<RealVector>* Mixture;
	Var<Profile>* Weight;
	Var<Real>* LogPopSize;
	Var<PosReal>* PopSize;
};


class MGMixBGCAAMutSelCodonSubMatrix : public MGMixAAMutSelCodonSubMatrix	{

	public:

	MGMixBGCAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, int inN) : 
			SubMatrix(instatespace->GetNstate()),
			MGMixAAMutSelCodonSubMatrix(instatespace, inNucMatrix, inN)	{
		bgc = 0;
		ebgc = 1;
	}

	protected:

	void SetBGC(double inbgc)	{
		bgc = inbgc;
		ebgc = exp(inbgc);
	}

	void ComputeStationaryNormalisationFactor()	{
		for (int i=0; i<N; i++)	{
			statnorm[i] = 0;
			for (int k=0; k<GetNstate(); k++) {
				double bgcfactor = 1;
				for (int pos=0; pos<3; pos++)	{
					int a = GetCodonPosition(pos,k);
					if ((a == 1) || (a == 2))	{
						bgcfactor *= ebgc;
					}
				}
				statnorm[i] += mMutStationary[k] * emixture[i][GetCodonStateSpace()->Translation(k)] * emu[GetCodonStateSpace()->Translation(k)] * bgcfactor;
			}
		}
	}


	double FStationary(int i, int refcodon)	{
		double total = 0;
		for (int k=0; k<GetNstate(); k++) {
			double bgcfactor = 1;
			for (int pos=0; pos<3; pos++)	{
				int a = GetCodonPosition(pos,k);
				if ((a == 1) || (a == 2))	{
					bgcfactor *= ebgc;
				}
			}
			total += mMutStationary[k] * emixture[i][GetCodonStateSpace()->Translation(k)] * emu[GetCodonStateSpace()->Translation(k)] * bgcfactor;
		}
		double bgcfactor = 1;
		for (int pos=0; pos<3; pos++)	{
			int a = GetCodonPosition(pos,refcodon);
			if ((a == 1) || (a == 2))	{
				bgcfactor *= ebgc;
			}
		}
		return mMutStationary[refcodon] * emixture[i][GetCodonStateSpace()->Translation(refcodon)] * emu[GetCodonStateSpace()->Translation(refcodon)] * bgcfactor / total;
	}

	double GetMeanStationary(int codon)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += weight[i] * FStationary(i,codon);
		}
		return total;
	}

	double FTransition(int i, int refcodon1, int refcodon2, double bgcdiff,double ebgcdiff)	{

		double* y = mixture[i];
		double* ey = emixture[i];

		int aa1 = GetCodonStateSpace()->Translation(refcodon1);
		int aa2 = GetCodonStateSpace()->Translation(refcodon2);

		double pi = FStationary(i,refcodon1);
		double s = y[aa2] + mu[aa2] - y[aa1] - mu[aa1] + bgcdiff;
		double es = ey[aa2] * emu[aa2] / ey[aa1] / emu[aa1] * ebgcdiff;
		double omega = 1;
		if (s < -100)	{
			omega = 0;
		}
		else if (s>100)	{
			omega = s;
		}
		else if (fabs(s) > 1e-12)	{
			omega = s * es / (es - 1);
		}
		if (omega < 0)	{
			cerr << "negative omega : " << omega << '\n';
			exit(1);
		}
		return pi * omega;
	}
		
	double GetMeanTransition(int codon1, int codon2, double bgcdiff, double ebgcdiff)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += weight[i] * FTransition(i,codon1, codon2, bgcdiff,ebgcdiff);
		}
		return total;
	}

	void ComputeStationary()	{

		// compute stationary probabilities under pure mutation model
		for (int i=0; i<GetNstate(); i++)	{
			mMutStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
		}

		// ComputeStationaryNormalisationFactor();

		double total = 0;
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] = GetMeanStationary(i);
			total += mStationary[i];
		}
		/*
		if (fabs(total - 1)>1e-8)	{
			cerr << "error in MGMix: total stationary not equal to unity: " << total << '\n';
		}
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] /= total;
		}
		*/
	}

	double Synomega(double s)	{
		double omega = 1;
		if (s < -100)	{
			omega = 0;
		}
		else if (s>100)	{
			omega = s;
		}
		else if (fabs(s) > 1e-12)	{
			omega = s / (1 - exp(-s));
		}
		return omega;
	}

	void  ComputeArray(int i)	{

		double total = 0;
		for (int j=0; j<GetNstate(); j++)	{
			if (i!=j)	{
				int pos = GetDifferingPosition(i,j);
				if ((pos != -1) && (pos != 3))	{
					int a = GetCodonPosition(pos,i);
					int b = GetCodonPosition(pos,j);
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
						double tmp = GetMeanTransition(i,j,bgcdiff,exp(bgcdiff));
						Q[i][j] *= tmp / mStationary[i];
						if (tmp < 0)	{
							cerr << "negative transition\n";
							cerr << tmp << '\n';
							exit(1);
						}
						if (mStationary[i] < 0)	{
							cerr << "negative stat\n";
							cerr << mStationary[i] << '\n';
							exit(1);
						}
					}
					else	{
						Q[i][j] *= Synomega(bgcdiff);
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


class RandomMGMixBGCAAMutSelCodonSubMatrix : public MGMixBGCAAMutSelCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGMixBGCAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inMu, VarArray<RealVector>* inMixture, Var<Profile>* inWeight, Var<Real>* inBGC) :
			SubMatrix(instatespace->GetNstate()),
			MGMixBGCAAMutSelCodonSubMatrix(instatespace,inmatrix,inMixture->GetSize()) ,
			RandomCodonSubMatrix(instatespace),
			matrix(inmatrix), Mu(inMu), Mixture(inMixture), Weight(inWeight), BGC(inBGC) {

		Register(matrix);
		Register(Mu);
		Register(Weight);
		Register(BGC);
		for (int i=0; i<Mixture->GetSize(); i++)	{
			Register(Mixture->GetVal(i));
		}
	} 

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetMu(Mu->GetArray());
		SetWeights(Weight->GetArray());
		SetBGC(BGC->val());
		for (int i=0; i<Mixture->GetSize(); i++)	{
			SetMixture(Mixture->GetVal(i)->GetArray(), i);
		}
	}

	protected:

	RandomSubMatrix* matrix;
	Var<RealVector>* Mu;
	VarArray<RealVector>* Mixture;
	Var<Profile>* Weight;
	Var<Real>* BGC;
};

#endif
