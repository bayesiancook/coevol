
#ifndef MGMEANAAMUTSELCODONSUBMATRIX_H
#define MGMEANAAMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
#include "CovMatrix.h"
#include "Var.h"
#include "GaussIntegral.h"

static const int default_level = 2;

class MGMeanAAMutSelCodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGMeanAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* inmu, Var<CovMatrix>* insigma, int inlevel = default_level) :
			SubMatrix(instatespace->GetNstate()),
			MGCodonSubMatrix(instatespace, inNucMatrix), activated(false) {

		x = 0;
		ey = 0;
		w = 0;
		y = 0;
		mu = new double[Naa];
		emu = new double[Naa];
		mMutStationary = new double[GetNstate()];
		level = default_level;
		LoadGrid();
		SetMu(inmu);
		SetSigma(insigma);
		NegativeTransitionCounter = 0;
	}

	~MGMeanAAMutSelCodonSubMatrix()	{
		delete[] mu;
		delete[] mMutStationary;
	}

	// useful?
	/*
	double	GetMu(int aa1)	{return mu[aa1];}
	double 	GetSigma() {return sigma;}
	*/

	int GetNegativeTransitionNumber() {return NegativeTransitionCounter;}

	void Activate()	{
		activated = true;
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
			emu[i] = exp(mu[i] - max);
		}
	}

	void	SetSigma(Var<CovMatrix>* insigma)	{
		sigma = insigma;
		gridflag = false;
	}

	void ComputeFStat()	{

		for (int i=0; i<N; i++)	{
			double total = 0;
			for (int codon=0; codon<GetNstate(); codon++)	{
				fstat[i][codon] = mMutStationary[codon] * ey[i][GetCodonStateSpace()->Translation(codon)] * emu[GetCodonStateSpace()->Translation(codon)];
				total += fstat[i][codon];
			}
			for (int codon=0; codon<GetNstate(); codon++)	{
				fstat[i][codon] /= total;
			}
		}
	}
	// x is a vector of log fitnesses over the 20 amino acids
	/*
	double FStationary(double* ey, int refcodon)	{
	
		double total = 0;
		for (int i=0; i<GetNstate(); i++) {
			total += mMutStationary[i] * ey[GetCodonStateSpace()->Translation(i)] * emu[GetCodonStateSpace()->Translation(i)];
		}

		return mMutStationary[refcodon] * ey[GetCodonStateSpace()->Translation(refcodon)] * emu[GetCodonStateSpace()->Translation(refcodon)] / total;
	}
	*/

	double GetMeanStationary(int codon)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += w[i] * fstat[i][codon];
			// total += w[i] * FStationary(ey[i],codon);
		}
		return total / exp(0.5 * Naa * log(Pi));
	}

	double FTransition(int i, int refcodon1, int refcodon2)	{

		int aa1 = GetCodonStateSpace()->Translation(refcodon1);
		int aa2 = GetCodonStateSpace()->Translation(refcodon2);

		double pi = fstat[i][refcodon1];
		// double pi = FStationary(ey,refcodon1);
		double s = y[i][aa2] + mu[aa2] - y[i][aa1] - mu[aa1];
		double es = ey[i][aa2] * emu[aa2] / ey[i][aa1] / emu[aa1];
		double omega = 1;
		if (s < -200)	{
			omega = 0;
		}
		else if (s < -100)	{
			omega = -s * es;
		}
		else if (s>100)	{
			omega = s;
		}
		else if (fabs(s) > 1e-6)	{
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
			total += w[i] * FTransition(i,codon1, codon2);
		}
		return total / exp(0.5 * Naa * log(Pi));
	}

	void ComputeStationary()	{

		if (activated)	{
			if (! gridflag)	{
				Update();
			}
			ComputeFStat();
			// compute stationary probabilities under pure mutation model
			double total = 0;
			for (int i=0; i<GetNstate(); i++)	{
				mMutStationary[i] = NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i));
				total += mMutStationary[i];
			}

			// renormalize stationary probabilities
			for (int i=0; i<GetNstate(); i++)	{
				mMutStationary[i]  /= total;
			}

			total = 0;
			for (int i=0; i<GetNstate(); i++)	{
				mStationary[i] = GetMeanStationary(i);
				total += mStationary[i];
				// cerr << i << '\t' << mStationary[i] << '\n';
			}
			for (int i=0; i<GetNstate(); i++)	{
				mStationary[i] /= total;
			}
		}
	}

	void  ComputeArray(int i)	{

		if (activated)	{
		if (! gridflag)	{
			Update();
		}
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
						if (tmp < 0)	{
							tmp = 1e-100;
							NegativeTransitionCounter++;
						}
						Q[i][j] *= GetMeanTransition(i,j) / mStationary[i];
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
	}

	void LoadGrid()	{
		if (x)	{
			for (int i=0; i<N; i++)	{
				delete[] x[i];
				delete[] ey[i];
				delete[] y[i];
				delete[] fstat[i];
			}
			delete[] x;
			delete[] ey;
			delete[] y;
			delete[] w;
			delete[] fstat;
		}

		ostringstream xfile;
		xfile <<  "herm_d20_level" << level << "_x.txt";
		ifstream xis(xfile.str().c_str());

		ostringstream wfile;
		wfile <<  "herm_d20_level" << level << "_w.txt";
		ifstream wis(wfile.str().c_str());
	
		int Nstate;
		xis >> N >> Nstate;
		if (Nstate != Naa)	{
			cerr << "error in SparseGaussIntegral: non matching dimension\n";
			exit(1);
		}
		cerr << N << '\t' << Nstate << '\n';

		x = new double*[N];
		ey = new double*[N];
		y = new double*[N];
		fstat = new double*[N];
		w = new double[N];
		for (int i=0; i<N; i++)	{
			x[i] = new double[Naa];
			ey[i] = new double[Naa];
			y[i] = new double[Naa];
			fstat[i] = new double[GetNstate()];
			for (int k=0; k<Naa; k++)	{
				xis >> x[i][k];
			}
			wis >> w[i];
		}
	}


	void Update()	{
		if (activated)	{
			double* aux = new double[Naa];
			double** u = sigma->GetEigenVect();
			double* diag = sigma->GetEigenVal();
			for (int i=0; i<N; i++)	{
				for (int j=0; j<Naa; j++)	{
					aux[j] = x[i][j] * sqrt(diag[j]) * sqrt(2);
				}
				for (int j=0; j<Naa; j++)	{
					double total = 0;
					for (int k=0; k<Naa; k++)	{
						total += u[j][k] * aux[k];
					}
					y[i][j] = total;
					if (isnan(total))	{
						cerr << "nan in MGMean::Update \n";
						cerr << *sigma << '\n' << '\n';
						for (int j=0; j<Naa; j++)	{
							cerr << aux[j] << '\t' << diag[j] << '\t' << x[i][j] << '\n';
						}
						exit(1);
						cerr << '\n';
						for (int j=0; j<Naa; j++)	{
							for (int k=0; k<Naa; k++)	{
								cerr << u[j][k] << '\t';
							}
							cerr << '\n';
						}
						cerr << '\n';
						exit(1);
					}
				}
				double max = y[i][0];
				for (int k =1; k<Naa; k++)	{
					if (max < y[i][k])	{
						max = y[i][k];
					}
				}
				for (int k=0; k<Naa; k++)	{
					ey[i][k] = exp(y[i][k] - max);
				}
			}
				
			
			delete[] aux;
			gridflag = true;
		}
	}

	private:

	int level;
	int N;
	double** x;
	double** fstat;
	double** y;
	double** ey;
	double* w;

	// data members

	double* mMutStationary;
	int refcodon;
	int refcodon1;
	int refcodon2;
	double* mu;
	double* emu;
	Var<CovMatrix>* sigma;
	int NegativeTransitionCounter;

	bool gridflag;
	bool activated;
};


class RandomMGMeanAAMutSelCodonSubMatrix : public MGMeanAAMutSelCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGMeanAAMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inMu, Var<CovMatrix>* inSigma) :
			SubMatrix(instatespace->GetNstate()),
			MGMeanAAMutSelCodonSubMatrix(instatespace,inmatrix,inMu->GetArray(), inSigma) ,
			RandomCodonSubMatrix(instatespace),
			matrix(inmatrix), Mu(inMu), Sigma(inSigma) {

		cerr << "register\n";
		Register(matrix);
		Register(Mu);
		Register(Sigma);
		
	} 

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetMu(Mu->GetArray());
		SetSigma(Sigma);
	}

	protected:

	RandomSubMatrix* matrix;
	Var<RealVector>* Mu;
	Var<CovMatrix>* Sigma;
};

#endif
