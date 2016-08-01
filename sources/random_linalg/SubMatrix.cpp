#include "linalg.h"
#include "SubMatrix.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

int SubMatrix::nuni = 0;
int SubMatrix::nunimax = 0;

// ---------------------------------------------------------------------------
//		 SubMatrix()
// ---------------------------------------------------------------------------
/*
SubMatrix::SubMatrix() : Nstate(0), normalise(false) {}


SubMatrix::SubMatrix(const SubMatrix& from) : Nstate(from.inNstate), normalise(from.innormalise)	{
	Create();
}

SubMatrix& SubMatrix::operator=(const SubMatrix& from)	{

	if (Nstate != from.Nstate)	{
		cerr << "error in SubMatrix::operator= : non matching dimensions : " << Nstate << " and " << from.Nstate << '\n';
		throw;
	}
	// copy?
	// Corrupt ?
	CorruptMatrix();
}
*/
SubMatrix::SubMatrix(int inNstate, bool innormalise) : Nstate(inNstate), normalise(innormalise)	{
	ndiagfailed = 0;
	Create();
}

void SubMatrix::Create()	{

	Q = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		Q[i] = new double[Nstate];
	}

	u = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		u[i] = new double[Nstate];
	}

	invu = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		invu[i] = new double[Nstate];
	}

	v = new double[Nstate];
	vi = new double[Nstate];

	mStationary = new double[Nstate];

	UniMu = 1;
	mPow = new double**[UniSubNmax];
	for (int n=0; n<UniSubNmax; n++)	{
		mPow[n] = 0;
		/*
		mPow[n] = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			mPow[n][i] = new double[Nstate];
		}
		*/
	}

	flagarray = new bool[Nstate];
	diagflag = false;
	statflag = false;
	for (int i=0; i<Nstate; i++)	{
		flagarray[i] = false;
	}		
	powflag = false;

}

// ---------------------------------------------------------------------------
//		 ~SubMatrix()
// ---------------------------------------------------------------------------


SubMatrix::~SubMatrix()	{

	for (int i=0; i<Nstate; i++)	{
		delete[] Q[i];
		delete[] u[i];
		delete[] invu[i];
	}
	delete[] Q;
	delete[] u;
	delete[] invu;

	if (mPow)	{
		for (int n=0; n<UniSubNmax; n++)	{
			if (mPow[n])	{
				for (int i=0; i<Nstate; i++)	{
					delete[] mPow[n][i];
				}
				delete[] mPow[n];
			}
		}
		delete[] mPow;
	}
	delete[] mStationary;
	delete[] flagarray;
	delete[] v;
	delete[] vi;


}

// ---------------------------------------------------------------------------
//		 void ScalarMul()
// ---------------------------------------------------------------------------

void	SubMatrix::ScalarMul(double e)	{

	for (int i=0; i< Nstate ; i++)	{
		v[i] *= e;
		vi[i] *= e;
	}
	UniMu *= e;
}


// ---------------------------------------------------------------------------
//		 Diagonalise()
// ---------------------------------------------------------------------------

int SubMatrix::Diagonalise()	{

	if (! ArrayUpdated())	{
		UpdateMatrix();
	}

	int nmax = 1000;
	double epsilon = 1e-20;
	int n = LinAlg::DiagonalizeRateMatrix(Q,mStationary,Nstate,v,u,invu,nmax,epsilon);
	bool failed = (n == nmax);
	
	/*
	double * w = new double[Nstate];
	int* iw = new int[Nstate];
	double** a = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		a[i] = new double[Nstate];
	}

	// copy Q into a :
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			a[i][j] = Q[i][j];
		}
	}

	// diagonalise a into v and u
	int failed = EigenRealGeneral(Nstate, a, v, vi, u, iw, w);
	if (failed)	{
		ndiagfailed ++;
	}

	// copy u into a :
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			a[i][j] = u[i][j];
		}
	}
	// invert a into invu
	InvertMatrix(a, Nstate, w, iw, invu);

	for (int i=0; i<Nstate; i++)	{
		delete[] a[i];
	}
	delete[] a;
	delete[] w;
	delete[] iw;
	*/

	diagflag = true;

	return failed;
}


// ---------------------------------------------------------------------------
//		 ComputeRate()
// ---------------------------------------------------------------------------

double SubMatrix::GetRate()	{

	if (! ArrayUpdated())	{
		UpdateStationary();
		for (int k=0; k<Nstate; k++)	{
			ComputeArray(k);
		}
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	double norm = 0;
	for (int i=0; i<Nstate-1; i++)	{
		for (int j=i+1; j<Nstate; j++)	{
			norm += mStationary[i] * Q[i][j];
		}
	}
	return 2 * norm;
}


double* SubMatrix::GetEigenVal() {
	if (! diagflag)	{
		Diagonalise();
	}
	return v;
}

double** SubMatrix::GetEigenVect() {
	if (! diagflag)	{
		Diagonalise();
	}
	return u;
}

double** SubMatrix::GetInvEigenVect() {
	if (! diagflag) Diagonalise();
	return invu;
}

// ---------------------------------------------------------------------------
//		 Update()
// ---------------------------------------------------------------------------

void SubMatrix::UpdateMatrix()	{
	UpdateStationary();
	for (int k=0; k<Nstate; k++)	{
		ComputeArray(k);
	}
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = true;
	}
	if (isNormalised())	{
		Normalise();
	}
	// CheckReversibility();
}



// ---------------------------------------------------------------------------
//		 Normalise()
// ---------------------------------------------------------------------------

void SubMatrix::Normalise()	{

	double norm = GetRate();
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			Q[i][j] /= norm;
		}
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Powers
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void SubMatrix::ActivatePowers()	{

	if (! powflag)	{
		if (! ArrayUpdated())	{
			UpdateMatrix();
		}

		UniMu = 0;
		for (int i=0; i<Nstate; i++)	{
			if (UniMu < fabs(Q[i][i]))	{
				UniMu = fabs(Q[i][i]);
			}
		}

		CreatePowers(0);
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				mPow[0][i][j] = 0;
			}
		}
		for (int i=0; i<Nstate; i++)	{
			mPow[0][i][i] = 1;
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				mPow[0][i][j] += Q[i][j] / UniMu;
				if (mPow[0][i][j] < 0)	{
					cerr << "error in SubMatrix::ComputePowers: negative prob : ";
					cerr << i << '\t' << j << '\t' << mPow[0][i][j] << '\n';
					cerr << "Nstate : " << Nstate << '\n';
					exit(1);
				}
			}
		}
		npow = 1;
		powflag = true;
	}
}

void SubMatrix::InactivatePowers()	{

	if (powflag)	{
		for (int n=0; n<UniSubNmax; n++)	{
			if (mPow[n])	{
				for (int i=0; i<Nstate; i++)	{
					delete[] mPow[n][i];
				}
				delete [] mPow[n];
				mPow[n] = 0;
			}
		}
		nunimax += npow;
		nuni++;

		npow = 0;
		powflag = false;
	}
}

void SubMatrix::CreatePowers(int n)	{

	if (! mPow[n])	{
		mPow[n] = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			mPow[n][i] = new double[Nstate];
		}
	}
}


double SubMatrix::GetUniformizationMu() {

	if (! powflag)	{
		ActivatePowers();
	}
	return UniMu;
}

double SubMatrix::Power(int n, int i, int j)	{

	if (! powflag)	{
		ActivatePowers();
	}
	if (!n)	{
		return (i == j);
	}
	if (n > UniSubNmax)	{
		return Stationary(j);
	}
	if (n > npow)	{
		ComputePowers(n);
	}
	return mPow[n-1][i][j];
}


void SubMatrix::ComputePowers(int N)	{

	if (! powflag)	{
		ActivatePowers();
	}
	if (N>npow)	{
		for (int n=npow; n<N; n++)	{
			CreatePowers(n);
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					double& t = mPow[n][i][j];
					t = 0;
					for (int k=0; k<Nstate; k++)	{
						t += mPow[n-1][i][k] * mPow[0][k][j];
					}
				}
			}
		}
		npow = N;
	}
}

void SubMatrix::ToStream(ostream& os)	{

	os << GetNstate() << '\n';
	os << "stationaries: \n";
	for (int i=0; i<GetNstate(); i++)	{
		os << Stationary(i) << '\t';
	}
	os << '\n';
	
	os << "rate matrix\n";
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			os << Q[i][j] << '\t';
		}
		os << '\n';
	}

	os << '\n';
	for (int i=0; i<GetNstate(); i++)	{
		os << v[i] << '\t';
	}
	os << '\n';
	
}

void SubMatrix::CheckReversibility()	{

	double max = 0;
	int imax = 0;
	int jmax = 0;
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=i+1; j<GetNstate(); j++)	{
			double tmp = fabs(Stationary(i) * Q[i][j] - Stationary(j) * Q[j][i]);
			if (max < tmp)	{
				max = tmp;
				imax = i;
				jmax = j;
			}
		}
	}
	if (max > 1e-6)	{
		cerr << "max irreversibility: " << max << '\n';
		cerr << imax << '\t' << jmax << '\t' << Stationary(imax) << '\t' << Q[imax][jmax] << '\t' << Stationary(jmax) << '\t' << Q[jmax][imax] << '\n';
		exit(1);
	}	
}

