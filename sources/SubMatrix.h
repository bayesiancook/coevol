
// SubMatrix:
// this class implements
// the instant rate matrix of a substitution process
// but only the mathematical aspects of it
//
// RandomSubMatrix:
// this class derives from SubMatrix and from Dnode
// therefore, it knows everything about substitution processes (as a SubMatrix)
// and at the same time, it can be inserted as a deterministic node in a probabilistic model (as a Dnode)
//
// If you need to define a new substitution process
// - derive a new class deriving (directly or indirectly) from SubMatrix
// in this class, implements the ComputeArray and ComputeStationary functions
// (these are the functions that construct the instant rates and the stationary probabilities of the matrix)
//
// - derive a new class deriving both from RandomSubMatrix, and from your new class
// in this class, implement SetParameters
// (this function is responsible for updating the internal parameters that the SubMatrix uses in ComputeArray And ComputeStationary,
// based on the values stored by the parent nodes)
//
//

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

class AbstractTransitionMatrix	{

	public:

	virtual ~AbstractTransitionMatrix() {}

	virtual void 		BackwardPropagate(const double* down, double* up, double length) = 0;
	virtual void 		ForwardPropagate(const double* up, double* down, double length) = 0;
	virtual const double*	GetStationary() = 0;
	virtual double		Stationary(int i) = 0;

	/*
	virtual int		DrawStationary() = 0;
	virtual int		DrawFiniteTime(int state) = 0;
	virtual void		GetFiniteTimeTransitionProb(int state, double* aux) = 0;
	*/

	virtual int		GetNstate() = 0;
	virtual void		CorruptMatrix() = 0;
	virtual double		operator()(int, int) = 0;
	virtual const double*	GetRow(int i) = 0;

	virtual bool check() {return true;}

	
};


class TransitionMatrix : public virtual AbstractTransitionMatrix	{

	public:

	TransitionMatrix(int inNstate)	{
		Nstate = inNstate;
		matflag = false;
		Create();
	}

	TransitionMatrix(int inNstate, double** inR, double *instationary)	{
		Nstate = inNstate;
		matflag = false;
		R = inR;
		stationary = instationary;
	}

	virtual ~TransitionMatrix()	{
		Delete();
	}

	int GetNstate() {return Nstate;}

	double** GetR()	{
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		return R;
	}

	const double* GetStationary()	{
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		return stationary;
	}

	const double* GetRow(int i)	{
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		return R[i];
	}

	double Stationary(int i)	{
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		return stationary[i];
	}

	double operator()(int i, int j)	{
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		return R[i][j];
	}

	virtual void 		BackwardPropagate(const double* down, double* up, double length) {
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		for (int i=0; i<Nstate; i++)	{
			double tmp = 0;
			for (int j=0; j<Nstate; j++)	{
				tmp += R[i][j] * down[j];
			}
			up[i] = tmp;
		}
		up[Nstate] = down[Nstate];

	}

	virtual void 		ForwardPropagate(const double* up, double* down, double length) {
		cerr << "in forward\n";
		exit(1);
		if (! matflag)	{
			ComputeArrayAndStat();
		}
		for (int i=0; i<Nstate; i++)	{
			double tmp = 0;
			for (int j=0; j<Nstate; j++)	{
				tmp += up[j] * R[j][i];
			}
			down[i] = tmp;
		}
	}

	void CorruptMatrix()	{
		matflag = false;
	}

	protected:

	void ComputeArrayAndStat()	{
		ComputeStationary();
		ComputeArray();
		matflag = true;
	}

	virtual void ComputeArray() = 0;
	virtual void ComputeStationary() = 0;

	virtual bool check() {
		bool ok = true;
		for(int i=0; i<Nstate; i++) {
			double total =0;
			for(int j=0; j<Nstate; j++) {
				if(R[i][j]<0) {
					//cerr << "Error : negative value in transition matrix : R[" << i << "][" << j << "] = " << R[i][j] << endl;
					ok = false;
				}
				total +=R[i][j];
			}
			if(total-1>10e-6 || total-1<-10e-6) {
				//cerr << "Error : a line in transition matrix does not sum to 1 : Line " << i << ", sum " << total << endl;
				ok = false;
			}
		}
		return ok;
	}

	void Create()	{
		stationary = new double[Nstate];
		// bkstationary = new double[Nstate];
		R = new double*[Nstate];
		// bkR = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			R[i] = new double[Nstate];
			// bkR[i] = new double[Nstate];
		}
	}

	void Delete()	{
		delete[] stationary;
		// delete[] bkstationary;
		for (int i=0; i<Nstate; i++)	{
			delete[] R[i];
			// delete[] bkR[i];
		}
		delete[] R;
		// delete[] bkR;
	}

	int Nstate;
	double* stationary;
	double** R;
	bool matflag;

};

class SubMatrix   : public virtual AbstractTransitionMatrix {


	protected:

	// these 2 pure virtual functions are the most essential component of the SubMatrix class
	// see GTRSubMatrix.cpp and CodonSubMatrix.cpp for examples

	// ComputeArray(int state) is in charge of computing the row of the rate matrix
	// corresponding to all possible rates of substitution AWAY from state
	//
	virtual void 		ComputeArray(int state) = 0;

	// ComputeStationary() is in charge of computing the vector of stationary probabilities (equilibirum frequencies)
	// of the substitution process
	virtual void 		ComputeStationary() = 0;

	public:

	static const int	UniSubNmax = 500;

	static int		nuni;
	static int		nunimax;

	static double		GetMeanUni() {return ((double) nunimax) / nuni;}

				SubMatrix(int Nstate, bool innormalise = false);
	virtual 		~SubMatrix();

	void			Create();

	double 			operator()(int, int);
	const double* 		GetRow(int i);

	virtual const double* 	GetStationary();
	double 			Stationary(int i);

	int 			GetNstate() {return Nstate;}

	double 			GetRate();
	void 			ScalarMul(double e);

	bool			isNormalised() {return normalise;}
	void			Normalise();

	virtual void 		CorruptMatrix();
	void 			UpdateMatrix();

	void			ActivatePowers();
	void			InactivatePowers();
	double 			Power(int n, int i, int j);
	double			GetUniformizationMu();

	double* 		GetEigenVal();
	double** 		GetEigenVect();
	double** 		GetInvEigenVect();


	virtual void		ToStream(ostream& os);
	void			CheckReversibility();

	int 			GetDiagStat() { return ndiagfailed;}

	virtual void 		BackwardPropagate(const double* down, double* up, double length);
	virtual void 		ForwardPropagate(const double* up, double* down, double length);
	// virtual void 		FiniteTime(int i0, double* down, double length);

	double**		GetQ() {return Q;}
	void			ComputeExponential(double range, double** expo);
	void			ApproachExponential(double range, double** expo, int prec = 1024);
	void			PowerOf2(double** y, int z);

	static double meanz;
	static double maxz;
	static double nz;


	protected:

	void 			UpdateRow(int state);
	void 			UpdateStationary();


	void 			ComputePowers(int n);
	void 			CreatePowers(int n);

	bool			ArrayUpdated();

	int 			Diagonalise();

	// data members

	bool powflag;
	bool diagflag;
	bool statflag;
	bool* flagarray;

	int Nstate;
	int npow;
	double UniMu;

	double*** mPow;

	// Q : the infinitesimal generator matrix
	double ** Q;

	// the stationary probabilities of the matrix
	double* mStationary;

	bool normalise;

	// an auxiliary matrix
	double** aux;

	protected:

	// v : eigenvalues
	// vi : imaginary part
	// u : the matrix of eigen vectors
	// invu : the inverse of u


	double ** u;
	double ** invu;
	double * v;
	double * vi;

	/*
	double ** expu;
	double ** expu2;
	int discn;
	void ComputeExponential(double length);
	*/

	int ndiagfailed;

	
};


//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------

inline double	SubMatrix::operator()(int i, int j)	{
	if (! flagarray[i])	{
		UpdateRow(i);
	}
	return Q[i][j];
}

inline const double* SubMatrix::GetRow(int i) {
	if (! flagarray[i])	{
		UpdateRow(i);
	}
	return Q[i];
}

inline const double* SubMatrix::GetStationary() {
	if (! statflag)	{
		UpdateStationary();
	}
	return mStationary;
}

inline double SubMatrix::Stationary(int i)	{
	if (! statflag)	{
		UpdateStationary();
	}
	return mStationary[i];
}


inline void SubMatrix::CorruptMatrix()	{
	diagflag = false;
	statflag = false;
	for (int k=0; k<Nstate; k++)	{
		flagarray[k] = false;
	}
	InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated()	{

	bool qflag = true;
	for (int k=0; k<Nstate; k++)	{
		qflag &= flagarray[k];
	}
	return qflag;
}

inline void SubMatrix::UpdateStationary()	{
	ComputeStationary();
	statflag = true;
}

inline void SubMatrix::UpdateRow(int state)	{

	if (isNormalised())	{
		UpdateMatrix();
	}
	else	{
		if (! statflag)	{
			UpdateStationary();
		}
		ComputeArray(state);
		flagarray[state] = true;
	}
}

inline void SubMatrix::BackwardPropagate(const double* up, double* down, double length)	{

	double** eigenvect = GetEigenVect();
	double** inveigenvect= GetInvEigenVect();
	double* eigenval = GetEigenVal();

	double* aux = new double[GetNstate()];

	for (int i=0; i<GetNstate(); i++)	{
		aux[i] = 0;
	}
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			aux[i] += inveigenvect[i][j] * up[j];
		}
	}

	for (int i=0; i<GetNstate(); i++)	{
		aux[i] *= exp(length * eigenval[i]);
	}

	for (int i=0; i<GetNstate(); i++)	{
		down[i] = 0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			down[i] += eigenvect[i][j] * aux[j];
		}
	}


	for (int i=0; i<GetNstate(); i++)	{
		if (std::isnan(down[i]))	{
			cerr << "error in back prop\n";
			for (int j=0; j<GetNstate(); j++)	{
				cerr << up[j] << '\t' << down[j] << '\t' << Stationary(j) << '\n';
			}
			exit(1);
		}
	}
	double maxup = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (up[k] <0)	{
			cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
			// down[k] = 0;
		}
		if (maxup < up[k])	{
			maxup = up[k];
		}
	}
	double max = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (down[k] <0)	{
			down[k] = 0;
		}
		if (max < down[k])	{
			max = down[k];
		}
	}
	if (maxup == 0)	{
		cerr << "error in backward propagate: null up array\n";
		exit(1);
	}
	if (max == 0)	{
		cerr << "error in backward propagate: null array\n";
		for (int k=0; k<GetNstate(); k++)	{
			cerr << up[k] << '\t' << down[k] << '\n';
		}
		cerr << '\n';
		exit(1);
	}
	down[GetNstate()] = up[GetNstate()];

	delete[] aux;
}

inline void SubMatrix::ForwardPropagate(const double* down, double* up, double length)	{

	double** eigenvect = GetEigenVect();
	double** inveigenvect= GetInvEigenVect();
	double* eigenval = GetEigenVal();

	double* aux = new double[GetNstate()];

	for (int i=0; i<GetNstate(); i++)	{
		aux[i] = 0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			aux[i] += down[j] * eigenvect[j][i];
		}
	}

	for (int i=0; i<GetNstate(); i++)	{
		aux[i] *= exp(length * eigenval[i]);
	}

	for (int i=0; i<GetNstate(); i++)	{
		up[i] = 0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			up[i] += aux[j] * inveigenvect[j][i];
		}
	}

	delete[] aux;
}


/*
inline void SubMatrix::ComputeExponential(double length)	{

	for (int i=0; i<GetNstate(); i++)	{
		if (! flagarray[i])	{
			UpdateRow(i);
		}
	}

	double t = length;
	for (int i=0; i<discn; i++)	{
		t /= 4.0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			expu[i][j] = t * Q[i][j];
		}
	}
	for (int i=0; i<GetNstate(); i++)	{
		expu[i][i] += 1.0;
	}

	for (int n=0; n<discn; n++)	{
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu[i][k] * expu[k][j];
				}
				expu2[i][j] = tmp;
			}
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu2[i][k] * expu2[k][j];
				}
				expu[i][j] = tmp;
			}
		}
	}
}

inline void SubMatrix::BackwardPropagate(const double* up, double* down, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		down[i] = 0;
	}
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			down[i] += expu[i][j] * up[j];
		}
	}

	for (int i=0; i<GetNstate(); i++)	{
		if (isnan(down[i]))	{
			cerr << "error in back prop\n";
			for (int j=0; j<GetNstate(); j++)	{
				cerr << up[j] << '\t' << down[j] << '\n';
			}
			exit(1);
		}
	}
	double maxup = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (up[k] <0)	{
			cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
			// down[k] = 0;
		}
		if (maxup < up[k])	{
			maxup = up[k];
		}
	}
	double max = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (down[k] <0)	{
			down[k] = 0;
		}
		if (max < down[k])	{
			max = down[k];
		}
	}
	if (maxup == 0)	{
		cerr << "error in backward propagate: null up array\n";
		exit(1);
	}
	if (max == 0)	{
		cerr << "error in backward propagate: null array\n";
		for (int k=0; k<GetNstate(); k++)	{
			cerr << up[k] << '\t' << down[k] << '\n';
		}
		cerr << '\n';
		exit(1);
	}
	down[GetNstate()] = up[GetNstate()];
}

inline void SubMatrix::FiniteTime(int i0, double* up, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		up[i] = expu[i0][i];
	}
}

inline void SubMatrix::ForwardPropagate(const double* down, double* up, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		up[i] = 0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			up[i] += down[j] * expu[j][i];
		}
	}
}
*/

#endif // SUBMATRIX_H
