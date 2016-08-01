
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
using namespace std;

class SubMatrix  	{


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

	void 			CorruptMatrix();
	void 			UpdateMatrix();

	void			ActivatePowers();
	void			InactivatePowers();
	double 			Power(int n, int i, int j);
	double			GetUniformizationMu();

	double* 		GetEigenVal();
	double** 		GetEigenVect();
	double** 		GetInvEigenVect();


	void			ToStream(ostream& os);
	void			CheckReversibility();

	int 			GetDiagStat() { return ndiagfailed;}


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

	private:

	// v : eigenvalues
	// vi : imaginary part
	// u : the matrix of eigen vectors
	// invu : the inverse of u


	double ** u;
	double ** invu;
	double * v;
	double * vi;

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

#endif // SUBMATRIX_H
