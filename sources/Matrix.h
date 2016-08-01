#ifndef MATRIX_H
#define MATRIX_H


#include "BaseType.h"
#include "Var.h"

class RealMatrix : public BaseType	{


	RealMatrix() : rank(0), mat(0) {}

	RealMatrix(int inrank)	: rank(inrank) {
		Allocate();
	}

	~RealMatrix()	{
		Deallocate();
	}



	private:

	void Allocate()	{
		mat = new double*[rank];
		for (int i=0; i<rank; i++)	{
			mat[i] = new double[rank];
		}
	}

	void Deallocate()	{
		if (mat)	{
			for (int i=0; i<rank; i++)	{
				delete[] mat[i];
			}
			delete[] mat;
		}
	}

	double** mat;

};

class InverseWishart: public virtual Rvar<void*> {

	InverseWishart(Var<PosRealVector>* indiagpsi, int inweight);

	InverseWishart(int rank, Var<PosReal>** indiagpsi, int inweight);

	protected:

	int weight;
	Var<PosRealVector>* diagpsi;

};

class ConjugateInverseWishart : public InverseWishart, public ConjugatePrior<void*>	{

	private:

	// sufficient statistics: an empirical covariance matrix

};



#endif
