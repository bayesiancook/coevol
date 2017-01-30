

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

#ifndef RANDOMSUBMATRIX_H
#define RANDOMSUBMATRIX_H

#include "SubMatrix.h"
#include "Var.h"
#include "RandomTypes.h"

// a RandomSubMatrix is just a SubMatrix (substitution matrix)
// wrapped up into a Dvar<void*>
// so that it can be a node in the graph of the model
// has some specific corrupt/update behaviors (just for computational efficiency reasons)


class RandomTransitionMatrix : public virtual AbstractTransitionMatrix, public virtual Dnode	{

	public:

	virtual ~RandomTransitionMatrix()	{
		SetName("transition matrix");	
	}

	void specialUpdate()	{
		SetParameters();
	}

	protected:

	// needs to be implemented in all random sub matrices
	// is called before any update
	// re-setting all the parameters (pointers and values)
	// by which the matrix can access information from its parents
	virtual void SetParameters() = 0;

	void localRestore()	{
		SetParameters();
		Dnode::localRestore();
	}

	void localCorrupt(bool bk)	{
		CorruptMatrix();
		Dnode::localCorrupt(bk);
	}
};

class RandomSubMatrix : public virtual RandomTransitionMatrix, public virtual SubMatrix {

	public:

	RandomSubMatrix(int Nstate, bool innormalise = false) : SubMatrix(Nstate, innormalise) {}

	virtual ~RandomSubMatrix()	{}

	int GetNstate() {return SubMatrix::GetNstate();}
};

// class RandomTransitionMatrix : public virtual RandomAbstractTransitionMatrix, public virtual

#endif // RANDOMSUBMATRIX_H

