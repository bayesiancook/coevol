
#ifndef GC3CODON
#define GC3CODON

#include "CodonSubMatrix.h"


// most codon matrices rely on a mutation process at the level of nucleotides
// thus we create a NucCodonSubMatrix, which takes a substitution matrix representing the nucleotide 4x4 mutation process
// this is still an abstract class
class Nuc3CodonSubMatrix : public CodonSubMatrix	{

	public:

	Nuc3CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix1, SubMatrix* inNucMatrix2, SubMatrix* inNucMatrix3, bool innormalise) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise){

		SetNucMatrix1(inNucMatrix1);
		SetNucMatrix2(inNucMatrix2);
		SetNucMatrix3(inNucMatrix3);
	}

	SubMatrix*		GetNucMatrix(int pos)	{
					if (pos == 0)	{
						return NucMatrix1;
					}
					else if (pos == 1)	{
						return NucMatrix2;
					}
					else if (pos == 2)	{
						return NucMatrix3;
					}
					cerr << "error in Nuc3CodonSubMatrix: " << pos << '\n';
					exit(1);
					return 0;
			}

	protected:

	void			SetNucMatrix1(SubMatrix* inmatrix)	{
		NucMatrix1 = inmatrix;
		if (NucMatrix1->GetNstate() != Nnuc)	{
			cerr << "error in CodonSubMatrix: underyling mutation process should be a 4x4 matrix\n";
			throw;
		}
	}

	void			SetNucMatrix2(SubMatrix* inmatrix)	{
		NucMatrix2 = inmatrix;
		if (NucMatrix2->GetNstate() != Nnuc)	{
			cerr << "error in CodonSubMatrix: underyling mutation process should be a 4x4 matrix\n";
			throw;
		}
	}

	void			SetNucMatrix3(SubMatrix* inmatrix)	{
		NucMatrix3 = inmatrix;
		if (NucMatrix3->GetNstate() != Nnuc)	{
			cerr << "error in CodonSubMatrix: underyling mutation process should be a 4x4 matrix\n";
			throw;
		}
	}

	SubMatrix*	 NucMatrix1;
	SubMatrix*	 NucMatrix2;
	SubMatrix*	 NucMatrix3;

};

class MG3CodonSubMatrix : public Nuc3CodonSubMatrix	{

	public:

	MG3CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix1, SubMatrix* inNucMatrix2, SubMatrix* inNucMatrix3, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			Nuc3CodonSubMatrix(instatespace,inNucMatrix1, inNucMatrix2, inNucMatrix3, innormalise) {}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void 			ComputeStationary();

};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MG3OmegaCodonSubMatrix : public MG3CodonSubMatrix	{

	public:

	MG3OmegaCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix1, SubMatrix* inNucMatrix2, SubMatrix* inNucMatrix3, double inomega, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			MG3CodonSubMatrix(instatespace, inNucMatrix1, inNucMatrix2, inNucMatrix3, innormalise) ,
			omega(inomega) {}


	double			GetOmega() {return omega;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmega(double inomega) {omega = inomega;}

	// data members

	double omega;
};

class RandomMG3OmegaCodonSubMatrix : public MG3OmegaCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMG3OmegaCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix1, RandomSubMatrix* inmatrix2, RandomSubMatrix* inmatrix3, Var<PosReal>* inRandomOmega, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MG3OmegaCodonSubMatrix(instatespace,inmatrix1,inmatrix2,inmatrix3,inRandomOmega->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix1(inmatrix1),
			matrix2(inmatrix2),
			matrix3(inmatrix3),
			RandomOmega(inRandomOmega) {

		Register(matrix1);
		Register(matrix2);
		Register(matrix3);
		Register(RandomOmega);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix1(matrix1);
		SetNucMatrix2(matrix2);
		SetNucMatrix3(matrix3);
		SetOmega(RandomOmega->val());
	}

	protected:

	RandomSubMatrix* matrix1;
	RandomSubMatrix* matrix2;
	RandomSubMatrix* matrix3;
	Var<PosReal>* RandomOmega;
};

#endif

