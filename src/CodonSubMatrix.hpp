
#ifndef CODONSUBMATRIX_H
#define CODONSUBMATRIX_H

#include "RandomSubMatrix.h"
#include "CodonStateSpace.h"

const double omegamin = 1e-10;

// a general class representing all codon matrices
// this is still an abstract class
class CodonSubMatrix : public virtual SubMatrix	{

	public:

	CodonSubMatrix(CodonStateSpace* instatespace, bool innormalise)	: SubMatrix(instatespace->GetNstate(),innormalise), statespace(instatespace) {}

	CodonStateSpace*	GetCodonStateSpace() {return statespace;}
	// this is just to avoid repeating "GetCodonStateSpace()->" all the time...
	// int			GetNstate() {return statespace->GetNstate();}
	bool			Synonymous(int codon1, int codon2) {return statespace->Synonymous(codon1,codon2);}
	int			GetCodonPosition(int pos, int codon) {return statespace->GetCodonPosition(pos,codon);}
	int			GetDifferingPosition(int codon1, int codon2) {return statespace->GetDifferingPosition(codon1,codon2);}

	protected:
	CodonStateSpace* statespace;
};


// most codon matrices rely on a mutation process at the level of nucleotides
// thus we create a NucCodonSubMatrix, which takes a substitution matrix representing the nucleotide 4x4 mutation process
// this is still an abstract class
class NucCodonSubMatrix : public virtual CodonSubMatrix	{

	public:

	NucCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, bool innormalise) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise){

		SetNucMatrix(inNucMatrix);
	}

	SubMatrix*		GetNucMatrix()	{
		return NucMatrix;
	}

	protected:

	void			SetNucMatrix(SubMatrix* inmatrix)	{
		NucMatrix = inmatrix;
		if (NucMatrix->GetNstate() != Nnuc)	{
			cerr << "error in CodonSubMatrix: underyling mutation process should be a 4x4 matrix\n";
			throw;
		}
	}

	SubMatrix*	 NucMatrix;

};

//
// this is still an abstract class
class RandomCodonSubMatrix : public virtual CodonSubMatrix, public virtual RandomSubMatrix {

	public:

	RandomCodonSubMatrix(CodonStateSpace* instatespace, bool innorm = false) :
			SubMatrix(instatespace->GetNstate(), innorm),
			CodonSubMatrix(instatespace, innorm),
			RandomSubMatrix(instatespace->GetNstate(), innorm) {}

};


// The Muse and Gaut codon substitution process
// The simplest codon model based on a pure nucleotide mutation process (with stops excluded)
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGCodonSubMatrix : public NucCodonSubMatrix	{

	public:

	MGCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			NucCodonSubMatrix(instatespace,inNucMatrix,innormalise) {

		synnucarray = new double*[Nnuc];
		nonsynnucarray = new double*[Nnuc];
		for (int i=0; i<Nnuc; i++)	{
			synnucarray[i] = new double[Nnuc];
			nonsynnucarray[i] = new double[Nnuc];
			for (int j=0; j<Nnuc; j++)	{
				synnucarray[i][j] = 0;
				nonsynnucarray[i][j] = 0;
			}
		}
		nucflag = false;
	}

	virtual ~MGCodonSubMatrix()	{
		for (int i=0; i<Nnuc; i++)	{
			delete[] synnucarray[i];
			delete[] nonsynnucarray[i];
		}
	}

	double** GetSynNucArray()	{
		if (! nucflag)	{
			ComputeNucArrays();
			nucflag = true;
		}
		return synnucarray;
	}

	double** GetNonSynNucArray()	{
		if (! nucflag)	{
			ComputeNucArrays();
			nucflag = true;
		}
		return nonsynnucarray;
	}

	double GetNormStat()	{
		if (! nucflag)	{
			ComputeNucArrays();
			nucflag = true;
		}
		return 1 - stopstat;
	}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void 			ComputeStationary();

	virtual void		ComputeNucArrays();

	virtual void CorruptMatrix()	{
		nucflag = false;
		SubMatrix::CorruptMatrix();
	}

	double**		synnucarray;
	double**		nonsynnucarray;
	double			stopstat;
	bool			nucflag;

};

// MGCodonSubMatrix is a SubMatrix
// to use it in a probabilistic model, we need to make a RandomSubMatrix object
//
class RandomMGCodonSubMatrix : public MGCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGCodonSubMatrix(instatespace,inmatrix, innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix) {

		Register(matrix);
		specialUpdate();
	}

	// before updating the matrix instant rates and stationary probabilities
	// we just need to check that we are pointing on the right nucleotide mutation process
	void SetParameters()	{
		SetNucMatrix(matrix);
	}

	protected:

	RandomSubMatrix* matrix;
};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmegaCodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmegaCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double inomega, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			omega(inomega) {}


	double			GetOmega() {return omega + omegamin;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmega(double inomega) {omega = inomega;}

	virtual void		ComputeNucArrays();

	void ToStream(ostream& os)	{
		os << "Omega : " << omega << '\n';
		os << "nuc matrix\n";
		GetNucMatrix()->ToStream(os);
		os << '\n';
		SubMatrix::ToStream(os);
	}

	// data members

	double omega;
};

class RandomMGOmegaCodonSubMatrix : public MGOmegaCodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmegaCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inRandomOmega, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGOmegaCodonSubMatrix(instatespace,inmatrix,inRandomOmega->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			RandomOmega(inRandomOmega) {

		Register(matrix);
		Register(RandomOmega);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetOmega(RandomOmega->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomOmega;
};

// this class implements the projection of a 61x61 codon substitution process
// onto a 20x20 amino-acid replacement process
// according to the formula:
//
// R_{ab} = [ \sum \pi_i Q_{ij} ] [ \sum_i \pi_i ]
//
// where a,b runs over amino-acids
// and i (resp. j) over all codons encoding for amino acid a (resp. b)
//
class AminoAcidReducedCodonSubMatrix : public virtual SubMatrix	{

	public:
	AminoAcidReducedCodonSubMatrix(CodonSubMatrix* incodonmatrix, bool innormalise = false) : SubMatrix(Naa,innormalise), codonmatrix(incodonmatrix) {
		aastatespace = new ProteinStateSpace();
	}


	CodonSubMatrix* GetCodonMatrix() {return codonmatrix;}
	CodonSubMatrix* GetCodonSubMatrix() {return codonmatrix;}
	CodonStateSpace* GetCodonStateSpace() {return GetCodonMatrix()->GetCodonStateSpace();}
	ProteinStateSpace* GetProteinStateSpace() {return aastatespace;}

	protected:

	void SetCodonMatrix(CodonSubMatrix* incodonmatrix)	{
		codonmatrix = incodonmatrix;
	}

	void 			ComputeArray(int state);
	void 			ComputeStationary();

	CodonSubMatrix* codonmatrix;
	ProteinStateSpace* aastatespace;

};

// the random class corresponding to AminoAcidReducedCodonSubMatrix
// thus, for any codon model, you can immediately instantiate a projected aminoacid replacement process
// and use it in a probabilistic model using a RandomAminoAcidReducedCodonSubMatrix
//
// it only takes the time to type the name...

class RandomAminoAcidReducedCodonSubMatrix : public AminoAcidReducedCodonSubMatrix, public RandomSubMatrix	{

	public:

	RandomAminoAcidReducedCodonSubMatrix(RandomCodonSubMatrix* inrandomcodonmatrix, bool innormalise = false) :
			SubMatrix(Naa, innormalise),
			AminoAcidReducedCodonSubMatrix(inrandomcodonmatrix,innormalise) ,
			RandomSubMatrix(Naa, innormalise) ,
			randomcodonmatrix(inrandomcodonmatrix)	{

		Register(randomcodonmatrix);
		specialUpdate();

	}

	void SetParameters()	{
		SetCodonMatrix(randomcodonmatrix);
	}

	protected:

	RandomCodonSubMatrix* randomcodonmatrix;
};


#endif


