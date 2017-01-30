#ifndef AMINOACIDOMEGASUBMATRIX_H
#define AMINOACIDOMEGASUBMATRIX_H

#include "RandomSubMatrix.h"
#include "GTRSubMatrix.h"
#include "SimilarityMatrix.h"
#include "CodonSubMatrix.h"
#include "SplitAAMatrix.h"

class AAOmegaSubMatrix : public GTRSubMatrix {
	public:

	AAOmegaSubMatrix (SimilarityMatrix* inmatrix,const double* inrelrate, const double* instat, double inomega, bool innormalise):
				SubMatrix (Naa, innormalise),
				GTRSubMatrix(Naa,inrelrate,instat,innormalise){
					mSimilarityMatrix = inmatrix;
					omega = inomega;
				};

	double GetOmega(){return omega;}

	protected:

	void setOmega(double inomega) {
		omega = inomega;
	}
	void	ComputeArray(int state);

	double omega;
	SimilarityMatrix* mSimilarityMatrix;

};

class RandomAAOmegaSubMatrix : public RandomSubMatrix, public AAOmegaSubMatrix {

	public:
 
	RandomAAOmegaSubMatrix(SimilarityMatrix* inmatrix, Var<Profile>* inrelrate, Var<Profile>* instat ,Var<PosReal>* inomega, bool innormalise = false):
				SubMatrix (Naa, innormalise),
				RandomSubMatrix (Naa,innormalise),
				AAOmegaSubMatrix(inmatrix,inrelrate->val().GetArray(), instat->val().GetArray(), inomega->val(), innormalise){

					stat = instat;
					relrate = inrelrate;
					randomOmega = inomega;
					rescaledrelrate = new double[relrate->GetDim()];
					Register(relrate);
					Register(stat);
					Register(randomOmega);
					specialUpdate();
				};

	protected :

	void SetParameters();

	private:
	Var<Profile>* stat;
	Var<PosReal>* randomOmega;
	Var<Profile>* relrate;
	double* rescaledrelrate;
};

class AADoubleOmegaSubMatrix : public MGCodonSubMatrix	{

	public:

	AADoubleOmegaSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, SimilarityMatrix* inmatrix, double inomegadNdS, double inomegaKrKc, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise){
				SetOmega(inomegadNdS, inomegaKrKc);
				mSimilarityMatrix = inmatrix;
				statespace = instatespace;
			};

	double		GetOmegadNdS() {return omegadNdS;}
	double		GetOmegaKrKc() {return omegaKrKc;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void	ComputeArray(int state);
	void	SetOmega(double inomegadNdS, double inomegaKrKc);

	// data members
	double omegadNdS;
	double omegaKrKc;
	SimilarityMatrix* mSimilarityMatrix;
	CodonStateSpace* statespace;
};

class RandomAADoubleOmegaSubMatrix : public AADoubleOmegaSubMatrix, public RandomCodonSubMatrix {

	public:

	RandomAADoubleOmegaSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inNucMatrix, Var<PosReal>* inRandomOmegadNdS, Var<PosReal>* inRandomOmegaKrKc, SimilarityMatrix* inMatrix, bool innormalise = false):
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			AADoubleOmegaSubMatrix(instatespace,inNucMatrix,inMatrix, inRandomOmegadNdS->val(), inRandomOmegaKrKc->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise){
				matrix = inNucMatrix;
				RandomOmegadNdS = inRandomOmegadNdS;
				RandomOmegaKrKc = inRandomOmegaKrKc;
				Register(matrix);
				Register(RandomOmegadNdS);
				Register(RandomOmegaKrKc);
				specialUpdate();
	};

	protected:

	void SetParameters();

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomOmegadNdS;
	Var<PosReal>* RandomOmegaKrKc;
};

class AAregSubMatrix : public virtual SubMatrix {

	public :

	AAregSubMatrix(const double* instat, double inomega, double* slope, double* intercept, bool innormalise = false ) : SubMatrix(Naa, innormalise){
				omega = inomega;
				mSlope = slope;
				mIntercept = intercept;
				CopyStationary(instat);
			};

	double GetOmega(){return omega;}

	double GetSlope(int i, int j) {return mSlope[index(i,j,GetNstate())];}
	double GetIntercept(int i, int j){return mIntercept[index(i,j,GetNstate())];}

	static int index(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}
 
	protected:

	void	setOmega(double inomega);
	void	ComputeArray(int state);
	void	ComputeStationary(){};
	void	CopyStationary (const double* instat);
	void 	CopySlopeIntercept(double* inslope, double* inintercept, int dimension);

	private:
	double omega;
	double* mIntercept;
	double* mSlope;
};

class RandomAAregSubMatrix : public RandomSubMatrix,  public AAregSubMatrix {

	public :

	RandomAAregSubMatrix (Var<Profile>* instat, Var<PosReal>* inomega, Var<RealVector>* inslope, Var<RealVector>* inintercept, bool innormalise = false ):
			SubMatrix (Naa, innormalise),
			RandomSubMatrix (Naa, innormalise),
			AAregSubMatrix(instat->val().GetArray(), inomega->val(), inslope->val().GetArray(), inintercept->val().GetArray(), innormalise){
				stat = instat;
				randomOmega = inomega;
				slope = inslope;
				intercept = inintercept;
				Register(stat);
				Register(slope);
				Register(intercept);
				Register(randomOmega);
			};

	protected:
	void SetParameters();

	private:
	Var<Profile>* stat;
	Var<PosReal>* randomOmega;
	Var<RealVector>* slope;
	Var<RealVector>* intercept;
};


// split pairs of amino-acids depending on
// whether they include C to T transitions or not
// and then, whether they are radical or not
class SplitAAOmegaSubMatrix : public GTRSubMatrix {
	public:

	SplitAAOmegaSubMatrix (SimilarityMatrix* inmatrix,const double* inrelrate, const double* instat, double intstv, double inomegats, double inomegatv, SplitAAMatrix* insplitmat, int insplittype, bool innormalise):
				SubMatrix (Naa, innormalise),
				GTRSubMatrix(Naa,inrelrate,instat,innormalise){
					mSimilarityMatrix = inmatrix;
					tstv = intstv;
					omegats = inomegats;
					omegatv = inomegatv;
					mSplitMat = insplitmat;
					splittype = insplittype;
				};

	double GetTsTv() {return tstv;}
	double GetOmegaTs(){return omegats;}
	double GetOmegaTv(){return omegatv;}

	protected:

	void setTsTv(double intstv) {
		tstv = intstv;
	}
	void setOmegaTs(double inomegats) {
		omegats = inomegats;
	}
	void setOmegaTv(double inomegatv) {
		omegatv = inomegatv;
	}
	void	ComputeArray(int state);

	double tstv;
	double omegats;
	double omegatv;
	SimilarityMatrix* mSimilarityMatrix;
	SplitAAMatrix* mSplitMat;
	int splittype;

};

class RandomSplitAAOmegaSubMatrix : public RandomSubMatrix, public SplitAAOmegaSubMatrix {

	public:
 
	RandomSplitAAOmegaSubMatrix(SimilarityMatrix* inmatrix, Var<Profile>* inrelrate, Var<Profile>* instat , Var<PosReal>* intstv, Var<PosReal>* inomegats, Var<PosReal>* inomegatv, SplitAAMatrix* insplitaaMatrix, int insplittype, bool innormalise = false):
				SubMatrix (Naa, innormalise),
				RandomSubMatrix (Naa,innormalise),
				SplitAAOmegaSubMatrix(inmatrix,inrelrate->val().GetArray(), instat->val().GetArray(), intstv->val(), 1 , inomegatv->val(), insplitaaMatrix, insplittype, innormalise){

					stat = instat;
					relrate = inrelrate;
					randomTsTv = intstv;
					randomOmegaTs = inomegats;
					randomOmegaTv = inomegatv;
					rescaledrelrate = new double[relrate->GetDim()];
					Register(relrate);
					Register(stat);
					Register(randomTsTv);
					Register(randomOmegaTs);
					Register(randomOmegaTv);
					specialUpdate();
				};

	protected :

	void SetParameters();

	private:
	Var<Profile>* stat;
	Var<PosReal>* randomTsTv;
	Var<PosReal>* randomOmegaTs;
	Var<PosReal>* randomOmegaTv;
	Var<Profile>* relrate;
	double* rescaledrelrate;
};

#endif






