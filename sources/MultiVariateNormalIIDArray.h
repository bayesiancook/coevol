
#ifndef "MVNORMIID_H"
#define "MVNORMIID_H"

#include "IID.h"
#include "ConjugateInverseWishart.h"

class ConjugateMultiVariateNormalIIDArray : public IIDArray<RealVector>	{

	public:
	// ConjugateMultiVariareNormalIIDArray(int insize, Var<CovMatrix>* inmatrix) : IIDArray<RealVector>(insize)	{
	ConjugateMultiVariareNormalIIDArray(int insize, MultiNormalSemiConjugate* insigma) : IIDArray<RealVector>(insize)	{
		sigma = insigma;
		Create();
	}

	protected:

	Rvar<RealVector>* CreateVal(int site)	{
		return new ConjugateMultiVariateNormal(sigma);
	}

	MultiNormalSemiConjugate* sigma;
};

#endif

