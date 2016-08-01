#ifndef WISHARRAY_H
#define WISHARRAY_H

class InverseWishartArray : public IIDArray<CovMatrix>	{
// class InverseWishartArray : public ValPtrArray<Rvar<CovMatrix> >	{

	public:

	InverseWishartArray(int insize, int indf, SigmaZero* insigmazero, SigmaZero* insigmazero0) : IIDArray<CovMatrix>(insize)	{
		df = indf;
		sigmazero = insigmazero;
		sigmazero0 = insigmazero0;
		Create();
	}

	protected:

	Rvar<CovMatrix>* CreateVal(int mat)	{
		if (!mat)	{
			return new InverseWishartMatrix(sigmazero,sigmazero->GetDim() + df);
		}
		return new InverseWishartMatrix(sigmazero0,sigmazero0->GetDim() + df);
	}

	private:

	int df;
	SigmaZero* sigmazero;
	SigmaZero* sigmazero0;

};

class DiagonalWishartArray : public IIDArray<CovMatrix>	{

	public:

	DiagonalWishartArray(int insize, int indf, SigmaZero* insigmazero, SigmaZero* insigmazero0) : IIDArray<CovMatrix> (insize)	{
		df = indf;
		sigmazero = insigmazero;
		sigmazero0 = insigmazero0;
		Create();
	}

	protected:

	Rvar<CovMatrix>* CreateVal(int mat)	{
		if (!mat)	{
			return new DiagonalCovMatrix(sigmazero,sigmazero->GetDim() + df);
		}
		return new DiagonalCovMatrix(sigmazero0,sigmazero0->GetDim() + df);
	}

	private:

	int df;
	SigmaZero* sigmazero;
	SigmaZero* sigmazero0;

};

class ConjugateInverseWishartArray : public IIDArray<CovMatrix>	{

	public:

	ConjugateInverseWishartArray(int insize, int indf, SigmaZero* insigmazero, SigmaZero* insigmazero0) : IIDArray<CovMatrix>(insize)	{
		df = indf;
		sigmazero = insigmazero;
		sigmazero0 = insigmazero0;
		Create();
	}

	protected:

	Rvar<CovMatrix>* CreateVal(int mat)	{
		if (!mat)	{
			return new ConjugateInverseWishart(sigmazero,sigmazero->GetDim() + df);
		}
		return new ConjugateInverseWishart(sigmazero0,sigmazero0->GetDim() + df);
	}

	private:

	int df;
	SigmaZero* sigmazero;
	SigmaZero* sigmazero0;

};

#endif

