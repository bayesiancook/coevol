
#ifndef NONREVCODON_H
#define NONREVCODON_H

#include "CodonSubMatrix.h"
#include "NonRevSubMatrix.h"

class NonRevCodonSubMatrix : public virtual NonRevSubMatrix, public virtual CodonSubMatrix	{

	public:

	NonRevCodonSubMatrix(CodonStateSpace* instatespace, const double* inrate, double inlength, bool innormalise = false, int indiscn = 10)	:
		SubMatrix(instatespace->GetNstate(),innormalise),
		NonRevSubMatrix(instatespace->GetNstate(),innormalise,indiscn),
		CodonSubMatrix(instatespace,innormalise)	{

			rate = inrate;

	}

	~NonRevCodonSubMatrix() {}


	protected:

	void 			ComputeStationary()	{}

	void SetRates(const double* inrate) {rate = inrate;}

	void ComputeArray(int cod1)	{

		double total = 0;
		for (int cod2=0; cod2<GetNstate(); cod2++)	{
			if (cod1!=cod2)	{
				int pos = GetDifferingPosition(cod1,cod2);
				if ((pos != -1) && (pos != 3))	{

					int i = GetCodonPosition(pos,cod1);
					int j = GetCodonPosition(pos,cod2);

					int index = 0;
					if ( ((i==0) && (j==3)) || ((i==3) && (j==0)) )	{
						// tmp = at2ta;
						index = 0;
					}
					if ( ((i==1) && (j==2)) || ((i==2) && (j==1)) )	{
						// tmp = cg2gc;
						index = 1;
					}
					if ( ((i==0) && (j==1)) || ((i==3) && (j==2)) )	{
						// tmp = at2cg;
						index = 2;
					}
					if ( ((i==1) && (j==0)) || ((i==2) && (j==3)) )	{
						// tmp = cg2at;
						index = 3;
					}
					if ( ((i==0) && (j==2)) || ((i==3) && (j==1)) ) {
						// tmp = at2gc;
						index = 4;
					}
					if ( ((i==1) && (j==3)) || ((i==2) && (j==0)) )	{
						// tmp = cg2ta;
						index = 5;
					}
					Q[cod1][cod2] = rate[index];
					if (! Synonymous(cod1,cod2))	{
						Q[cod1][cod2] *= rate[index + 6];
					}
					total += Q[cod1][cod2];
				}
				else	{
					Q[cod1][cod2] = 0;
				}
			}
		}
		Q[cod1][cod1] = -total;
	}

	const double* rate;

};

class RandomNonRevCodonSubMatrix : public virtual RandomNonRevSubMatrix, public virtual NonRevCodonSubMatrix	{

	public:

	RandomNonRevCodonSubMatrix(CodonStateSpace* instatespace, Var<PosReal>* inLength, Var<PosRealVector>* inRate, bool innormalise = false, int indiscn = 10)	:

		SubMatrix(instatespace->GetNstate(), innormalise),
		NonRevSubMatrix(instatespace->GetNstate(), innormalise, indiscn),
		RandomSubMatrix(instatespace->GetNstate(), innormalise) ,
		RandomNonRevSubMatrix(instatespace->GetNstate(), innormalise, indiscn),
		CodonSubMatrix(instatespace,innormalise),
		NonRevCodonSubMatrix(instatespace,inRate->GetArray(),innormalise,indiscn)	{

			Rate = inRate;
			Length = inLength;
			Register(Rate);
			Register(Length);
			specialUpdate();
	}


	protected:

	void SetParameters()	{
		SetRates(Rate->GetArray());
		SetLength(Length->val());
	}

	Var<PosRealVector>* Rate;
	Var<PosReal>* Length;
};

class RootCodonSubMatrix : public virtual CodonSubMatrix	{


	public:

	RootCodonSubMatrix(CodonStateSpace* instatespace, double ingc, bool innormalise = false) :
		SubMatrix(instatespace->GetNstate(),innormalise),
		CodonSubMatrix(instatespace,innormalise)	{
			gc = ingc;
			stat = new double[Nnuc];
	}

	~RootCodonSubMatrix()	{
		delete[] stat;
	}

	protected:

	void ComputeStationary()	{
		stat[0] = stat[3] = (1 - gc) / 2;
		stat[1] = stat[2] = gc / 2;
		double total = 0;
		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i] = stat[GetCodonPosition(0,i)] * stat[GetCodonPosition(1,i)] * stat[GetCodonPosition(2,i)];
			total += mStationary[i];
		}

		for (int i=0; i<GetNstate(); i++)	{
			mStationary[i]  /= total;
		}
	}

	void ComputeArray(int cod1)	{
	}

	double gc;
	double* stat;
};

class RandomRootCodonSubMatrix : public virtual RandomSubMatrix, public virtual RootCodonSubMatrix	{

	public:

	RandomRootCodonSubMatrix(CodonStateSpace* instatespace, Var<UnitReal>* inGC, bool innormalise = false)	:
		SubMatrix(instatespace->GetNstate(), innormalise),
		RandomSubMatrix(instatespace->GetNstate(), innormalise),
		CodonSubMatrix(instatespace,innormalise),
		RootCodonSubMatrix(instatespace,inGC->val(),innormalise)	{

			GC = inGC;
			Register(GC);
	}

	protected:

	void SetParameters()	{
		gc = GC->val();
	}

	Var<UnitReal>* GC;

};

#endif

