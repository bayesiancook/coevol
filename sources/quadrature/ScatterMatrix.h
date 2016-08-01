#ifndef SCATTERMATRIX_H
#define SCATTERMATRIX_H

#include "CovMatrix.h"
#include "BaseType.h"
#include "ValArray.h"

class ScatterMatrix : public Dvar<CovMatrix>	{

	public:

	ScatterMatrix(VarArray<RealVector>* iniidarray)	{
		iidarray = iniidarray;
		setval(CovMatrix(iidarray->GetVal(0)->GetDim()));
		bkvalue = CovMatrix(iidarray->GetVal(0)->GetDim());
		for (int p=0; p<iidarray->GetSize(); p++)	{
			Register(iidarray->GetVal(p));
		}
		specialUpdate();
	}
		
	void specialUpdate()	{
		ScatterizeOnZero(iidarray);
	}

	private:

	VarArray<RealVector>* iidarray;
};

#endif

