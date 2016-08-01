
#ifndef RENPOS_H
#define RENPOS_H

class RenormalizedPosRealVector : public Dvar<Profile>	{

	public:

	RenormalizedPosRealVector(Var<PosRealVector>* inup)	{
		int dimension = inup->GetDim();
		setval(Profile(dimension));
		bkvalue = Profile(dimension);
		up = inup;
		Register(up);
		specialUpdate();
	}

	double GetEntropy()	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot -= (*this)[i] * log((*this)[i]);
		}
		return tot;
	}

	protected:

	void specialUpdate()	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot += (*up)[i];
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = (*up)[i] / tot;
		}
	}

	Var<PosRealVector>* up;

};

#endif

