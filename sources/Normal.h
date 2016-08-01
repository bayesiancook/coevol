#ifndef NORMAL_H
#define NORMAL_H


#include "RandomTypes.h"
#include "Var.h"
#include "BiologicalSequences.h"
#include "Move.h"


// iid normal with same mean and variance

class IIDUniform: public virtual Rvar<RealVector>	{

	public:

	IIDUniform(Var<Real>* inroot, int dim, double inmax = 100)	{
		root = inroot;
		Register(root);
		max = inmax;
		setval(RealVector(dim));
		bkvalue = RealVector(dim);
		ClampVector = new bool[dim];
		for (int i=0; i<dim; i++)	{
			ClampVector[i] = false;
		}
		Sample();
	}

	~IIDUniform()	{
		delete[] ClampVector;
	}

	void ClampAt(double inval, int index){
		val()[index] = inval;
		ClampVector[index] = true;
	}

	void ClampAtZero()	{
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = 0;
			bkvalue[i] = 0;
		}
		Clamp();
	}

	double logProb()	{
		return 0;
	}

	virtual double	Move(double tuning, int n)	{
		if (! isClamped())	{
			// Metropolis Hastings here
			Corrupt(true);
			if ((n<=0) || (n > dim))	{
				n = dim;
			}
			int* indices = new int[n];
			Random::DrawFromUrn(indices,n,dim);
			for (int i=0; i<n; i++)	{
				if (! ClampVector[indices[i]])	{
					double& d = vec[indices[i]];
					d += tuning * (Random::Uniform() - 0.5);
					while (fabs(d) > max)	{
						if (d > max)	{
							d = 2*max - d;
						}
						if (d < -max)	{
							d = -2*max - d;
						}
					}
				}
			}
			delete[] indices;
			double deltaLogProb = Update();
			double logRatio = deltaLogProb;
			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
			}
			return (double) accepted;
		}
		return 1;
	}

	/*
	double PiecewiseTranslationMove(double tuning, int index, int n)	{
		if (! isClamped())	{
			Corrupt(true);
			double logHastings = ProposePiecewiseTranslationMove(tuning,index,n);
			double deltaLogProb = Update();
			double logRatio = deltaLogProb + logHastings;
			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
				return 0;
			}
		}
		return 1;
	}

	double ProposePiecewiseTranslationMove(double tuning, int index, int n)	{
		double u = tuning * (Random::Uniform() - 0.5);
		return PiecewiseTranslation(u,index,n);
	}

	double PiecewiseTranslation(double u, int index, int n)	{
		for ( int i=0; i <n; i++){
			(*this)[index + i] += u;
		}
		return 0;
	}
	*/


	double GetMean()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += (*this)[i];
		}
		return total / GetDim();
	}

	double GetVar()	{
		double mean = 0;
		for (int i=0; i<GetDim(); i++)	{
			mean += (*this)[i];
		}
		double var = 0;
		for (int i=0; i<GetDim(); i++)	{
			var += (*this)[i] * (*this)[i];
		}
		mean /= GetDim();
		var /= GetDim();
		var -= mean * mean;
		return var;
	}

	protected:

	void drawSample()	{
		for (int i=0; i<GetDim(); i++)	{
			if (! ClampVector[i])	{
				(*this)[i] = Random::sNormal() * 2 * max - max;
			}
		}
	}
	private:

	bool* ClampVector;
	double max;
	Var<Real>* root;

};

class IIDNormal : public virtual Rvar<RealVector>	{

	public:

	IIDNormal(int dim) {
		setval(RealVector(dim));
		bkvalue = RealVector(dim);
		ClampVector = new bool[dim];
		for (int i=0; i<dim; i++)	{
			ClampVector[i] = false;
		}
	}

	IIDNormal(int dim, Var<Real>* inmean, Var<PosReal>*  invariance)	{
		setval(RealVector(dim));
		bkvalue = RealVector(dim);
		mean = inmean;
		meanvector = 0;
		variance = invariance;
		Register(mean);
		Register(variance);
		ClampVector = new bool[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			ClampVector[i] = false;
		}
		Sample();
	}

	IIDNormal(Var<RealVector>* inmeanvector, Var<PosReal>*  invariance)	{
		setval(RealVector(inmeanvector->GetDim()));
		bkvalue = RealVector(inmeanvector->GetDim());
		meanvector = inmeanvector;
		mean = 0;
		variance = invariance;
		Register(meanvector);
		Register(variance);
		ClampVector = new bool[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			ClampVector[i] = false;
		}
		Sample();
	}

	~IIDNormal()	{
		delete[] ClampVector;
	}

	void ClampAtZero()	{
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = 0;
			bkvalue[i] = 0;
		}
		Clamp();
	}

	void ClampAt(double inval, int index){
		val()[index] = inval;
		ClampVector[index] = true;
	}

	/*
	bool isClamped(int index)	{
		return ClampVector[index];
	}
	*/

	double logProb()	{
		double total = 0;
		if (mean)	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - mean->val();
				total += tmp * tmp;
			}
		}
		else	{
			for (int i=0; i<GetDim(); i++)	{
				double tmp = (*this)[i] - (*meanvector)[i];
				total += tmp * tmp;
			}
		}
		return -0.5 * GetDim() * log(2 * Pi * variance->val()) -0.5 * total / variance->val();
	}

	double ProposeMove(double tuning)	{
		int choose = (int) (3 * Random::Uniform());
		if (choose==0)	{
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					val()[i] += tuning * (Random::Uniform() - 0.5);
				}
			}
		}
		else if (choose==1)	{
			double u = tuning * (Random::Uniform() - 0.5);
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					val()[i] += u;
				}
			}
		}
		else	{
			int i = (int) (GetDim() * Random::Uniform());
			if(!ClampVector[i]){
				val()[i] += tuning * (Random::Uniform() - 0.5);
			}
		}
		return 0;
	}

	double PiecewiseTranslationMove(double tuning, int index, int n)	{
		if (! isClamped())	{
			Corrupt(true);
			double logHastings = ProposePiecewiseTranslationMove(tuning,index,n);
			double deltaLogProb = Update();
			double logRatio = deltaLogProb + logHastings;
			bool accepted = (log(Random::Uniform()) < logRatio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
				return 0;
			}
		}
		return 1;
	}

	double ProposePiecewiseTranslationMove(double tuning, int index, int n)	{
		double u = tuning * (Random::Uniform() - 0.5);
		return PiecewiseTranslation(u,index,n);
	}

	double PiecewiseTranslation(double u, int index, int n)	{
		for ( int i=0; i <n; i++){
			if (! ClampVector[i])	{
				(*this)[index + i] += u;
			}
		}
		return 0;
	}


	double GetMean()	{
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			total += (*this)[i];
		}
		return total / GetDim();
	}

	double GetVar()	{
		double mean = 0;
		for (int i=0; i<GetDim(); i++)	{
			mean += (*this)[i];
		}
		double var = 0;
		for (int i=0; i<GetDim(); i++)	{
			var += (*this)[i] * (*this)[i];
		}
		mean /= GetDim();
		var /= GetDim();
		var -= mean * mean;
		return var;
	}

	protected:

	void drawSample()	{
		if (mean)	{
			for (int i=0; i<GetDim(); i++)	{
				if (! ClampVector[i])	{
					(*this)[i] = Random::sNormal() * sqrt(variance->val()) + mean->val();
				}
			}
		}
		else	{
			for (int i=0; i<GetDim(); i++)	{
				if (! ClampVector[i])	{
					(*this)[i] = Random::sNormal() * sqrt(variance->val()) + (*meanvector)[i];
				}
			}
		}
	}

	bool* ClampVector;
	Var<Real>* mean;
	Var<RealVector>* meanvector;
	Var<PosReal>* variance;

};


class IIDAddition : public Dvar<RealVector> {

	public :

	IIDAddition(Var<RealVector>* ina, Var<RealVector>* inb)	{
		a = ina;
		b = inb;
		Register(a);
		Register(b);
		specialUpdate();
	}

	void specialUpdate()	{
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = (*a)[i] + (*b)[i];
		}
	}

	private:

	Var<RealVector>* a;
	Var<RealVector>* b;

};

class Addition : public Dvar<Real>	{

	public:

	Addition(Var<Real>* ina, Var<Real>* inb)	{
		a = ina;
		b = inb;
		Register(a);
		Register(b);
		specialUpdate();
	}

	void specialUpdate()	{
		setval(a->val() + b->val());
	}

	private:
	Var<Real>* a;
	Var<Real>* b;

};

#endif

