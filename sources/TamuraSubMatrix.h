

#ifndef TAMURASUBMATRIX_H
#define TAMURASUBMATRIX_H

#include "RandomSubMatrix.h"
#include "BiologicalSequences.h"

class TamuraSubMatrix : public virtual SubMatrix	{

	public:

				TamuraSubMatrix(double intstv, double ingc, bool innormalise = false) : SubMatrix(Nnuc,innormalise)	{
					tstv = intstv;
					gc = ingc;
				}

				~TamuraSubMatrix() {};

	protected:

	bool isTransition(int s1, int s2)	{
		return ((s1 == 0) && (s2==2)) || ((s1==1) && (s2==3) || ((s1==2) && (s2==0)) || ((s1==3) && (s2==1)));
	}

	void 			ComputeArray(int state)	{
					double total = 0;
					for (int k=0; k<Nnuc; k++)	{
						if (k != state)	{
							double rate = 1.0;
							if ((k==0) || (k==3))	{
								rate *= (1 - gc);
							}
							else {
								rate *= gc;
							}
							if (isTransition(k,state))	{
								rate *= tstv;
							}
							Q[state][k] = rate;
							total += rate;
						}
					}
					Q[state][state] = -total;
				}

	void 			ComputeStationary() {
					// compute stationary probabilities
					double total = 0;
					for (int i=0; i<Nnuc; i++)	{
						if ((i==1) || (i==2))	{
							mStationary[i] = 0.5 * gc;
						}
						else	{
							mStationary[i] = 0.5 * (1 - gc);
						}
						total += mStationary[i];
					}
					if (fabs(total -1) > 1e-6)	{
						cerr << "error in Tamura::compute stat\n";
						exit(1);
					}
				}

	double tstv;
	double gc;
};


class TamuraRandomSubMatrix : public RandomSubMatrix, public TamuraSubMatrix  {

	public:

	TamuraRandomSubMatrix(Var<PosReal>* inTsTv, Var<UnitReal>* inGC, bool innormalise = false) : SubMatrix(Nnuc, innormalise), RandomSubMatrix(Nnuc, innormalise) , TamuraSubMatrix(0.5,inGC->val(), innormalise)	{
		TsTv = inTsTv;
		GC = inGC;
		if (TsTv)	{
			Register(TsTv);
		}
		Register(GC);
		SetParameters();
		specialUpdate();
	}

	protected:

	void SetParameters()	{
		if (TsTv)	{
			tstv = TsTv->val();
		}
		else	{
			tstv = 0.5;
		}
		gc = GC->val();
	}

	private:
	Var<PosReal>* TsTv;
	Var<UnitReal>* GC;
};


#endif

