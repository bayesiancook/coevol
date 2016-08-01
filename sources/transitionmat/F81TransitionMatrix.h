
#include "SubMatrix.h"

class F81TransitionMatrix : public virtual TransitionMatrix	{

	public:

	F81TransitionMatrix(int inNstate, const double* stat, double inlength)	{

		flag = false;
		Nstate = inNstate;
		length = inlength;
		Create();
		for (int i=0; i<Nstate; i++)	{
			stationary[i] = stat[i];
		}
	}

	F81TransitionMatrix(const F81TransitionMatrix& from)	{

		flag = false;
		Nstate = from.Nstate;
		length = from.length;
		Create();
		for (int i=0; i<Nstate; i++)	{
			stationary[i] = from.stationary[i];
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				R[i][j] = from.R[i][j];
			}
		}
		// flag = from.flag;
	}

	~F81TransitionMatrix()	{
		Delete();
	}

	F81TransitionMatrix& operator=(const F81TransitionMatrix& from)	{
		if (Nstate != from.Nstate)	{
			cerr << "error in f81: non matching dim\n";
			exit(1);
		}
		length = from.length;
		for (int i=0; i<Nstate; i++)	{
			stationary[i] = from.stationary[i];
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				R[i][j] = from.R[i][j];
			}
		}
		// flag = from.flag;
	}

	int GetNstate() {return Nstate;}

	double** GetR()	{
		if (! flag)	{
			ComputeArray();
		}
		return R;
	}

	protected:

	void CorruptMatrix()	{
		flag = false;
	}

	void Create()	{
		stationary = new double[Nstate];
		// bkstationary = new double[Nstate];
		R = new double*[Nstate];
		// bkR = new double*[Nstate];
		for (int i=0; i<Nstate; i++)	{
			R[i] = new double[Nstate];
			// bkR[i] = new double[Nstate];
		}
	}

	void Delete()	{
		delete[] stationary;
		// delete[] bkstationary;
		for (int i=0; i<Nstate; i++)	{
			delete[] R[i];
			// delete[] bkR[i];
		}
		delete[] R;
		// delete[] bkR;
	}

	/*
	void BackUpMatrix()	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				bkR[i][j] = R[i][j];
			}
		}
	}

	void RestoreMatrix()	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				R[i][j] = bkR[i][j];
			}
		}
	}
	*/

	void ComputeArray()	{

		double expo = exp(-length);
		for (int i=0; i<Nstate; i++)	{
			double total = 0;
			for (int j=0; j<Nstate; j++)	{
				if (i!=j)	{
					R[i][j] = (1 - expo) * stationary[j];
					total += R[i][j];
				}
			}
			R[i][i] = (1 - expo) * stationary[i] + expo;
		}
		flag = true;
	}

	virtual void 		BackwardPropagate(const double* down, double* up, double length) {
		if (! flag)	{
			ComputeArray();
		}
		for (int i=0; i<Nstate; i++)	{
			double tmp = 0;
			for (int j=0; j<Nstate; j++)	{
				tmp += R[i][j] * down[j];
			}
			up[i] = tmp;
		}

	}

	virtual void 		ForwardPropagate(const double* up, double* down, double length) {
		if (! flag)	{
			ComputeArray();
		}
		for (int i=0; i<Nstate; i++)	{
			double tmp = 0;
			for (int j=0; j<Nstate; j++)	{
				tmp += up[j] * R[j][i];
			}
			down[i] = tmp;
		}
	}

	const double* GetStationary()	{
		return stationary;
	}

	const double* GetRow(int i)	{
		if (! flag)	{
			ComputeArray();
		}
		return R[i];
	}

	double Stationary(int i)	{
		return stationary[i];
	}

	double operator()(int i, int j)	{
		if (! flag)	{
			ComputeArray();
		}
		return R[i][j];
	}

	int Nstate;
	double* stationary;
	double** R;
	double length;
	bool flag;

};

class RandomF81TransitionMatrix : public virtual F81TransitionMatrix, public virtual RandomTransitionMatrix {

	public:

	RandomF81TransitionMatrix(Var<Profile>* inRandomStat, Var<PosReal>* inRandomLength) : F81TransitionMatrix(inRandomStat->GetDim(), inRandomStat->GetArray(), 1.0)	{
		RandomStat = inRandomStat;
		RandomLength = inRandomLength;
		Register(RandomStat);
		Register(RandomLength);
		specialUpdate();
	}

	~RandomF81TransitionMatrix() {}

	protected:

	void specialUpdate()	{
		SetParameters();
	}

	void SetParameters()	{
		if (RandomLength)	{
			length = RandomLength->val();
		}
		else	{
			length = 0;
		}
		for (int i=0; i<Nstate; i++)	{
			stationary[i] = (*RandomStat)[i];
		}
	}

	void localRestore()	{
		SetParameters();
		Dnode::localRestore();
	}

	void localCorrupt(bool bk)	{
		CorruptMatrix();
		Dnode::localCorrupt(bk);
	}

	Var<Profile>* RandomStat;
	Var<PosReal>* RandomLength;

};

