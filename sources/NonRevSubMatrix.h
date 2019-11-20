
#ifndef NONREVSUBMATRIX_H
#define NONREVSUBMATRIX_H

#include "SubMatrix.h"

class NonRevSubMatrix : public virtual SubMatrix	{

	public:

	NonRevSubMatrix(int inNstate, bool innormalise = false, int indiscn = 10) : SubMatrix(inNstate,innormalise), discn(indiscn), expoflag(false) {
		expu = new double*[GetNstate()];
		expu2 = new double*[GetNstate()];
		for (int i=0; i<GetNstate(); i++)	{
			expu[i] = new double[GetNstate()];
			expu2[i] = new double[GetNstate()];
		}
		length = 0;
	}

	virtual ~NonRevSubMatrix()	{
		for (int i=0; i<GetNstate(); i++)	{
			delete[] expu[i];
			delete[] expu2[i];
		}
		delete[] expu;
		delete[] expu2;
	}

	void SetLength(double inlength)	{
		// expoflag = false;
		length = inlength;
	}

	virtual void ComputeExponential(double length);

	public:

	virtual void BackwardPropagate(const double* down, double* up, double length);
	virtual void ForwardPropagate(const double* up, double* down, double length);

	protected:

	void CorruptMatrix()	{
		SubMatrix::CorruptMatrix();
		expoflag = false;
	}

	double ** expu;
	double ** expu2;
	double length;
	int discn;
	bool expoflag;

};

class RandomNonRevSubMatrix : public virtual NonRevSubMatrix, public virtual RandomSubMatrix	{

	public:

	RandomNonRevSubMatrix(int Nstate, bool innormalise = false, int indisc = 10) :
		SubMatrix(Nstate, innormalise),
		NonRevSubMatrix(Nstate, innormalise, indisc),
		RandomSubMatrix(Nstate, innormalise) {}

	virtual ~RandomNonRevSubMatrix()	{}

};

void NonRevSubMatrix::ComputeExponential(double inlength)	{

	if (! expoflag)	{
	expoflag = true;
	if (length != inlength)	{
		cerr << "error : non matching length in non rev matrix\n";
		cerr << length << '\t' << inlength << '\n';
		exit(1);
	}

	for (int i=0; i<GetNstate(); i++)	{
		if (! flagarray[i])	{
			UpdateRow(i);
		}
	}

	double t = length;
	for (int i=0; i<discn; i++)	{
		t /= 4.0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			expu[i][j] = t * Q[i][j];
		}
	}
	for (int i=0; i<GetNstate(); i++)	{
		expu[i][i] += 1.0;
	}

	for (int n=0; n<discn; n++)	{
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu[i][k] * expu[k][j];
				}
				expu2[i][j] = tmp;
			}
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu2[i][k] * expu2[k][j];
				}
				expu[i][j] = tmp;
			}
		}
	}
	}
}

void NonRevSubMatrix::BackwardPropagate(const double* up, double* down, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		down[i] = 0;
	}
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			down[i] += expu[i][j] * up[j];
		}
	}

	for (int i=0; i<GetNstate(); i++)	{
		if (std::isnan(down[i]))	{
			cerr << "error in back prop\n";
			cerr << "non rev\n";
			for (int j=0; j<GetNstate(); j++)	{
				cerr << up[j] << '\t' << down[j] << '\n';
			}
			for (int j=0; j<GetNstate(); j++)	{
				for (int k=0; k<GetNstate(); k++)	{
					cerr << expu[j][k] << '\t';
				}
				cerr << '\n';
			}
			cerr << '\n';
			for (int j=0; j<GetNstate(); j++)	{
				for (int k=0; k<GetNstate(); k++)	{
					cerr << Q[j][k] << '\t';
				}
				cerr << '\n';
			}
			cerr << "length : " << length << '\n';
			exit(1);
		}
	}
	double maxup = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (up[k] <0)	{
			cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
			// down[k] = 0;
		}
		if (maxup < up[k])	{
			maxup = up[k];
		}
	}
	double max = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (down[k] <0)	{
			down[k] = 0;
		}
		if (max < down[k])	{
			max = down[k];
		}
	}
	if (maxup == 0)	{
		cerr << "error in backward propagate: null up array\n";
		exit(1);
	}
	if (max == 0)	{
		cerr << "error in backward propagate: null array\n";
		for (int k=0; k<GetNstate(); k++)	{
			cerr << up[k] << '\t' << down[k] << '\n';
		}
		cerr << '\n';
		exit(1);
	}
	down[GetNstate()] = up[GetNstate()];
}

/*
inline void NonRevSubMatrix::FiniteTime(int i0, double* up, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		up[i] = expu[i0][i];
	}
}
*/

void NonRevSubMatrix::ForwardPropagate(const double* down, double* up, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		up[i] = 0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			up[i] += down[j] * expu[j][i];
		}
	}
}

#endif

/*
inline void SubMatrix::ComputeExponential(double length)	{

	for (int i=0; i<GetNstate(); i++)	{
		if (! flagarray[i])	{
			UpdateRow(i);
		}
	}

	double t = length;
	for (int i=0; i<discn; i++)	{
		t /= 4.0;
	}

	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			expu[i][j] = t * Q[i][j];
		}
	}
	for (int i=0; i<GetNstate(); i++)	{
		expu[i][i] += 1.0;
	}

	for (int n=0; n<discn; n++)	{
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu[i][k] * expu[k][j];
				}
				expu2[i][j] = tmp;
			}
		}
		for (int i=0; i<GetNstate(); i++)	{
			for (int j=0; j<GetNstate(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetNstate(); k++)	{
					tmp += expu2[i][k] * expu2[k][j];
				}
				expu[i][j] = tmp;
			}
		}
	}
}

inline void SubMatrix::BackwardPropagate(const double* up, double* down, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		down[i] = 0;
	}
	for (int i=0; i<GetNstate(); i++)	{
		for (int j=0; j<GetNstate(); j++)	{
			down[i] += expu[i][j] * up[j];
		}
	}

	for (int i=0; i<GetNstate(); i++)	{
		if (std::isnan(down[i]))	{
			cerr << "error in back prop\n";
			for (int j=0; j<GetNstate(); j++)	{
				cerr << up[j] << '\t' << down[j] << '\n';
			}
			exit(1);
		}
	}
	double maxup = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (up[k] <0)	{
			cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
			// down[k] = 0;
		}
		if (maxup < up[k])	{
			maxup = up[k];
		}
	}
	double max = 0;
	for (int k=0; k<GetNstate(); k++)	{
		if (down[k] <0)	{
			down[k] = 0;
		}
		if (max < down[k])	{
			max = down[k];
		}
	}
	if (maxup == 0)	{
		cerr << "error in backward propagate: null up array\n";
		exit(1);
	}
	if (max == 0)	{
		cerr << "error in backward propagate: null array\n";
		for (int k=0; k<GetNstate(); k++)	{
			cerr << up[k] << '\t' << down[k] << '\n';
		}
		cerr << '\n';
		exit(1);
	}
	down[GetNstate()] = up[GetNstate()];
}

inline void SubMatrix::FiniteTime(int i0, double* up, double length)	{

	ComputeExponential(length);
	for (int i=0; i<GetNstate(); i++)	{
		up[i] = expu[i0][i];
	}
}
*/
