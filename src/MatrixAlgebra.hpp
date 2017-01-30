
#ifndef MATAL_H
#define MATAL_H

#include "CovMatrix.h"

class MatrixAlgebra	{

	public:

	MatrixAlgebra() {}
	~MatrixAlgebra() {}

	protected:

	double** MatrixCreate(int m, int n)	{
		double** tmp = new double*[m];
		for (int i=0; i<m; i++)	{
			tmp[i] = new double[n];
		}
		return tmp;
	}

	void MatrixDelete(double** p, int m)	{
		for (int i=0; i<m; i++)	{
			delete[] p[i];
		}
		delete[] p;
	}

	void MatrixProduct(double** p, double** q, double** r, int l, int m, int n)	{
		for (int i=0; i<l; i++)	{
			for (int j=0; j<n; j++)	{
				double tmp = 0;
				for (int k=0; k<m; k++)	{
					tmp += p[i][k] * q[k][j];
				}
				r[i][j] = tmp;
			}
		}
	}

	void MatrixSum(double** p, double** q, double** r, int m, int n)	{
		for (int i=0; i<m; i++)	{
			for (int j=0; j<n; j++)	{
				r[i][j] = p[i][j] + q[i][j];
			}
		}
	}

	void MatrixAdd(double** p, double** q, int m, int n)	{
		for (int i=0; i<m; i++)	{
			for (int j=0; j<n; j++)	{
				p[i][j] += q[i][j];
			}
		}
	}

	void MatrixSet(double** p, double** q, int m, int n, int offsetm = 0, int offsetn=0)	{
		for (int i=0; i<m; i++)	{
			for (int j=0; j<n; j++)	{
				if (q)	{
					p[i][j] = q[i+offsetm][j+offsetn];
				}
				else	{
					p[i][j] = 0;
				}
			}
		}
	}

	void MatrixSetIdentity(double** p, int m)	{
		for (int i=0; i<m; i++)	{
			for (int j=0; j<m; j++)	{
				p[i][j] = (i==j);
			}
		}
	}

	void MatrixScalarProduct(double** p, double lambda, int m, int n)	{
		for (int i=0; i<m; i++)	{
			for (int j=0; j<n; j++)	{
				p[i][j] *= lambda;
			}
		}
	}
};

#endif

