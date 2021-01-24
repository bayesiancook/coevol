#ifndef COVMATRIX_H
#define COVMATRIX_H

#include "linalg.h"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

using namespace std;

#include "RandomTypes.h"
#include "ValArray.h"

class CovMatrix : public BaseType, public Multiplicative {

	protected:


	// data members

	int dim;

	// value : the infinitesimal generator matrix
	double ** value;

	bool orthonormal;
	bool diagflag;


	static int maxn;
	static int ndiag;
	static double maxerror;
	static double toterror;

	public:

	static double GetMeanDiagError()	{
		if (ndiag)	{
			return toterror / ndiag;
		}
		return 0;
	}

	static double GetMaxDiagError()	{
		return maxerror;
	}

	static double GetMaxDiagN()	{
		return maxn;
	}

	CovMatrix() : dim(0) , value(0) , orthonormal(false)  {}

	CovMatrix(int indim) : orthonormal(false) {
		dim = indim;
		CreateAllMatrix();
		diagflag = false;
	}

	CovMatrix(double** Q, int indim) : orthonormal(false) {
		dim = indim;
		CreateAllMatrix();
		for (int i=0; i<GetDim(); i++)	{
			value[i][i] = Q[i][i];
			for (int j=0; j<i; j++)	{
				value[i][j] =  Q[i][j];
				value[j][i] = value[i][j];
			}
		}
		diagflag = false;
	}

	CovMatrix(const CovMatrix& from) : orthonormal(false) {
		dim = from.GetDim();
		CreateAllMatrix();
		diagflag = false;
		*this = from;
	}

	void SetOrthoNormal(bool in)	{
		orthonormal = in;
	}

	void CreateAllMatrix(){
		value = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			value[i] = new double[GetDim()];
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}

		invvalue = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			invvalue[i] = new double[GetDim()];
		}

		u = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			u[i] = new double[GetDim()];
		}

		invu = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			invu[i] = new double[GetDim()];
		}

		v = new double[GetDim()];
		logv = new double[GetDim()];
	}


	~CovMatrix() {

		for (int i=0; i<GetDim(); i++)	{
			delete[] value[i];
			delete[] invvalue[i];
			delete[] u[i];
			delete[] invu[i];
		}
		delete[] value;
		delete[] invvalue;
		delete[] u;
		delete[] invu;

		delete[] v;
		delete[] logv;
	}




	CovMatrix&	operator=(const CovMatrix& from)	{
		if (!dim)	{
			dim = from.GetDim();
			CreateAllMatrix();
			diagflag = false;
			}
		if (dim != from.dim)	{
			cerr << "error : non matching dimenstion for matrix\n";
			cerr << GetDim() << '\t' << from.GetDim() << '\n';
			exit(1);
		}
		for (int i=0; i<GetDim(); i++){
			value[i][i] =from.value[i][i];
			for  (int j=0; j<i; j++){
				value[i][j] =from.value[i][j];
				value[j][i] =from.value[i][j];
			}
		}
		diagflag = false;
		return *this;
	}

	void setVal(double** val, int inDim) {
		if (!dim)	{
			dim = inDim;
			CreateAllMatrix();
			diagflag = false;
			}
		if (dim != inDim)	{
			cerr << "error : non matching dimenstion for matrix\n";
			cerr << GetDim() << '\t' << inDim << '\n';
			exit(1);
		}
		for (int i=0; i<GetDim(); i++){
			value[i][i] =val[i][i];
			for  (int j=0; j<i; j++){
				value[i][j] =val[i][j];
				value[j][i] =val[i][j];
			}
		}
		diagflag = false;

	}

	/*
	int GetNVals()	{
		return dim*(dim+1)/2;
	}

	void GetVals(double* array, int& index)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i; j<GetDim(); j++)	{
				array[index++] = val[i][j];
			}
		}
	}

	void SetVals(double* array, int& index)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i; j<GetDim(); j++)	{
				val[i][j] = val[j][i] = array[index++];
			}
		}
	}
	*/

	CovMatrix&	operator+=(const CovMatrix from)	{
				for (int i=0; i<GetDim(); i++){
					value[i][i] += from.value[i][i];
					for  (int j=0; j<i; j++){
						value[i][j] += from.value[i][j];
						value[j][i] += from.value[i][j];
					}
				}
				return *this;
				diagflag = false;
			}

	double*		operator[](int i) {
				double* v = value[i];
				return v;
			}

	double	operator()(int i, int j) const {
		return value[i][j];
	}

	CovMatrix& operator*=(double scal) {
		for (int i=0; i<GetDim(); i++){
			value[i][i] *= scal;
			for  (int j=0; j<i; j++){
				value[i][j] *= scal;
				value[j][i] *= scal;
			}
		}
		diagflag = false;
		return *this;
	}

	CovMatrix operator*(double scal) {

		CovMatrix cm = *this;
		cm*=scal;
		return cm;
	}


	virtual double	ProposeMove(double tuning)	{
		int j = (int) (GetDim() * (GetDim()+1) * 0.5 * Random::Uniform());
		int i=0;
		while(i+1<=j){
			i++;
			j -= i;
		}

		if (i == j)	{
			double infbound = 0;
			for (int k=0; k<GetDim(); k++)	{
				if (k != i)	{
					double tmp = value[i][k] * value[i][k] / value[k][k];
					if (infbound < tmp)	{
						infbound = tmp;
					}
				}
			}
			double m = tuning * (Random::Uniform() - 0.5);
			value[i][i] += m;
			if (value[i][i] < infbound)	{
				value[i][i] = 2 * infbound - value[i][i];
			}
		}
		else	{
			double bound = sqrt(value[i][i] * value[j][j]);
			double m = tuning * (Random::Uniform() - 0.5);
			value[i][j] += m;
			while ((value[i][j] < -bound) || (value[i][j] > bound))	{
				if (value[i][j] < -bound)	{
					value[i][j] = -2 * bound - value[i][j];
				}
				if (value[i][j] > bound)	{
					value[i][j] = 2 * bound - value[i][j];
				}
			}
			value[j][i] = value[i][j];
		}
		return 0;
	}

	/*
	virtual double	ProposeMove(double tuning)	{
		int j = (int) (GetDim() * (GetDim()+1) * 0.5 * Random::Uniform());
		int i=0;
		while(i+1<=j){
			i++;
			j -= i;
		}
		double m = tuning * (Random::Uniform() - 0.5);
		// a diagonal value (variance)
		if(i == j){
			double borneInf = -400;
			double c = 0;
			for (i=0; i< GetDim() ; i++) {
				if(i!=j){
					if(i<j){
						c = value[i][j];
					}
					if(j<i){
						c = value[i][j];
					}
					c = abs( c / sqrt(value[i][i]));
					if(c > borneInf){
						borneInf = c;
					}
				}
			}
			i = j;
			double d = value[i][j] - m - borneInf;
			if(d < 0){
				m -= 2 * d;
			}
		}
		// a covariance
		else{
			double borne = sqrt(value[i][i]) * sqrt(value[j][j]) - abs(value[i][j]);
			while (abs(m) > borne){
				if (m < -borne){
					m = - 2 * borne - m;
				}
				if (m > borne){
					m = 2 * borne - m;
				}
			}
		}
		value[j][i] += m;
		value[i][j] = value[j][i];
		return 0;
	}
	*/

	virtual void drawVal(double* vec){
		double* randomvalues = new double[GetDim()];
		for (int i=0; i<GetDim(); i++) {
			randomvalues[i] =  Random::sNormal() * sqrt(GetEigenVal()[i]);
			vec[i]=0;
		}
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				vec[j] +=  randomvalues[i] * GetEigenVect()[j][i];
			}
		}
		delete[] randomvalues;
	}

	void drawValInv(double* vec){
		double* randomvalues = new double[GetDim()];
		for (int i=0; i<GetDim(); i++) {
			randomvalues[i] =  Random::sNormal() / sqrt(GetEigenVal()[i]);
			vec[i]=0;
		}
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				vec[j] +=  randomvalues[i] * GetEigenVect()[j][i];
			}
		}
		delete[] randomvalues;
	}

	double logValProb(double* dval)	{

		double tXSX = 0;
		for (int i=0; i<GetDim() ; i++) {
			tXSX += GetInvMatrix()[i][i] * dval[i] * dval[i];
			for (int j=0; j<i ; j++) {
				tXSX += 2 * GetInvMatrix()[i][j] * dval[j] * dval[i];
			}
		}
		return -0.5 * (GetLogDeterminant() + tXSX);
	}

	private:


	double ** invvalue;

	// v : eigenvalues
	// u : the matrix of eigen vectors
	// invu : the inverse of u

	double ** u;
	double ** invu;
	double * v;
	double * logv;

	public:


	double* GetEigenVal() {
		if (! diagflag) Diagonalise();
		return v;
	}

	double* GetLogEigenVal() {
		if (! diagflag) Diagonalise();
		return logv;
	}

	double** GetEigenVect() {
		if (! diagflag)	{
			Diagonalise();
		}

		return u;
	}

	double GetDeterminant() {
		double ret = 1;
		for (int i=0; i<GetDim(); i++) {
			ret *= GetEigenVal()[i];
		}
		return ret;
	}

	double GetLogDeterminant() {
		double ret = 0;
		for (int i=0; i<GetDim(); i++) {
			ret += GetLogEigenVal()[i];
		}
		/*
		if (isnan(ret))	{
			cerr << "covmatrix det: nan\n";
			exit(1);
		}
		*/
		return ret;
	}

	double** GetInvEigenVect() {
		if (! diagflag) Diagonalise();
		return invu;
	}

	int GetDim()const {
		return dim;
	}

	void Reset()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
	}

	void SetIdentity()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			value[i][i] = 1;
		}
	}

	void Project(int index, double** m)	{

		int k = 0;
		for (int i=0; i<GetDim(); i++)	{
			if (i != index)	{
				int l = 0;
				for (int j=0; j<GetDim(); j++)	{
					if (j != index)	{
						m[k][l] = value[i][j] - value[i][index] * value[j][index] / value[index][index];
						l++;
					}
				}
				k++;
			}
		}
	}


	/*
	void ProjectOnFirstComponent(double* sol)	{

		// sol is a vector of dimension GetDim()
		// first component is residual error
		// all other are correlation coefficients

		double* w = new double[GetDim() - 1];
		int* iw = new int[GetDim() - 1];
		double** a = new double*[GetDim() - 1];
		double** b = new double*[GetDim() - 1];
		for (int i=0; i<GetDim()-1; i++)	{
			a[i] = new double[GetDim()-1];
			b[i] = new double[GetDim()-1];
			for (int j=0; j<GetDim()-1; j++)	{
				a[i][j] = value[i+1][j+1];
			}
		}

		InvertMatrix(a, GetDim()-1, w, iw, b);

		// apply b on vector

		for (int i=0; i<GetDim()-1; i++)	{
			double tot = 0;
			for (int j=0; j<GetDim()-1; j++)	{
				tot += b[i][j] * value[0][j+1];
			}
			sol[i+1] = tot;
		}

		double tot = value[0][0];
		for (int i=1; i<GetDim(); i++)	{
			tot += sol[i] * sol[i] * value[i][i];
			tot -= 2 * sol[i] * value[0][i];
		}
		for (int i=1; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += 2 * sol[i] * sol[j] * value[i][j];
			}
		}
		sol[0] = tot;

		delete[] w;
		delete[] iw;
		for (int i=0; i<GetDim()-1; i++)	{
			delete[] b[i];
			delete[] a[i];
		}
		delete[] b;
		delete[] a;
	}
	*/

	/*
	void Add(const CovMatrix& from)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] += from.value[i][j];
			}
		}
	}
	*/

	int ScalarMultiplication(double scal)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] *= scal;
			}
		}
		diagflag = false;
		return GetDim();
	}

	double** GetMatrix() {
		return value;
	}

	double** GetMatrix() const{
		return value;
	}

	double** GetInvMatrix() {
		if (! diagflag) Diagonalise();
		return invvalue;
	}


	bool isPositiveDefine(){
		if (! diagflag) Diagonalise();
		bool r = true;
		for (int i=0; i<GetDim(); i++){
			if(GetEigenVal()[i] <= 1e-6){
				r = false;
			}
		}
		return r;
	}

	void CorruptDiag()	{
		diagflag = false;
	}

	double GetMax()	{
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = fabs(value[i][j]);
				if (max < tmp)	{
					max = tmp;
				}
			}
		}
		return max;
	}

	//Set the matrix to it s inverse //loook si diagflag
	int Invert(){
		double** a = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			a[i] = new double[GetDim()];
		}

		// copy value into a :
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				a[i][j] = value[i][j];
			}
		}

		// invert a into value
		// InvertMatrix(a, GetDim(), w, iw, value);
		double logdet = LinAlg::Gauss(a,GetDim(),value);

		// cerr << "check inverse : " << CheckInverse() << '\n';
		for (int i=0; i<GetDim(); i++)	{
			delete[] a[i];
		}
		delete[] a;
		diagflag = false;
		if (std::isinf(logdet))	{
			cerr << "error in cov matrix: non invertible\n";
			return 1;
			exit(1);
		}
		return 0;
	}

	double CheckInverse()	{
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += value[i][k] * GetInvMatrix()[k][j];
				}
				if (i == j)	{
					tot --;
				}
				if (max < fabs(tot))	{
					max = fabs(tot);
				}
			}
		}
		return max;
	}


	int Diagonalise()	{
		int nmax = 1000;
		double epsilon = 1e-10;

		int n = LinAlg::DiagonalizeSymmetricMatrix(value,dim,nmax,epsilon,v,u);
		if (maxn < n)	{
			maxn = n;
		}
		bool failed = (n == nmax);
		if (failed)	{
			cerr << "diag failed\n";
			cerr << n << '\n';
			for (int i=0; i<dim; i++)	{
				for (int j=0; j<dim; j++)	{
					cerr << value[i][j] << '\t';
				}
				cerr << '\n';
			}
			exit(1);
		}

		// normalise u
		for (int i=0; i<dim; i++)	{
			double total = 0;
			for (int j=0; j<dim; j++)	{
				total += u[j][i] * u[j][i];
			}
		}
		// u-1 = tu
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				invu[j][i] = u[i][j];
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			logv[i] = log(v[i]);
		}

		LinAlg::Gauss(value,GetDim(),invvalue);

		diagflag = true;
		double tmp = CheckDiag();
		if (maxerror < tmp)	{
			maxerror = tmp;
		}
		toterror += tmp;
		ndiag ++;
		return failed;
	}

	/*
	int Diagonalise()	{

		int nmax = 1000;
		double epsilon = 1e-10;
		bool failed = false;
		double eps = 1e-8;

		int effdim = 0;
		int index[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			index[i] = 0;
			for (int j=0; j<GetDim(); j++)	{
				if (fabs(value[i][j]) > eps)	{
					index[i] = 1;
				}
			}
			if (index[i])	{
				effdim++;
			}
		}
		if (effdim < GetDim())	{

			cerr << "effdim < dim in covmatrix\n";
			exit(1);

			double** a = new double*[effdim];
			double** inva = new double*[effdim];
			double** b = new double*[effdim];
			double* w = new double[effdim];
			for (int i=0; i<effdim; i++)	{
				a[i] = new double[effdim];
				inva[i] = new double[effdim];
				b[i] = new double[effdim];
			}
			int k = 0;
			for (int i=0; i<GetDim(); i++)	{
				if (index[i])	{
					int l = 0;
					for (int j=0; j<GetDim(); j++)	{
						if (index[j])	{
							a[k][l] = value[i][j];
							l++;
						}
					}
					k++;
				}
			}

			int n = LinAlg::DiagonalizeSymmetricMatrix(a,effdim,nmax,epsilon,w,b);

			if (maxn < n)	{
				maxn = n;
			}
			failed = (n == nmax);
			if (failed)	{
				cerr << "diag failed\n";
				cerr << "in reduced code : " << effdim << '\n';
				cerr << (*this) << '\n';
				// oStream(cerr);
				exit(1);
			}

			LinAlg::Gauss(a,effdim,inva);

			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					invvalue[i][j] = 0;
					u[i][j] = 0;
				}
			}
			for (int i=0; i<GetDim(); i++)	{
				u[i][i] = 1;
			}
			k = 0;
			for (int i=0; i<GetDim(); i++)	{
				if (index[i])	{
					int l = 0;
					for (int j=0; j<GetDim(); j++)	{
						if (index[j])	{
							invvalue[i][j] = inva[k][l];
							u[i][j] = b[k][l];
							l++;
						}
					}
					k++;
				}
			}
			k = 0;
			for (int i=0; i<GetDim(); i++)	{
				if (index[i])	{
					v[i] = w[k];
					k++;
				}
			}

			for (int i=0; i<effdim; i++)	{
				delete[] a[i];
				delete[] inva[i];
				delete[] b[i];
			}
			delete[] a;
			delete[] inva;
			delete[] b;
			delete[] w;
		}
		else	{

			int n = LinAlg::DiagonalizeSymmetricMatrix(value,dim,nmax,epsilon,v,u);
			if (maxn < n)	{
				maxn = n;
			}
			failed = (n == nmax);
			if (failed)	{
				return failed;
				cerr << "diag failed\n";
				cerr << (*this) << '\n';
				// oStream(cerr);
				exit(1);
			}


			LinAlg::Gauss(value,GetDim(),invvalue);
		}

		// normalise u
		for (int i=0; i<dim; i++)	{
			double total = 0;
			for (int j=0; j<dim; j++)	{
				total += u[j][i] * u[j][i];
			}
			if (fabs(total - 1) > 1e-7)	{
				cerr << "in CovMatrix::Diagonalize: normalization error\n";
				cerr << total << '\n';
				exit(1);
			}
		}
		// u-1 = tu
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				invu[j][i] = u[i][j];
			}
		}

		// double logdet = 0;
		for (int i=0; i<GetDim(); i++)	{
			logv[i] = log(v[i]);
			// logdet += logv[i];
		}

		diagflag = true;
		double tmp = CheckDiag();
		if (maxerror < tmp)	{
			maxerror = tmp;
		}
		toterror += tmp;
		ndiag ++;
		return failed;
	}
	*/

	double CheckDiag()	{
		double** a = new double*[GetDim()];
		double** b = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			a[i] = new double[GetDim()];
			b[i] = new double[GetDim()];
		}

		/*
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				cerr << GetMatrix()[i][j] << '\t';
			}
			cerr << '\n';
		}
		cerr << '\n';

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				cerr << u[i][j] << '\t';
			}
			cerr << '\n';
		}
		cerr << '\n';

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				cerr << invu[i][j] << '\t';
			}
			cerr << '\n';
		}
		cerr << '\n';
		*/

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += invu[i][k] * GetMatrix()[k][j];
				}
				a[i][j] = tot;
			}
		}

		double max = 0;

	// 	cerr << "check diag\n";
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += a[i][k] * u[k][j];
				}
				b[i][j] = tot;
				if (i != j)	{
					if (max < fabs(tot))	{
						max = fabs(tot);
					}
				}
	//			cerr << b[i][j] << '\t';
			}
	//		cerr << '\n';
		}

		/*
		cerr << "check inverse\n";
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += invu[i][k] * u[k][j];
				}
				a[i][j] = tot;
				cerr << a[i][j] << '\t';
			}
			cerr << '\n';
		}
		*/

		for (int i=0; i<GetDim(); i++)	{
			delete[] a[i];
			delete[] b[i];
		}
		delete[] a;
		delete[] b;

		return max;
	}

	/*
	double CheckInverse()	{
		double** a = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			a[i] = new double[GetDim()];
		}
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += GetInvMatrix()[i][k] * GetMatrix()[k][j];
				}
				a[i][j] = tot;
				if (max < fabs(tot))	{
					max = fabs(tot);
				}
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			delete[] a[i];
		}
		delete[] a;
		return max;
	}
	*/

	void ScatterizeOnZero(double** inValues, int Nval){
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				value[i][j] = 0;
				for (int l=0; l<Nval; l++) {
					value[i][j] += inValues[l][i] * inValues[l][j];
				}
			}
		}
		diagflag = false;
	}

	void ScatterizeOnZero(VarArray<RealVector>* inValues){
		int Nval = inValues->GetSize();
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				value[i][j] = 0;
				for (int l=0; l<Nval; l++) {
					value[i][j] += (*inValues->GetVal(l))[i] * (*inValues->GetVal(l))[j];
				}
			}
		}
		diagflag = false;
	}

	void PrintCorrelationCoefficients(ostream& os)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i==j) 	{
					os << "\t-";
				}
				else	{
					os << '\t' << GetMatrix()[i][j] / sqrt(GetMatrix()[i][i] * GetMatrix()[j][j]);
				}
			}
			os << '\n';
		}
	}

	void PrintEigenVectors(ostream& os)	{
		os << "val";
		for (int j=0; j<GetDim(); j++)	{
			os << '\t' << j;
		}
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			os << v[i] << '\t';
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << u[i][j];
			}
			os << '\n';
		}
		os << '\n';
		os << "inverse eigenvector matrix\n";
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << invu[i][j];
			}
			os << '\n';
		}
		os << '\n';

		// proportion of variance explained
		double** prop = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			prop[i] = new double[GetDim()];
		}
		for (int i=0; i<GetDim(); i++)	{
			double total = 0;
			for (int j=0; j<GetDim(); j++)	{
				prop[i][j] = invu[j][i] * invu[j][i] * v[j];
				// prop[i][j] = u[i][j] * u[i][j] * v[j];
				// prop[j][i] = invu[j][i] * invu[j][i] * v[j];
				total += prop[i][j];
			}
			for (int j=0; j<GetDim(); j++)	{
				prop[i][j] /= total;
			}
		}
		os << "proportion of variance explained\n";
		for (int i=0; i<GetDim(); i++)	{
			os << i;
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << prop[i][j];
			}
			os << '\n';
		}
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			delete[] prop[i];
		}
		delete[] prop;

	}

	friend ostream& operator<<(ostream& os, const CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				os << '\t' << r.GetMatrix()[i][j];
			}
			os << '\n';
		}
		return os;
	}

	friend istream& operator>>(istream& is, CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				is >> r.value[i][j];
			}
		}
		return is;
	}

};

#endif
