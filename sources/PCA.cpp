

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace std;

#include "linalg.h"

class Scatter 	{


	int dim;

	double* mean;
	double** value;
	double** u;
	double** invu;
	double* v;
	double* vi;
	double max1;
	double max2;
	int imax1;
	int imax2;

	public:

	Scatter(int indim, int inN, double** inobs)	{
		dim = indim;
		Create();
		Reset();
		ComputeMean(inobs,inN);
		ComputeScatterMatrix(inobs,inN);
		// Diagonalise();
	}

	int GetDim()const {
		return dim;
	}

	double* GetEigenVal() {
		return v;
	}

	double** GetEigenVect() {
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
			ret += log(GetEigenVal()[i]);
		}
		return ret;
	}

	double** GetInvEigenVect() {
		return invu;
	}

	double* GetMean() {return mean;}
	double** GetVar() {return value;}

	// from of dimension dim
	// to of dimension 2
	void GetPCACoordinate(double* from, double* to)	{
		double* tmp = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			tmp[i] = from[i] - mean[i];
		}
		double* newx = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			double tot = 0;
			for (int j=0; j<GetDim(); j++)	{
				tot += invu[i][j] * tmp[j];
			}
			newx[i] = tot;
		}
		to[0] = newx[imax1];
		to[1] = newx[imax2];
		delete[] tmp;
		delete[] newx;
	}


	// finalmean and finalvar are of dimension 2 and 2x2
	void PCATransform(Scatter* inref, double* finalmean, double** finalvar)	{

		// newmean = inref->invu * (mean - inref->mean);

		double* tmp = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			tmp[i] = mean[i] - inref->mean[i];
		}
		double* newmean = new double[GetDim()];
		double** invU = inref->GetInvEigenVect();
		for (int i=0; i<GetDim(); i++)	{
			double tot = 0;
			for (int j=0; j<GetDim(); j++)	{
				tot += invU[i][j] * tmp[j];
			}
			newmean[i] = tot;
		}
		finalmean[0] = newmean[inref->imax1];
		finalmean[1] = newmean[inref->imax2];

		// newvar = inref->invu * var * inref->u;
		double** U = inref->GetEigenVect();
		double** a = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			a[i] = new double[GetDim()];
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += value[i][k] * U[k][j];
				}
				a[i][j] = tot;
			}
		}
		double** b = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			b[i] = new double[GetDim()];
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += invU[i][k] * a[k][j];
				}
				b[i][j] = tot;
			}
		}

		finalvar[0][0] = b[inref->imax1][inref->imax1];
		finalvar[0][1] = b[inref->imax1][inref->imax2];
		finalvar[1][0] = b[inref->imax2][inref->imax1];
		finalvar[1][1] = b[inref->imax2][inref->imax2];

		for (int i=0; i<GetDim(); i++)	{
			delete[] b[i];
			delete[] a[i];
		}
		delete[] b;
		delete[] a;
		delete[] newmean;
		delete[] tmp;
	}

	private:

	void Create(){
		mean = new double[GetDim()];
		value = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			value[i] = new double[GetDim()];
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
		vi = new double[GetDim()];
	}

	void Reset()	{
		for (int i=0; i<GetDim(); i++)	{
			mean[i] = 0;
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
	}

	void ComputeMean(double** inValues, int Nval)	{

		for (int i=0; i<GetDim(); i++)	{
			mean[i] = 0;
			for (int j=0; j<Nval; j++)	{
				mean[i] += inValues[j][i];
			}
			mean[i] /= Nval;
		}
	}

	void ComputeScatterMatrix(double** inValues, int Nval){
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				value[i][j] = 0;
				for (int l=0; l<Nval; l++) {
					value[i][j] += (inValues[l][i] - mean[i]) * (inValues[l][j] - mean[j]);
				}
				value[i][j] /= (Nval - 1);
			}
		}
	}

	public:

	int Diagonalise() {

		double * w = new double[GetDim()];
		int* iw = new int[GetDim()];
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

		// diagonalise a into v and u
		int failed = EigenRealGeneral(GetDim(), a, v, vi, u, iw, w);
		if (failed)	{
			cerr << "diagonalisation failed\n";
			exit(1);
		}

		max1 = 0;
		max2 = 0;
		imax1 = -1;
		imax2 = -1;
		for (int i=0; i<GetDim(); i++)	{
			if (imax1 == -1)	{
				max1 = v[i];
				imax1 = i;
			}
			else	{
				if (imax2 == -1)	{
					max2 = v[i];
					imax2 = i;
				}
				if (max1 < v[i])	{
					max2 = max1;
					imax2 = imax1;
					max1 = v[i];
					imax1 = i;
				}
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			cerr << v[i] << '\n';
		}
		cerr << '\n';
		cerr << imax1 << '\t' << v[imax1] << '\t' << max1 << '\n';
		cerr << imax2 << '\t' << v[imax2] << '\t' << max2 << '\n';
		cerr << '\n';

		for (int i=0; i<dim; i++)	{
			double total = 0;
			for (int j=0; j<dim; j++)	{
				total += u[j][i] * u[j][i];
			}
			total = sqrt(total);
			for (int j=0; j<dim; j++)	{
				u[j][i] /= total;
			}
		}
		// u-1 = tu
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				invu[j][i] = u[i][j];
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			delete[] a[i];
		}
		delete[] a;
		delete[] w;
		delete[] iw;

		return failed;
	}
};



int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int size;
	int Ntaxa;
	int Nstate;

	is >> size >> Ntaxa >> Nstate;

	double** obs = new double*[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		obs[i] = new double[Nstate];
		for (int j=0; j<Nstate; j++)	{
			is >> obs[i][j];
		}
	}

	double*** pred = new double**[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		pred[i] = new double*[size];
		for (int j=0; j<size; j++)	{
			pred[i][j] = new double[Nstate];
		}
	}

	for (int i=0; i<size; i++)	{
		for (int j=0; j<Ntaxa; j++)	{
			for (int k=0; k<Nstate; k++)	{
				is >> pred[j][i][k];
			}
		}
	}

	cerr << "observed\n";
	Scatter* Obs = new Scatter(Nstate,Ntaxa,obs);
	cerr << "diagonalise\n";
	Obs->Diagonalise();
	cerr << "ok\n\n";

	double** pcaobs = new double*[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		pcaobs[i] = new double[2];
		Obs->GetPCACoordinate(obs[i],pcaobs[i]);
	}

	double** predmean = new double*[Ntaxa];
	double*** predvar = new double**[Ntaxa];
	double*** pcapred = new double**[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		pcapred[i] = new double*[size];
		for (int j=0; j<size; j++)	{
			pcapred[i][j] = new double[2];
			Obs->GetPCACoordinate(pred[i][j],pcapred[i][j]);
		}
	}
	for (int i=0; i<Ntaxa; i++)	{
		predmean[i] = new double[2];
		predvar[i] = new double*[2];
		for (int j=0; j<2; j++)	{
			predvar[i][j] = new double[2];
		}
	}

	cerr << "pred\n";
	Scatter** Pred = new Scatter*[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		cerr << i << '\t' << Ntaxa << '\n';
		Pred[i] = new Scatter(Nstate,size,pred[i]);
		Pred[i]->PCATransform(Obs,predmean[i],predvar[i]);
	}

	// make R script

	string name = argv[1];
	ofstream os((name + ".R").c_str());

	os << " a <- matrix(c(";
	for (int j=0; j<2; j++)	{
		for (int i=0; i<Ntaxa; i++)	{
			os << pcaobs[i][j];
			if (! ((i==Ntaxa-1) && (j==1)))	{
				os << ',';
			}
		}
	}
	os << ")," << Ntaxa << "," << 2 << ")\n";

	os << "plot(a,col=\"red\",lwd=4)\n";

	os << "b <- matrix(c(";
	for (int j=0; j<2; j++)	{
		for (int i=0; i<Ntaxa; i++)	{
			for (int k=0; k<size; k++)	{
				os << pcapred[i][k][j];
				if (! ((i==Ntaxa-1) && (k==size-1) && (j==1)))	{
					os << ',';
				}
			}
		}
	}
	os << ")," << Ntaxa*size << "," << 2 << ")\n";

	os << "points(b,cex=0.2)\n";

	for (int i=0; i<Ntaxa; i++)	{
		for (int k=0; k<size; k++)	{
			os << "segments(" << pcaobs[i][0] << ',' << pcaobs[i][1] << ',' << pcapred[i][k][0] << ',' << pcapred[i][k][1] << ",lwd=0.2,lt=2)\n";
		}
	}
}
