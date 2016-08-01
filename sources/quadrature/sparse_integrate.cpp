#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "linalg.h"

using namespace std;


int main(int argc, char* argv[])	{

	ifstream covis(argv[1]);

	int Nstate = 20;
	/*
	int level = atoi(argv[2]);
	int Nstate;
	int N;
	ostringstream xfile = "herm_d20_level" << level << "_x.txt";
	ifstream xis(xfile.str().c_str());

	ostringstream wfile = "herm_d20_level" << level << "_w.txt";
	ifstream wis(wfile.str().c_str());

	xis >> N >> Nstate;
	cerr << N << '\t' << Nstate << '\n';

	double** x = new double*[N];
	double* p = new double[N];
	for (int i=0; i<N; i++)	{
		for (int k=0; k<Nstate; k++)	{
			xis >> x[i][k];
		}
		wis >> p[i];
	}
	*/

	double* mean = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		covis >> mean[k];
	}
	
	double** covar = new double*[Nstate];
	for (int k=0; k<Nstate; k++)	{
		covar[k] = new double[Nstate];
		for (int l=0; l<Nstate; l++)	{
			covis >> covar[k][l];
			cerr << covar[k][l] << '\t';
		}
		cerr << '\n';
	}

	double* v = new double[Nstate];
	double* vi = new double[Nstate];
	double* v2 = new double[Nstate];

	double** u = new double*[Nstate];
	double** invu = new double*[Nstate];
	double** u2 = new double*[Nstate];
	for (int k=0; k<Nstate; k++)	{
		u[k] = new double[Nstate];
		invu[k] = new double[Nstate];
		u2[k] = new double[Nstate];
	}

	double * w = new double[Nstate];
	int* iw = new int[Nstate];
	double** a = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		a[i] = new double[Nstate];
	}

	// copy covar into a :
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			a[i][j] = covar[i][j];
		}
	}

	// diagonalise a into v and u
	EigenRealGeneral(Nstate, a, v, vi, u, iw, w);

	for (int i=0; i<Nstate; i++)	{
		cerr << v[i] << '\t' << vi[i] << '\n';
	}

	// normalise u
	for (int i=0; i<Nstate; i++)	{
		double total = 0;
		for (int j=0; j<Nstate; j++)	{
			total += u[j][i] * u[j][i];
		}
		total = sqrt(total);
		for (int j=0; j<Nstate; j++)	{
			u[j][i] /= total;
		}
	}
	/*
	// check
	// copy u into a :
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += covar[i][k] * u[k][j];
			}
			a[i][j] = total;
		}
	}
	for (int i=0; i<Nstate; i++)	{
		for (int j=0; j<Nstate; j++)	{
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += u[k][i] * a[k][j];
			}
			cerr << ((double) ((int) (1000* total))) / 1000 << '\t';
		}
		cerr << '\n';
	}
	*/
	

	for (int i=0; i<Nstate; i++)	{
		delete[] a[i];
	}
	delete[] a;
	delete[] w;
	delete[] iw;

	
}

