#include "GaussIntegral.h"

const int Nstate = 20;

double f(double* x)	{

	double s = (x[2] - x[1]);
	if (s > 100)	{
		return  s * exp(-s);
	}
	if (s < -100)	{
		return  s;
	}

	if (fabs(s) < 1e-8)	{
		return  1;
	}
		
	return s / (exp(s) - 1);

};


int main(int argc, char* argv[])	{

	int N = atoi(argv[1]);
	int K = atoi(argv[2]);
	double delta = atof(argv[3]);
	int level = atoi(argv[4]);
	ifstream covis(argv[5]);

	double* mean = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		covis >> mean[k];
	}
	
	double** covar = new double*[Nstate];
	for (int k=0; k<Nstate; k++)	{
		covar[k] = new double[Nstate];
		for (int l=0; l<Nstate; l++)	{
			covis >> covar[k][l];
			// covar[k][l] = (k==l) ? 1 : 0 ;
			covar[k][l] /= 10;
		}
	}

	Const<CovMatrix>* Covar = new Const<CovMatrix>(CovMatrix(covar,Nstate));
	cerr << "logdet  " << Covar->GetLogDeterminant() << '\n';
	Covar->SetOrthoNormal(true);

	ThermoGaussIntegral* thermo = new ThermoGaussIntegral(Covar,&f,N,K,delta);
	thermo->ComputeIntegral();
	cout << "integral is     : " << thermo->GetIntegral() << '\n';

	ImportanceGaussIntegral* is = new ImportanceGaussIntegral(Covar,&f,N);
	is->ComputeIntegral();
	cout << "integral is     : " << is->GetIntegral() << '\n';

	SparseGaussIntegral* sparse = new SparseGaussIntegral(Covar,&f,level);
	sparse->ComputeIntegral();
	cout << "sparse estimate : " << sparse->GetIntegral() << '\n';
}
