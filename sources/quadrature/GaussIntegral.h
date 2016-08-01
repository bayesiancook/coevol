
#include "InverseWishartMatrix.h"

class GaussIntegral	{

	public:

	GaussIntegral(Var<CovMatrix>* incovar, double (*inf)(double*)) : covar(incovar), f(inf), val(0) {
		covar->SetOrthoNormal(true);
	}

	virtual ~GaussIntegral() {}
	// prepare all fields that will be shared across several numerical integrations based on the same covariance matrix
	// thus, should be called whenever the covariance matrix has changed
	virtual void Update() {}

	void ResetFunction(double (*inf)(double*)) {f = inf;}

	// assumes that all necessary updates in case of a change of the covariance matrix have already been done
	// computes the integral
	virtual void ComputeIntegral() = 0;
	double GetIntegral() {return val;}

	int GetDim() {return covar->GetDim();}

	protected:

	Var<CovMatrix>* covar;
	double (*f)(double*);
	double val;

};

class MultiGauss : public Rvar<RealVector>	{

	public:

	MultiGauss(Var<CovMatrix>* insigma, double (*inf)(double*))	{
		sigma = insigma;
		f = inf;
		beta = 1;
		setval(RealVector(insigma->GetDim()));
		Register(sigma);
	}
	
	~MultiGauss(){
	}

	void SetBeta(double inbeta)	{
		beta = inbeta;
		Corrupt(false);
		Update();
	}

	void drawSample(){
		sigma->drawVal(GetArray());
	}

	double ProposeMove(double tuning){
		if (Random::Uniform() < 0.5)	{
			for ( int i=0; i <GetDim(); i++){
				val()[i] += tuning * (Random::Uniform() - 0.5);
			}
		}
		else	{
			int i = (int) (GetDim() * Random::Uniform());
			val()[i] += tuning * (Random::Uniform() - 0.5);
		}
		return 0;
	}

	double logProb() {
		double tXSX = 0;
		double* dval = GetArray();
		for (int i=0; i<GetDim() ; i++) {
			tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
			for (int j=0; j<i ; j++) {
				tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
			}
		}
		double d = -0.5 * (sigma->GetLogDeterminant() + tXSX);
		d += beta * log((*f)(dval));
		return d;
	}


	protected:

	Var<CovMatrix>* sigma;
	double (*f) (double*);
	double beta;

};

class ImportanceGaussIntegral : public GaussIntegral	{

	public:

	ImportanceGaussIntegral(Var<CovMatrix>* incovar, double (*inf)(double*), int inN) : GaussIntegral(incovar, inf), N(inN) {}

	void ComputeIntegral()	{
		MultiGauss* x = new MultiGauss(covar,f);
		x->SetBeta(1);
		x->Sample();
		covar->FullCorrupt();
		covar->FullUpdate();
		
		double total = 0;

		for (int n=0; n<N; n++)	{
			x->Sample();
			total += (*f)(x->GetArray());
		}
		total /= N;
		val = total;
	}

	private:
	int N;
};

class ThermoGaussIntegral : public GaussIntegral	{

	public:

	ThermoGaussIntegral(Var<CovMatrix>* incovar, double (*inf)(double*), int inN, int inK, double indelta) : GaussIntegral(incovar, inf), N(inN), K(inK), delta(indelta) {}

	void ComputeIntegral()	{
		MultiGauss* x = new MultiGauss(covar,f);
		x->SetBeta(0);
		x->Sample();
		covar->FullCorrupt();
		covar->FullUpdate();
		
		double total = 0.5 * log((*f)(x->GetArray()));
		double success = 0;

		for (int n=1; n<N; n++)	{
			double beta = ((double) n) / N;
			x->SetBeta(beta);
			covar->FullCorrupt();
			covar->FullUpdate();
			for (int k=0; k<K; k++)	{
				success += x->Move(delta);
			}
			total += log((*f)(x->GetArray()));
		}
		x->SetBeta(1.0);
		for (int k=0; k<K; k++)	{
			success += x->Move(delta);
		}
		total += 0.5 * log((*f)(x->GetArray()));

		total /= N;
		val = exp(total);
		cerr << "success rate : " << success / N / K << '\n';
	}

	private:
	int N;
	int K;
	double delta;
};


class SparseGaussIntegral : public GaussIntegral	{

	public:

	SparseGaussIntegral(Var<CovMatrix>* incovar, double (*inf)(double*),int inlevel) : GaussIntegral(incovar, inf), level(inlevel) , N(0), x(0), y(0), w(0) {
		LoadGrid();
	}

	void LoadGrid()	{
		if (x)	{
			for (int i=0; i<N; i++)	{
				delete[] x[i];
				delete[] y[i];
			}
			delete[] x;
			delete[] y;
			delete[] w;
		}

		ostringstream xfile;
		xfile <<  "herm_d20_level" << level << "_x.txt";
		ifstream xis(xfile.str().c_str());

		ostringstream wfile;
		wfile <<  "herm_d20_level" << level << "_w.txt";
		ifstream wis(wfile.str().c_str());
	
		int Nstate;
		xis >> N >> Nstate;
		if (Nstate != GetDim())	{
			cerr << "error in SparseGaussIntegral: non matching dimension\n";
			exit(1);
		}
		cerr << N << '\t' << Nstate << '\n';

		x = new double*[N];
		y = new double*[N];
		w = new double[N];
		for (int i=0; i<N; i++)	{
			x[i] = new double[GetDim()];
			y[i] = new double[GetDim()];
			for (int k=0; k<GetDim(); k++)	{
				xis >> x[i][k];
			}
			wis >> w[i];
		}

		cerr << "grid loaded\n";
	}


	void Update()	{
		double* aux = new double[GetDim()];
		double** u = covar->GetEigenVect();
		// double** invu = covar->GetInvEigenVect();
		double* diag = covar->GetEigenVal();
		for (int i=0; i<N; i++)	{
			// cerr << w[i] << '\t';
			for (int j=0; j<GetDim(); j++)	{
				// aux[j] = x[i][j] / sqrt(diag[j]);
				// aux[j] = x[i][j] * 2 / sqrt(diag[j]);
				aux[j] = x[i][j] * sqrt(diag[j]) * sqrt(2);
			}
			for (int j=0; j<GetDim(); j++)	{
				double total = 0;
				for (int k=0; k<GetDim(); k++)	{
					total += u[j][k] * aux[k];
					// total += invu[j][k] * aux[k];
				}
				y[i][j] = total;
				// cerr << total << '\t';
			}
			// cerr << '\n';
			/*
			for (int j=0; j<GetDim(); j++)	{
				y[i][j] = x[i][j];
			}
			*/
		}
		delete[] aux;
	}

	void ComputeIntegral()	{
		Update();
		double total = 0;
		for (int i=0; i<N; i++)	{
			total += w[i] * (*f)(y[i]);
			// cerr << w[i] << '\t' << total << '\n';
		}
		val = total / exp(20 * log(sqrt(3.1415926535)));
	}

	
	private:

	int level;
	int N;
	double** x;
	double** y;
	double* w;
};

