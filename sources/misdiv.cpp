

#include "Random.h"

int n;
int ntaxa;
const int dim = 2;

double logg(double lambda, double mu, double t)	{
	double ret = 0;
	if (mu < 1e-6)	{
		ret = -lambda*t - 2*log(lambda);
	}
	else if (lambda > mu)	{
		double d = (lambda - mu) * t;
		double e = 0;
		if (d < 250)	{
			e = exp(-d);
		}
		ret = - d - 2 * log(lambda - mu*e);
	}
	else	{
		double d = (mu - lambda) * t;
		double e = 0;
		if (d < 250)	{
			e = exp(-d);
		}
		ret = - d - 2 * log(mu - lambda*e);
	}

	if (isnan(ret))	{
		cerr << "internal log prob is nan\n";
		exit(1);
	}

	return ret;
}

double logprob(double* theta, double* age)	{

	double lambda = theta[0];
	double mu = theta[1];

	double ret = 0;

	double t = age[ntaxa-2];
	if (mu < 1e-6)	{
		ret += (3*ntaxa-4)*log(lambda) - lambda * t;
	}
	else if (lambda > mu)	{
		double d = (lambda - mu) * t;
		double e = 0;
		if (d < 250)	{
			e = exp(-d);
		}
		ret += (ntaxa-1) * log(lambda) - log(mu) + (2*ntaxa-2)*log(lambda - mu) + log(log(lambda / (lambda - mu*e)));
	}
	else	{
		double d = (mu - lambda) * t;
		double e = 0;
		if (d < 250)	{
			e = exp(-d);
		}
		ret += (ntaxa-1) * log(lambda) - log(mu) + (2*ntaxa-2)*log(mu - lambda) + log(log(lambda * e / (mu - lambda*e)));
	}

	for (int j=0; j<ntaxa-1; j++)	{
		ret += logg(lambda,mu,age[j]);
	}

	if (isnan(ret))	{
		cerr << "log prob is nan\n";
		cerr << "lambda : " << lambda << '\n';
		cerr << "mu     : " << mu << '\n';
		exit(1);
	}
	return ret;
}

double meanlogprob(double* theta0, double** theta, double** age, double* logprobs, double& effsize)	{

	double max = 0;
	for (int i=0; i<n; i++)	{
		logprobs[i] = logprob(theta0,age[i]) - logprob(theta[i],age[i]);
		if ((!i) || (max < logprobs[i]))	{
			max = logprobs[i];
		}
	}
	double mean = 0;
	double var = 0;
	for (int i=0; i<n; i++)	{
		double tmp = exp(logprobs[i] - max);
		mean += tmp;
		var += tmp*tmp;
	}
	var /= mean * mean;
	effsize = 1.0 / var;
	mean /= n;
	return log(mean) + max;
}

int main(int argc, char* argv[])	{

	// load file

	ifstream is(argv[1]);
	double delta = atof(argv[2]);
	int every = atoi(argv[3]);
	int until = atoi(argv[4]);
	int fixmu = atoi(argv[5]);

	is >> n >> ntaxa;
	double** theta = new double*[n];
	double** age = new double*[n];
	double* logprobs = new double[n];
	double* meantheta = new double[dim];
	for (int j=0; j<dim; j++)	{
		meantheta[j] = 0;
	}
	for (int i=0; i<n; i++)	{
		theta[i] = new double[dim];
		age[i] = new double[ntaxa-1];
		for (int j=0; j<dim; j++)	{
			is >> theta[i][j];
			if (theta[i][j])	{
				meantheta[j] += log(theta[i][j]);
			}
		}
		for (int j=0; j<ntaxa-1; j++)	{
			is >> age[i][j];
		}
	}
	for (int j=0; j<dim; j++)	{
		meantheta[j] /= n;
		if (meantheta[j])	{
			meantheta[j] = exp(meantheta[j]);
		}
	}


	double* current = new double[dim];
	double* next= new double[dim];
	for (int j=0; j<dim; j++)	{
		current[j] = meantheta[j];
	}
	if (fixmu)	{
		current[1] = 0;
	}

	int nstep = 0;
	int cycle = 0;
	while ((until == -1) || (nstep < until))	{
		double effsize;
		double effsize2;
		double currentlogprob = meanlogprob(current,theta,age,logprobs,effsize);
		for (int j=0; j<dim; j++)	{
			next[j] = current[j] * exp(delta*(Random::Uniform() - 0.5));
		}
		if (fixmu)	{
			next[1] = 0;
		}
		double nextlogprob = meanlogprob(next,theta,age,logprobs,effsize2);
		if (nextlogprob > currentlogprob)	{
			for (int j=0; j<dim; j++)	{
				current[j] = next[j];
			}
			currentlogprob = nextlogprob;
			effsize = effsize2;
		}

		if (!cycle)	{
			cout << nstep;
			cout << '\t' << currentlogprob;
			for (int j=0; j<dim; j++)	{
				cout << '\t' << current[j];
			}
			cout << '\t' << effsize;
			cout << '\n';
		}
		cycle ++;
		if (cycle == every)	{
			cycle = 0;
			nstep++;
		}
	}
}


