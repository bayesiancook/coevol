#include "BaseType.h"

const double Profile::MIN = 1e-20;

double	Profile::ProposeMove(double tuning, int n)	{ // n==0dirichlet resampling, otherwise, vase communiquants
	double ret = 0;
	if (!n)	{ // dirichlet
		double* oldprofile = new double[dim];
		for (int i=0; i<dim; i++)	{
			oldprofile[i] = profile[i];
		}
		double total = 0;
		for (int i=0; i<dim; i++)	{
			profile[i] = Random::sGamma(tuning*oldprofile[i]);
			if (profile[i] == 0)	{
				cerr << "error in dirichlet resampling : 0 \n";
				exit(1);
			}
			total += profile[i];
		}

		double logHastings = 0;
		for (int i=0; i<dim; i++)	{
			profile[i] /= total;

			logHastings += - Random::logGamma(tuning*oldprofile[i]) + Random::logGamma(tuning*profile[i])
						-  (tuning*profile[i] -1.0) * log(oldprofile[i]) + (tuning * oldprofile[i] -1.0) * log(profile[i]);
		}

		delete[] oldprofile;
		return logHastings;
	}
	else	{
		if (2*n > dim)	{
			n = dim / 2;
		}
		int* indices = new int[2*n];
		Random::DrawFromUrn(indices,2*n,dim);
		for (int i=0; i<n; i++)	{
			int i1 = indices[2*i];
			int i2 = indices[2*i+1];
			double tot = profile[i1] + profile[i2];
			double x = profile[i1];

			// double h = tuning * (Random::Uniform() - 0.5);
			double h = tot * tuning * (Random::Uniform() - 0.5);
			/*
			int c = (int) (h / (2 * tot));
			h -= c*2*tot;
			*/
			x += h;
			while ((x<0) || (x>tot))	{
				if (x<0)	{
					x = -x;
				}
				if (x>tot)	{
					x = 2*tot - x;
				}
			}
			profile[i1] = x;
			profile[i2] = tot - x;
		}
		delete[] indices;
	}
	return ret;
}

