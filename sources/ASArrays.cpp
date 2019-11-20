
#include "Random.h"

class ASArrays	{

	double T;
	double dt;
	double da;
	int NT;
	int Na;
	double amin;
	double amax;

	public:

	ASArrays(double inT, double indt, double inamin, double inamax, double inda)	{
		T = inT;
		dt = indt;
		amin = inamin;
		amax = inamax;
		da = inda;
		NT = (int) (T/dt);
		Na = (int) ((amax - amin) / da);
	}

	void Check(double a)	{

		double max = 0;
		for (int j=0; j<=NT; j++)	{
			double t = dt * j;

			double total = 0;
			double kstep = 0.001;
			// double kstep = 0.001 / (1+t);
			// double kstep = 0.001 / (1+exp(a*log(t)));
			double kmax = 1000;
			// double kmax = exp(log(100)/a);
			// cerr << j << '\t' << NT << '\t' << kmax << '\t' << kstep << '\t' << kmax / kstep << '\n';
			int count = 0;
			for (double k = 0; k<=kmax; k+=kstep)	{
				if (!k)	{
					total += 0.5 * cos(k*t) * exp(-exp(a*log(k))/a);
				}
				else	{
					total += cos(k*t) * exp(-exp(a*log(k))/a);
				}
				count++;
			}
			// total /= count;
			total *= kstep;
			total /= Pi;
			double check = 0;
			if (a == 1.0)	{
				check = logCauchyDensity(1.0,t);
			}
			else if (a == 2.0)	{
				check = logNormalDensity(1.0,t);
			}
			double error = fabs(log(total) - check);
			if (max < error)	{
				max = error;
			}
			cerr << t << '\t' << error << '\t' << log(total) << '\t' << check << '\n';
		}
		cerr << '\n';
		cerr << "max error : " << max << '\n';
		exit(1);
	}

	void ComputeArrays(ofstream os)	{

		os << "ASArray\n";
		os << T << '\t' << dt << '\t' << amin << '\t' << amax << '\t' << da << '\n';
		for (int i=0; i<Na; i++)	{
			double a = amin + da * i;
			os << a;
			cerr << i << '\t' << Na << '\n';
			for (int j=0; j<=NT; j++)	{
				double t = dt * j;

				double total = 0;
				double kstep = 0.01 / t;
				double kmax = exp(log(100)/a);
				for (double k = 0; k<=kmax; k+=kstep)	{
					total += cos(k*t) * exp(-exp(a*log(k)));
				}
				total /= Pi;
				if (total <= 0)	{
					cerr << "error in lnpdf: negative density\n";
					cerr << a << '\t' << t << '\t' << total << '\n';
				}
				double logtotal = log(total);
				if (std::isnan(logtotal))	{
					cerr << "error in lnpdf: nan\n";
					exit(1);
				}
				if (std::isinf(logtotal))	{
					cerr << "error in lnpdf: inf\n";
					exit(1);
				}
				os << '\t' << logtotal;
			}
			os << '\n';
		}
	}

	double logNormalDensity(double scale, double x)	{
		return -0.5 * log(2 * Pi) - log(scale) - 0.5 * x * x / scale / scale;
	}

	double logCauchyDensity(double scale, double x)	{
		return -log (Pi * scale * (1 + x*x/scale/scale));
	}


};

int main(int argc, char* argv[])	{

	double T = atof(argv[1]);
	double dt = atof(argv[2]);
	double a = atof(argv[3]);
	double amin = 1.0;
	double amax = 2.0;
	double da = 0.01;

	/*
	double amin = atof(argv[3]);
	double amax = atof(argv[4]);
	double da = atof(argv[5]);
	*/

	ASArrays as(T,dt,amin,amax,da);

	cerr << "check\n";
	as.Check(a);
	cerr << '\n';
}


