#include "Random.h"

int main(int argc, char* argv[])	{

	int N = atoi(argv[1]);
	double theta = atof(argv[2]);
	double logprob[N+1];
	double prob[N+1];
	double x = 1;
	double total = 0;
	double logfact[N+1];
	double l = 0;
	logfact[0] = 0;
	for (int i=1; i<=N; i++)	{
		l += log(i);
		logfact[i] = l;
	}
	for (int i=0; i<=N; i++)	{
		logprob[i] = (N-i) * log(1-theta) + i*log(theta) + logfact[N] - logfact[i] - logfact[N-i];
	//	cerr << logprob[i] << '\t' << exp(logprob[i]) << '\n';
		prob[i] = exp(logprob[i]);
		total += prob[i];
	}

	cerr << "total prob : " << total << '\n';

	int i = -1;
	double inf = 0;
	while (inf < 0.025)	{
		i++;
		inf += prob[i];
	}
	cerr << "min: " << i-1 << '\n';

	i = N+1;
	double sup = 0;
	while (sup < 0.025)	{
		i--;
		sup += prob[i];
	}

	cerr << "max: " << i+1 << '\n';

}

