

// still to be implemented: relations between log kappa1, log kappa2 and mean s
    /*
	double NeutralityIndexFactor(int nind, int jmax)	{


        double denom = 0;
        for (int k=1; k<2*nind; k++)    {
            denom += 1.0 / k;
        }

		double total = 0;
        for (int j=2; j<=jmax; j++)	{
            double z = boost::math::zeta(j) / j;
            double num = 0;
            for (int k=1; k<2*nind; k++)    {
                num += exp(Random::logGamma(2*nind+1) - Random::logGamma(k+1) - Random::logGamma(2*nind-k+1) + Random::logGamma(k) + Random::logGamma(2*nind-k+j) - Random::logGamma(2*nind+j));
            }
            total += z * num / denom;
        }
		cerr << "NI factor : " << total << '\n';
		cerr << "NI(0.2)   : " << 1 + 0.2*total << '\n';
		return total;
	}
    */

    /*
    int i(2);
    double k(0);
    double tempo(0);
    do {
        tempo = (boost::math::zeta (i) * factorielle(1999) * factorielle(i) * factorielle(2001)) / (factorielle(2000 + i) * factorielle(1999) * i);
        k += tempo;
        i++;
    }while (tempo > pow(10, -7));

    K = new Const<PosReal>(k);
    */

class ConstrainedLogKappa2 : public Dvar<Real>   {

    public:

    ConstrainedLogKappa2(Var<Real>* inbeta, Var<Real>* inlogkappa1) {
        beta = inbeta;
        logkappa1 = inlogkappa1;
        Register(beta);
        Register(logkappa1);
        specialUpdate();
    }

    void specialUpdate()    {
        cerr << "in constrained log kappa2\n";
        cerr << "still to be implemented\n";
        exit(1);
    }

    protected:
    Var<Real>* beta;
    Var<Real>* logkappa1;
};
