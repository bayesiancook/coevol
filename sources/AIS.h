
#ifndef AIS_H
#define AIS_H

class AISController : public Mnode {

	public:

	public:

	AISController() : nais(0) {}

	void AISRegister(Var<PosReal>* p, double li, double lf)	{
		param.push_back(p);
		p->Register(this);
		loginit.push_back(li);
		logfinal.push_back(lf);
		nais++;
	}

	double AISSet(double p)	{
		Corrupt(true);
		for (int i=0; i<nais; i++)	{
			param[i]->setval((1-p)*loginit[i] + p*logfinal[i]);
			// param[i]->setval(exp(p*loginit[i] + (1-p)*logfinal[i]));
		}
		double logratio = Update();
		return logratio;
	}

	protected:

	int nais;
	vector<double> loginit;
	vector<double> logfinal;
	vector<Var<PosReal>*> param;

};


#endif
