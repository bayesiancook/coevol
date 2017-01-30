#ifndef GLOBALSCALE_H
#define GLOBALSCALE_H

class GlobalScalingFunction : public Dvar<void> {

	public:

	GlobalScalingFunction() {}
	virtual ~GlobalScalingFunction() {}

	virtual double GetScalingFactor(double time1, double time2) = 0;
};


class BurstScalingFunction : public GlobalScalingFunction	{

	public:

	BurstScalingFunction(Var<PosReal>* int0, Var<PosReal>* infactor, Var<PosReal>* inrate)	{

		t0 = int0;
		factor = infactor;
		rate = inrate;
		Register(t0);
		Register(factor);
		Register(rate);
	}

	double GetScalingFactor(double time1, double time2)	{

		if (time1 < time2)	{
			cerr << "error in getscalingfactor: time1 < time2 : " << time1 << '\t' << time2 << '\n';
			exit(1);
		}

		double v = time1 - time2;
		if (time2 < t0->val())	{
			double tmax = time1 - t0->val();
			if (time1 > t0->val())	{
				tmax = 0;
			}
			double tmin = time2 - t0->val();
			// factor * int_tmin^tmax exp(-r(t0-t)) dt = factor / r * (exp(r tmax) - exp(r tmin));
			v += factor->val() / rate->val() * (exp(rate->val() * tmax) - exp(rate->val() * tmin));
		}

		return v / (time1 - time2);
	}

	void specialUpdate()	{
		// nothing to do
	}

	protected:

	Var<PosReal>* t0;
	Var<PosReal>* factor;
	Var<PosReal>* rate;
};

#endif

