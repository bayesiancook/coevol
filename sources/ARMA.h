
#ifndef ARMA_H
#define ARMA_H

class ARMAGammaWhiteNoiseBranchMean : public Rvar<PosReal>	{

	public:

	ARMAGammaWhiteNoiseBranchMean(Var<PosReal>* inalpha, Var<PosReal>* intime, Var<PosReal>* inupval, Var<PosReal>* inrho, Var<PosReal>* intime2)	{

		alpha = inalpha;
		time = intime;
		upval = inupval;
		rho = inrho;
		time2 = intime2;
		Register(alpha);
		Register(time);
		Register(upval);
		Register(rho);
		Register(time2);
		Sample();
	}

	double logProb()	{
		double a = alpha->val() * time->val();
		if (!upval)	{
			return a * log(a) - Random::logGamma(a) + (a-1) * log(this->val()) - a * this->val();
		}
	}

	protected:

	void drawSample()	{
		double a = alpha->val() * time->val();
		setval(Random::sGamma(a));
	}

	private:

	Var<PosReal>* alpha;
	Var<PosReal>* time;

};

class GammaWhiteNoiseProcess : public BranchProcess<PosReal>	{

	public:

	GammaWhiteNoiseProcess(LengthTree* intree, Var<PosReal>* inalpha)  : BranchProcess<PosReal>(intree->GetTree()) {
		SetWithRoot(false);
		alpha = inalpha;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~GammaWhiteNoiseProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return tree->GetTree();}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	protected:

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			this->GetBranchVal(link->GetBranch())->setval(1.0);
			RecursiveReset(link->Out());
		}
	}

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new GammaWhiteNoiseBranchMean(alpha,tree->GetBranchVal(link->GetBranch()));
	}

	private:

	Var<PosReal>* alpha;
	LengthTree* tree;

};

#endif

